# -*- coding: utf-8 -*-
import json
from datetime import datetime
from django.db import models
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError, ObjectDoesNotExist
from metapub import pubmedcentral, CrossRef, PubMedFetcher
from django.conf import settings

#TODO: remove this alltogether
from Bio import Entrez
Entrez.email = "talley_lambert@hms.harvard.edu"

User = settings.AUTH_USER_MODEL


def get_pmid_from_doi(doi):
    CR = CrossRef()
    PMF = PubMedFetcher()

    pubmedID = pubmedcentral.get_pmid_for_otherid(doi)
    if not pubmedID:
        # wasn't found based on DOI
        # try crossref:
        results = CR.query(doi)
        if len(results) == 1:
            ids = PMF.pmids_for_query(results[0]['title'])
            if len(ids) == 1:
                pubmedID = ids[0]
            else:
                ids = PMF.pmids_for_citation(**results[0]['slugs'])
                if len(ids) == 1 and not ids[0] == 'NOT_FOUND':
                    pubmedID = ids[0]
    return pubmedID


class Author(models.Model):
    initials = models.CharField(max_length=4)
    last_name = models.CharField(max_length=80)
    publications = models.ManyToManyField('Reference', through='ReferenceAuthor')

    @property
    def protein_contributions(self):
        return set([p for ref in self.publications.all() for p in ref.protein_primary_reference.all()])

    @property
    def first_authorships(self):
        return [p.reference for p in self.referenceauthor_set.filter(author_idx=0)]

    @property
    def last_authorships(self):
        return [p.reference for p in self.referenceauthor_set.filter(author_idx=0)]

    def __str__(self):
        return self.last_name + ' ' + self.initials


class Reference(models.Model):
    pmid = models.CharField(max_length=15, unique=True, blank=True, null=True, verbose_name="PMID")
    doi = models.CharField(max_length=50, unique=True, blank=True, null=True, verbose_name="DOI")
    added_by = models.ForeignKey(User, related_name='entries', blank=True, null=True)
    updated_by = models.ForeignKey(User, related_name='entries_modifiers', blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    title = models.CharField(max_length=512, blank=True)
    journal = models.CharField(max_length=512, blank=True)
    pages = models.CharField(max_length=20, blank=True)
    volume = models.CharField(max_length=10, blank=True, null=True)
    #_author_list = models.CharField(max_length=300, blank=True, db_column="author_list")
    authors = models.ManyToManyField("Author", through='ReferenceAuthor')
    #pubdate = models.CharField(max_length=64, blank=True)
    pubdate = models.DateTimeField(blank=True, null=True)
    so = models.CharField(max_length=128, blank=True, verbose_name=u'SO')
    ref = models.CharField(max_length=512, blank=True, null=True)

    @property
    def first_author(self):
        return self.referenceauthor_set.get(author_idx=0).author

    @classmethod
    def create(cls, pmid=None, doi=None, **kwargs):
        """creation method that allows passing either pmid or doi and the other will be looked up"""
        if not (pmid or doi):
            raise ValueError("Must provide at least pmid or doi when creating Reference")
        if doi and not pmid:
            pmid = get_pmid_from_doi(doi)
        if pmid and not doi:
            doi = pubmedcentral.get_doi_for_otherid(pmid)
        ref = cls(pmid=pmid, doi=doi, **kwargs)
        return ref

    def populate_from_pubmed(self):
        if self.pmid:
            pubmed_record = Entrez.read(Entrez.esummary(db='pubmed', id=self.pmid, retmode='xml'))
            if len(pubmed_record):
                pubmed_record = pubmed_record[0]
                if not self.doi:
                    try:
                        self.doi = pubmed_record['DOI']
                    except AttributeError:
                        self.doi = None
                self.title = pubmed_record['Title']
                self.journal = pubmed_record['Source']
                self.pages = pubmed_record['Pages']
                self.volume = pubmed_record['Volume']
                # FIXME: this is going to break
                try:
                    self.pubdate = datetime.strptime(pubmed_record['PubDate'], '%Y %b %d')
                except ValueError:
                    try:
                        self.pubdate = datetime.strptime(pubmed_record['PubDate'], '%Y %b')
                    except ValueError:
                        self.pubdate = datetime.strptime(pubmed_record['PubDate'].split(' ')[0], '%Y')

                self.so = pubmed_record['SO']
                try:
                    self.ref = pubmed_record['References']
                except IndexError:
                    self.ref = None
                super(Reference, self).save()
                for idx, auth in enumerate(pubmed_record['AuthorList']):
                    authsplit = auth.split(' ')
                    last_name = authsplit[0] if len(authsplit) else None
                    initials = authsplit[1] if len(authsplit) > 1 else None

                    author, created = Author.objects.get_or_create(
                        initials=initials,
                        last_name=last_name,
                    )
                    authmemb = ReferenceAuthor(
                        reference=self,
                        author=author,
                        author_idx=idx,
                    )
                    authmemb.save()
        elif self.doi:
            CR = CrossRef()
            query = CR.query(self.doi)
            if len(query) == 1:
                query = query[0]
                self.title = query.get('title')
                self.journal = query['slugs'].get('jtitle')
                self.pages = "-".join([query['slugs'].get('spage'), query['slugs'].get('epage')])
                self.volume = query['slugs'].get('volume')
                self.pubdate = datetime.strptime(query.get('year'), '%Y')
            else:
                return

    def clean(self):
        super(Reference, self).clean()
        if self.pmid is None and self.doi is None:
            raise ValidationError('Reference must have either pmid or doi')

    def save(self, *args, **kwargs):
        self.populate_from_pubmed()
        super(Reference, self).save(*args, **kwargs)

    @property
    def citation(self):
        try:
            return "{} ({})".format(self.first_author.last_name, self.pubdate.strftime('%Y'))
        except ObjectDoesNotExist:
            if self.pmid:
                return "PMID: {}".format(self.pmid)
            elif self.doi:
                return "DOI: {}".format(self.doi)

    def __repr__(self):
        return "Reference(pmid={})".format(self.pmid)

    def __str__(self):
        try:
            return self.citation
        except Exception:
            return super(Reference, self).__str__()


class ReferenceAuthor(models.Model):
    reference = models.ForeignKey(Reference, on_delete=models.CASCADE)
    author = models.ForeignKey(Author, on_delete=models.CASCADE)
    author_idx = models.PositiveSmallIntegerField()

    def __str__(self):
        return "<AuthorMembership: {} in PMID: {}>".format(self.author, self.reference.pmid)

    class Meta:
        ordering = ['author_idx']