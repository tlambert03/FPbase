# -*- coding: utf-8 -*-
from datetime import datetime
from django.db import models
from django.contrib.auth import get_user_model
from django.core.exceptions import ObjectDoesNotExist
from django.core.validators import MinValueValidator, MaxValueValidator
from model_utils.models import TimeStampedModel
from .helpers import doi_lookup, name_to_initials
from django.urls import reverse

User = get_user_model()


class Author(models.Model):
    family = models.CharField(max_length=80)
    given = models.CharField(max_length=80)
    initials = models.CharField(max_length=10)
    publications = models.ManyToManyField('Reference', through='ReferenceAuthor')

    @property
    def protein_contributions(self):
        return set([p for ref in self.publications.all() for p in ref.primary_proteins.all()])

    @property
    def first_authorships(self):
        return [p.reference for p in self.referenceauthor_set.filter(author_idx=0)]

    @property
    def last_authorships(self):
        return [p.reference for p in self.referenceauthor_set.all() if p.author_idx == p.author_count-1]

    def save(self, *args, **kwargs):
        self.initials = name_to_initials(self.initials)
        super().save(*args, **kwargs)

    def get_absolute_url(self):
        return reverse("references:author-detail", args=[self.id])

    def full_name(self):
        if self.given:
            return "{} {}".format(self.given, self.family)
        else:
            return "{} {}".format(self.initials, self.family)

    def __repr__(self):
        return "Author(family='{}', initials='{}')".format(self.family, self.initials)

    def __str__(self):
        return self.family + ' ' + self.initials

    class Meta:
        unique_together = (('family', 'initials'),)


class Reference(TimeStampedModel):
    doi = models.CharField(max_length=50, unique=True, blank=False, verbose_name="DOI")
    pmid = models.CharField(max_length=15, unique=True, null=True, blank=True, verbose_name="Pubmed ID")
    title = models.CharField(max_length=512, blank=True)
    journal = models.CharField(max_length=512, blank=True)
    pages = models.CharField(max_length=20, blank=True)
    volume = models.CharField(max_length=10, blank=True, default='')
    issue = models.CharField(max_length=10, blank=True, default='')
    year = models.PositiveIntegerField(
            validators=[
                MinValueValidator(1960),
                MaxValueValidator(datetime.now().year)],
            help_text="Use the following format: <YYYY>")
    authors = models.ManyToManyField("Author", through='ReferenceAuthor')

    created_by = models.ForeignKey(User, related_name='reference_author', blank=True, null=True)
    updated_by = models.ForeignKey(User, related_name='reference_modifier', blank=True, null=True)

    @property
    def first_author(self):
        try:
            return self.referenceauthor_set.get(author_idx=0).author
        except Exception:
            return None

    @property
    def ordered_authors(self):
        return [ra.author for ra in self.referenceauthor_set.all()]

    @property
    def protein_secondary_reference(self):
        return self.proteins.exclude(id__in=self.primary_proteins.all())

    def clean(self):
        if self.doi:
            self.doi = self.doi.strip()

    def save(self, *args, **kwargs):
        info = doi_lookup(self.doi)
        authors = info.pop('authors')
        authorlist = []
        for author in authors:
            auth, _ = Author.objects.get_or_create(
                initials=name_to_initials(author['given']),
                family=author['family'],
                defaults={'given': author['given'].replace('.', '')}
            )
            authorlist.append(auth)
        for k, v in info.items():
            setattr(self, k, v)
        super().save(*args, **kwargs)
        for idx, author in enumerate(authorlist):
            authmemb = ReferenceAuthor(
                reference=self,
                author=author,
                author_idx=idx,
            )
            authmemb.save()

    @property
    def citation(self):
        try:
            if self.authors.count() == 0:
                if self.title:
                    return "{}...".format(self.title[:30])
                else:
                    return "{}".format(self.doi)
            elif self.authors.count() > 2:
                middle = ' et al. '
            elif self.authors.count() == 2:
                secondauthor = self.referenceauthor_set.get(author_idx=1).author.family
                middle = ' & {} '.format(secondauthor)
            else:
                middle = ' '

            return "{}{}({})".format(self.first_author.family, middle, self.year)
        except ObjectDoesNotExist:
                return "doi: {}".format(self.doi)

    def get_absolute_url(self):
        return reverse("references:reference-detail", args=[self.id])

    def __repr__(self):
        return "Reference(doi={})".format(self.doi)

    def __str__(self):
        try:
            return self.citation
        except Exception:
            return super(Reference, self).__str__()


class ReferenceAuthor(models.Model):
    reference = models.ForeignKey(Reference, on_delete=models.CASCADE)
    author = models.ForeignKey(Author, on_delete=models.CASCADE)
    author_idx = models.PositiveSmallIntegerField()

    @property
    def author_count(self):
        return self.reference.authors.count()

    def __str__(self):
        return "<AuthorMembership: {} in doi: {}>".format(self.author, self.reference.doi)

    class Meta:
        ordering = ['author_idx']
