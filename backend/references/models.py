from datetime import datetime

import pytz
from django.contrib.auth import get_user_model
from django.core.exceptions import ObjectDoesNotExist
from django.core.validators import (
    MaxLengthValidator,
    MaxValueValidator,
    MinLengthValidator,
    MinValueValidator,
)
from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel
from proteins.validators import validate_doi

from .helpers import doi_lookup, name_to_initials

User = get_user_model()


class Author(TimeStampedModel):
    family = models.CharField(max_length=80)
    given = models.CharField(max_length=80)
    initials = models.CharField(max_length=10)
    publications = models.ManyToManyField("Reference", through="ReferenceAuthor")

    @property
    def protein_contributions(self):
        return {p for ref in self.publications.all() for p in ref.primary_proteins.all()}

    @property
    def first_authorships(self):
        return [p.reference for p in self.referenceauthor_set.filter(author_idx=0)]

    @property
    def last_authorships(self):
        return [p.reference for p in self.referenceauthor_set.all() if p.author_idx == p.author_count - 1]

    def save(self, *args, **kwargs):
        self.initials = name_to_initials(self.initials)
        super().save(*args, **kwargs)

    def get_absolute_url(self):
        return reverse("references:author-detail", args=[self.id])

    def full_name(self):
        if self.given:
            return f"{self.given} {self.family}".title()
        else:
            return f"{self.initials} {self.family}".title()

    def __repr__(self):
        return f"Author(family='{self.family}', given='{self.given}'')"

    def __str__(self):
        return f"{self.family} {self.initials}"

    class Meta:
        unique_together = (("family", "initials"),)


class Reference(TimeStampedModel):
    doi = models.CharField(
        max_length=50,
        unique=True,
        blank=False,
        verbose_name="DOI",
        validators=[validate_doi],
    )
    pmid = models.CharField(max_length=15, unique=True, null=True, blank=True, verbose_name="Pubmed ID")
    title = models.CharField(max_length=512, blank=True)
    journal = models.CharField(max_length=512, blank=True)
    pages = models.CharField(max_length=20, blank=True)
    volume = models.CharField(max_length=10, blank=True, default="")
    issue = models.CharField(max_length=10, blank=True, default="")
    firstauthor = models.CharField(max_length=100, blank=True, default="")
    citation = models.CharField(max_length=256, blank=True, default="")
    date = models.DateField(null=True, blank=True)
    year = models.PositiveIntegerField(
        validators=[
            MinLengthValidator(4),
            MaxLengthValidator(4),
            MinValueValidator(1960),
            MaxValueValidator(datetime.now(pytz.utc).year + 1),
        ],
        help_text="YYYY",
    )
    authors = models.ManyToManyField("Author", through="ReferenceAuthor")
    summary = models.CharField(max_length=512, blank=True, help_text="Brief summary of findings")

    created_by = models.ForeignKey(
        User,
        related_name="reference_author",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )
    updated_by = models.ForeignKey(
        User,
        related_name="reference_modifier",
        blank=True,
        null=True,
        on_delete=models.CASCADE,
    )

    class Meta:
        ordering = ("date",)

    def get_authors(self):
        return self.authors.order_by("referenceauthor")

    @property
    def first_author(self):
        try:
            return self.get_authors()[0]
        except Exception:
            return None

    @property
    def protein_secondary_reference(self):
        return self.proteins.exclude(id__in=self.primary_proteins.all())

    def get_citation(self, authorlist):
        try:
            if len(authorlist) == 0:
                if self.title:
                    return f"{self.title[:30]}..."
                else:
                    return f"{self.doi}"
            elif len(authorlist) > 2:
                middle = " et al. "
            elif len(authorlist) == 2:
                secondauthor = authorlist[1].family
                middle = f" & {secondauthor} "
            else:
                middle = " "
            self.firstauthor = authorlist[0].family
            return f"{self.firstauthor}{middle}({self.year})"
        except ObjectDoesNotExist:
            return f"doi: {self.doi}"

    def get_absolute_url(self):
        return reverse("references:reference-detail", args=[self.id])

    def __repr__(self):
        return f"Reference(doi={self.doi})"

    def __str__(self):
        try:
            return self.citation
        except Exception:
            return super().__str__()

    def clean(self):
        if self.doi:
            self.doi = self.doi.strip()

    def save(self, skipdoi=False, *args, **kwargs):
        if self.doi:
            self.doi = self.doi.lower()
        if not skipdoi:
            info = doi_lookup(self.doi)
            authors = info.pop("authors")
            authorlist = []
            for author in authors:
                auth, _ = Author.objects.get_or_create(
                    initials=name_to_initials(author["given"]),
                    family=author["family"],
                    defaults={"given": author["given"].replace(".", "")},
                )
                authorlist.append(auth)
            for k, v in info.items():
                setattr(self, k, v)
            self.citation = self.get_citation(authorlist)
        super().save(*args, **kwargs)
        if not skipdoi:
            ReferenceAuthor.objects.filter(reference_id=self.id).delete()
            for idx, author in enumerate(authorlist):
                authmemb = ReferenceAuthor(reference=self, author=author, author_idx=idx)
                authmemb.save()

    def prot_secondary(self):
        return [p.name for p in self.protein_secondary_reference]

    def prot_primary(self):
        return [p.name for p in self.primary_proteins.all()]

    def _excerpts(self):
        return [exc.content for exc in self.excerpts.all()]

    def url(self):
        return self.get_absolute_url()


class ReferenceAuthor(models.Model):
    reference = models.ForeignKey(Reference, on_delete=models.CASCADE)
    author = models.ForeignKey(Author, on_delete=models.CASCADE)
    author_idx = models.PositiveSmallIntegerField()

    class Meta:
        ordering = ["author_idx"]

    def __str__(self):
        return f"<AuthorMembership: {self.author} in doi: {self.reference.doi}>"

    @property
    def author_count(self):
        return self.reference.authors.count()
