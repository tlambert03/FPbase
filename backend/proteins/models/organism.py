from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING
from urllib.parse import quote

from django.core.exceptions import ValidationError
from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel

from proteins.extrest import entrez
from proteins.models.mixins import Authorable

if TYPE_CHECKING:
    from proteins.models import Protein


class Organism(Authorable, TimeStampedModel):
    """A class for the parental organism (species) from which the protein has been engineered"""

    # Attributes
    id = models.PositiveIntegerField(
        primary_key=True, verbose_name="Taxonomy ID", help_text="NCBI Taxonomy ID"
    )  # genbank protein accession number
    scientific_name = models.CharField(max_length=128, blank=True)
    division = models.CharField(max_length=128, blank=True)
    common_name = models.CharField(max_length=128, blank=True)
    species = models.CharField(max_length=128, blank=True)
    genus = models.CharField(max_length=128, blank=True)
    rank = models.CharField(max_length=128, blank=True)

    if TYPE_CHECKING:
        proteins: models.QuerySet[Protein]

    def __str__(self):
        return self.scientific_name

    class Meta:
        verbose_name = "Organism"
        ordering = ["scientific_name"]

    def get_absolute_url(self):
        return reverse("proteins:organism-detail", args=[self.pk])

    def save(self, *args, **kwargs):
        if info := entrez.get_organism_info(self.id):
            self.__dict__.update(info)

        super().save(*args, **kwargs)

    def url(self):
        return self.get_absolute_url()

    @property
    def italicize(self) -> bool:
        """Whether the name is a genus/species name (italicized by convention)."""
        name = self.scientific_name
        return (
            bool(name)
            and name[0].isupper()
            and self.rank in ("", "genus", "species", "subspecies")
        )


COMMONS_UPLOAD = "https://upload.wikimedia.org/wikipedia/commons"


class OrganismPhoto(models.Model):
    """A photograph of an organism, from Wikimedia Commons or uploaded.

    Commons photos are copied into our own storage with `manage.py fetch_organism_photos`;
    Wikimedia rate-limits hotlinked images, so we only link to Commons until then.
    """

    organism = models.OneToOneField(Organism, on_delete=models.CASCADE, related_name="photo")
    commons_file = models.CharField(
        max_length=255,
        blank=True,
        help_text="Wikimedia Commons file name, without the 'File:' prefix. "
        "Leave blank to use an uploaded image instead.",
    )
    image = models.ImageField(
        upload_to="organism_photos/",
        blank=True,
        width_field="width",
        height_field="height",
        help_text="The image we serve. Filled from the Commons file by "
        "`manage.py fetch_organism_photos`, or upload one.",
    )
    thumbnail = models.ImageField(upload_to="organism_photos/thumbs/", blank=True, editable=False)
    width = models.PositiveIntegerField(null=True, blank=True, help_text="Image width (px)")
    height = models.PositiveIntegerField(null=True, blank=True, help_text="Image height (px)")
    author = models.CharField(max_length=255)
    license = models.CharField(max_length=64, help_text="e.g. 'CC BY-SA 4.0' or 'Public domain'")
    license_url = models.URLField(blank=True)
    source_url = models.URLField(blank=True, help_text="Page describing the original image.")
    pictured = models.CharField(
        max_length=128,
        blank=True,
        help_text="Scientific name of the pictured organism, if it is not exactly this taxon "
        "(e.g. a named species standing in for an 'sp.' record).",
    )
    pictured_note = models.CharField(
        max_length=64,
        blank=True,
        help_text="How `pictured` relates, e.g. 'same genus' or 'synonym'.",
    )
    alt = models.CharField(max_length=255, blank=True, help_text="Description for screen readers.")

    class Meta:
        verbose_name = "Organism photo"

    def __str__(self) -> str:
        return f"Photo of {self.organism}"

    def clean(self) -> None:
        if not (self.commons_file or self.image):
            raise ValidationError("Provide either a Wikimedia Commons file or an uploaded image.")
        if self.commons_file and not self.width:
            raise ValidationError({"width": "Give the original width of the Commons file."})

    @property
    def is_commons(self) -> bool:
        return bool(self.commons_file)

    @property
    def title(self) -> str:
        """Title of the work (CC 2.x/3.0 attribution asks for it); the Commons file name."""
        return self.commons_file.removeprefix("File:").rsplit(".", 1)[0].strip()

    @property
    def short_credit(self) -> str:
        return f"Photo: {self.author} · {self.license}"

    @property
    def credit_source_url(self) -> str:
        if self.source_url:
            return self.source_url
        if self.commons_file:
            return f"https://commons.wikimedia.org/wiki/File:{_commons_name(self.commons_file)}"
        return ""

    def commons_url(self, width: int | None = None) -> str:
        """URL of the Commons file, as a thumbnail no wider than `width` when possible.

        `width` must be one of Wikimedia's standard thumbnail sizes.
        """
        name = _commons_name(self.commons_file)
        digest = hashlib.md5(name.encode()).hexdigest()
        quoted = quote(name)
        if width is None or (self.width and width >= self.width):
            return f"{COMMONS_UPLOAD}/{digest[0]}/{digest[:2]}/{quoted}"
        return f"{COMMONS_UPLOAD}/thumb/{digest[0]}/{digest[:2]}/{quoted}/{width}px-{quoted}"

    def size_at(self, width: int) -> tuple[int, int] | None:
        """(width, height) of the image when served at most `width` wide."""
        if not (self.width and self.height):
            return None
        w = min(width, self.width)
        return w, round(self.height * w / self.width)

    # template conveniences
    @property
    def src(self) -> str:
        if self.image:
            return self.image.url
        return self.commons_url(960) if self.commons_file else ""

    @property
    def src_size(self) -> tuple[int, int] | None:
        return self.size_at(960)

    @property
    def thumb_src(self) -> str:
        if self.thumbnail:
            return self.thumbnail.url
        if self.image:
            return self.image.url
        return self.commons_url(330) if self.commons_file else ""

    @property
    def thumb_size(self) -> tuple[int, int] | None:
        return self.size_at(330)

    @property
    def is_portrait(self) -> bool:
        return bool(self.width and self.height and self.height > self.width)

    @property
    def alt_text(self) -> str:
        return self.alt or f"Photograph of {self.pictured or self.organism.scientific_name}"


def _commons_name(file: str) -> str:
    return file.removeprefix("File:").strip().replace(" ", "_")
