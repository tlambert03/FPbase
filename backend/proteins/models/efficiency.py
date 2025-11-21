from typing import TYPE_CHECKING

from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.db.models import F, Max, OuterRef, Q, Subquery
from model_utils.models import TimeStampedModel

from proteins.models.fluorophore import Fluorophore
from proteins.util.efficiency import oc_efficiency_report

if TYPE_CHECKING:
    from proteins.models import OpticalConfig  # noqa F401


class OcFluorEffQuerySet(models.QuerySet):
    def outdated(self):
        fluor_objs = Fluorophore.objects.filter(id=OuterRef("object_id"))
        spectra_mod = fluor_objs.annotate(latest_spec=Max("spectra__modified")).values("latest_spec")[:1]

        fluor_mod = fluor_objs.values("modified")[:1]
        return self.annotate(
            fluor_mod=Subquery(fluor_mod),
            spec_mod=Subquery(spectra_mod),
        ).filter(
            Q(modified__lt=F("fluor_mod")) | Q(modified__lt=F("spec_mod")) | Q(modified__lt=F("oc__modified")),
        )


class OcFluorEff(TimeStampedModel):
    oc = models.ForeignKey["OpticalConfig"]("OpticalConfig", on_delete=models.CASCADE)
    limit = models.Q(app_label="proteins", model__in=("state", "dyestate"))
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE, limit_choices_to=limit)
    object_id = models.PositiveIntegerField()
    fluor = GenericForeignKey("content_type", "object_id")
    fluor_name = models.CharField(max_length=100, blank=True)
    ex_eff = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Excitation Efficiency",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )
    ex_eff_broad = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Excitation Efficiency (Broadband)",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )
    em_eff = models.FloatField(
        null=True,
        blank=True,
        verbose_name="Emission Efficiency",
        validators=[MinValueValidator(0), MaxValueValidator(1)],
    )
    brightness = models.FloatField(null=True, blank=True)
    objects = OcFluorEffQuerySet.as_manager()

    class Meta:
        unique_together = ("oc", "content_type", "object_id")

    def clean(self):
        allowed = list(ContentType.objects.filter(self.limit).values_list("id", "model"))
        allowed_ids = {item[0] for item in allowed}
        if self.content_type_id not in allowed_ids:
            models_list = ", ".join(sorted({item[1] for item in allowed}))
            raise ValidationError(f"ContentType for OcFluorEff.fluor must be in: {models_list}")

    def update_effs(self):
        rep = oc_efficiency_report(self.oc, [self.fluor]).get(self.fluor.slug, {})
        self.ex_eff = rep.get("ex")
        self.ex_eff_broad = rep.get("ex_broad")
        self.em_eff = rep.get("em")
        self.brightness = rep.get("bright")

    @property
    def outdated(self):
        oc_modified = getattr(self.oc, "modified", None)
        fluor_modified = getattr(self.fluor, "modified", None)
        spect_modified = None
        if hasattr(self.fluor, "spectra"):
            spect_modified = self.fluor.spectra.aggregate(latest=Max("modified")).get("latest")

        for candidate in (oc_modified, fluor_modified, spect_modified):
            if candidate and self.modified < candidate:
                return True
        return False

    def save(self, *args, **kwargs):
        if self.pk is None or self.outdated:
            self.update_effs()
        if not self.fluor_name:
            self.fluor_name = str(self.fluor)
        super().save(*args, **kwargs)

    def __repr__(self):
        return f"<OcFluorEff: {self.oc.name} with {self.fluor.slug}>"
