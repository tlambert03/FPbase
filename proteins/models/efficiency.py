from django.contrib.contenttypes.fields import GenericForeignKey
from django.contrib.contenttypes.models import ContentType
from model_utils.models import TimeStampedModel
from django.db import models
from django.core.exceptions import ValidationError
from django.core.validators import MaxValueValidator, MinValueValidator

from ..util.efficiency import oc_efficiency_report


class OcFluorEff(TimeStampedModel):
    oc = models.ForeignKey('OpticalConfig', on_delete=models.CASCADE)
    limit = models.Q(app_label='proteins', model='state') | models.Q(app_label='proteins', model='dye')
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE,
                                     limit_choices_to=limit)
    object_id = models.PositiveIntegerField()
    fluor = GenericForeignKey('content_type', 'object_id')
    fluor_name = models.CharField(max_length=100, blank=True)
    ex_eff = models.FloatField(null=True, blank=True, verbose_name="Excitation Efficiency",
                               validators=[MinValueValidator(0), MaxValueValidator(1)])
    ex_eff_broad = models.FloatField(null=True, blank=True, verbose_name="Excitation Efficiency (Broadband)",
                                     validators=[MinValueValidator(0), MaxValueValidator(1)])
    em_eff = models.FloatField(null=True, blank=True, verbose_name="Emission Efficiency",
                               validators=[MinValueValidator(0), MaxValueValidator(1)])
    brightness = models.FloatField(null=True, blank=True)

    class Meta:
        unique_together = ("oc", "content_type", 'object_id')

    def clean(self):
        if (self.content_type_id not in
                ContentType.objects.filter(self.limit).values_list('id', flat=True)):
            mods = ContentType.objects.filter(self.limit).values_list('model', flat=True)
            raise ValidationError('ContentType for OcFluorEff.fluor must be in: {}'
                                  .format(",".join(mods)))

    def update_effs(self):
        rep = oc_efficiency_report(self.oc, [self.fluor]).get(self.fluor.slug, {})
        self.ex_eff = rep.get('ex')
        self.ex_eff_broad = rep.get('ex_broad')
        self.em_eff = rep.get('em')
        self.brightness = rep.get('bright')

    def save(self, *args, **kwargs):
        if (self.pk is None or
                (self.modified < self.oc.modified) or
                (self.modified < self.fluor.modified)):
            self.update_effs()
        if not self.fluor_name:
            self.fluor_name = str(self.fluor)
        super().save(*args, **kwargs)
