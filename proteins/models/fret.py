# -*- coding: utf-8 -*-
from django.db import models
from model_utils.models import TimeStampedModel
from references.models import Reference
from .mixins import Authorable


class FRETpair(Authorable, TimeStampedModel):
    # relational class for FRET pairs to hold attributes about the pair

    # Attributes
    radius   = models.FloatField(blank=True, null=True)

    # Relations
    donor    = models.ForeignKey('Protein', null=False, blank=False, verbose_name='donor', related_name='FK_FRETdonor_protein', on_delete=models.CASCADE)
    acceptor = models.ForeignKey('Protein', null=False, blank=False, verbose_name='acceptor', related_name='FK_FRETacceptor_protein', on_delete=models.CASCADE)

    pair_references = models.ManyToManyField(Reference, related_name='FK_FRETpair_reference', blank=True)  # any additional papers that reference the FRET pair

    @property
    def name(self):
        return self.donor.name + '-' + self.acceptor.name

    @property
    def spectral_overlap(self):
        accEx  = self.acceptor.default_state.ex_spectra
        accEC  = self.acceptor.default_state.ext_coeff
        donEm  = self.donor.default_state.em_spectra
        # donQY  = self.donor.default_state.qy
        donCum = sum(donEm.y)
        minAcc = accEx.min_wave
        maxAcc = accEx.max_wave
        minEm  = donEm.min_wave
        maxEm  = donEm.max_wave

        startingwave = int(max(minAcc, minEm))
        endingwave = int(min(maxAcc, maxEm))

        A = accEx.wave_value_pairs()
        D = donEm.wave_value_pairs()
        overlap = [(pow(wave, 4) * A[wave] * accEC * D[wave] / donCum) for wave in range(startingwave, endingwave + 1)]

        return sum(overlap)

    def forsterDist(self, n=1.33, k=2/3):
        return .2108 * (pow((k) * (pow(n, -4) * self.spectral_overlap), (1. / 6.)))

    def __str__(self):
        return self.name

    class Meta:
        verbose_name = 'FRET Pair'
        unique_together = (('donor', 'acceptor'),)
