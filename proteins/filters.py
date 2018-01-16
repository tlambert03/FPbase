import django_filters
from django_filters import rest_framework as filters
from django import forms
from .models import Protein
from Bio import Seq, Alphabet
from .validators import cdna_sequence_validator


class ProteinFilterForm(forms.Form):
    def clean_seq__cdna_contains(self):
        seq = self.cleaned_data['seq__cdna_contains'].replace(' ', '').replace('\n', '').upper()
        cdna_sequence_validator(seq)
        return seq

    def clean_seq__icontains(self):
        seq = self.cleaned_data['seq__icontains'].replace(' ', '').replace('\n', '').upper()
        return seq


class ProteinFilter(filters.FilterSet):
    spectral_brightness = django_filters.NumberFilter(
        name='spectral_brightness',
        method='get_specbright',
        help_text='fold brightness relative to spectral neighbors')
    spectral_brightness__gt = django_filters.NumberFilter(
        name='spectral_brightness',
        method='get_specbright_gt', lookup_expr='gt',
        help_text='fold brightness relative to spectral neighbors')
    spectral_brightness__lt = django_filters.NumberFilter(
        name='spectral_brightness',
        method='get_specbright_lt', lookup_expr='lt',
        help_text='fold brightness relative to spectral neighbors')
    seq__cdna_contains = django_filters.CharFilter(
        name='seq',
        method='translate_cdna', lookup_expr='cdna_contains',
        help_text='cDNA sequence (in frame)')

    class Meta:
        model = Protein
        form = ProteinFilterForm
        order_by = 'default_state__em_max'
        fields = {
            'name': ['icontains', 'iendswith', 'istartswith', 'iexact', ],
            'seq': ['icontains', 'iendswith', 'istartswith', 'cdna_contains'],
            'default_state__ex_max': ['around', 'range', 'lte', 'gte', 'exact'],
            'default_state__em_max': ['around', 'range', 'lte', 'gte', 'exact'],
            'default_state__lifetime': ['gte',  'lte', 'range', 'exact'],
            'default_state__maturation': ['gte',  'lte', 'range', 'exact'],
            'default_state__ext_coeff': ['gte',  'lte', 'range', 'exact'],
            'default_state__qy': ['gte',  'lte', 'range', 'exact'],
            'default_state__brightness': ['gte',  'lte', 'range', 'exact'],
            'default_state__pka': ['gte',  'lte', 'range', 'exact'],
            'default_state__bleach_measurements__rate': ['gte',  'lte', 'range', 'exact'],
            'agg': ['exact'],
            #'status': ['exact'],
            'switch_type': ['exact'],
            'parent_organism': ['exact'],
            'primary_reference__year': ['gte', 'gt', 'lt', 'lte', 'range', 'exact'],
            'spectral_brightness': ['gt', 'lt'],
        }
        operators = {
            'lt': 'is less than',
            'gt': 'is greater than',
            'lte': 'is less than or equal to',
            'gte': 'is greater than or equal to',
            'around': 'is around',
            'exact': 'is',
            'iexact': 'is',
            'range': 'is between',
            'contains': 'contains (case sensitive)',
            'icontains': 'contains',
            'iendswith': 'ends with',
            'istartswith': 'starts with',
            'cdna_contains': 'cDNA could contain'
        }
        labels = {
            'default_state__ex_max': 'Excitation Maximum',
            'default_state__em_max': 'Emission Maximum',
            'default_state__lifetime': 'Lifetime (ns)',
            'default_state__maturation': 'Maturation (min)',
            'default_state__ext_coeff': 'Extinction Coefficient',
            'default_state__qy': 'Quantum Yield',
            'default_state__brightness': 'Brightness',
            'default_state__pka': 'pKa',
            'seq': 'Sequence',
            'agg': 'Oligomerization',
            'primary_reference__year': 'Year published',
            'default_state__bleach_measurements__rate': 'Photostability (s)',
        }

    def get_specbright(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.default_state and P.default_state.local_brightness == value]

    def get_specbright_lt(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.default_state and P.default_state.local_brightness < value]

    def get_specbright_gt(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.default_state and P.default_state.local_brightness > value]

    def translate_cdna(self, queryset, name, value):
        coding_dna = Seq.Seq(value, Alphabet.IUPAC.unambiguous_dna)
        return queryset.filter(seq__icontains=coding_dna.translate())

