import django_filters
from django_filters import rest_framework as filters
from django import forms
from .models import Protein, State
from Bio import Seq, Alphabet
from .validators import cdna_sequence_validator


class StateFilter(filters.FilterSet):
    ex_spectra = django_filters.BooleanFilter(name='ex_spectra', lookup_expr='isnull')
    em_spectra = django_filters.BooleanFilter(name='em_spectra', lookup_expr='isnull')

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

    def get_specbright(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.local_brightness == value]

    def get_specbright_lt(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.local_brightness < value]

    def get_specbright_gt(self, queryset, name, value):
        qsALL = list(queryset.all())
        return [P for P in qsALL if P.local_brightness > value]

    class Meta:
        model = State
        order_by = 'em_max'
        fields = {
            'name': ['icontains', 'iendswith', 'istartswith', 'iexact', ],
            'ex_max': ['around', 'range', 'lte', 'gte', 'exact'],
            'em_max': ['around', 'range', 'lte', 'gte', 'exact'],
            'lifetime': ['gte',  'lte', 'range', 'exact'],
            'maturation': ['gte',  'lte', 'range', 'exact'],
            'ext_coeff': ['gte',  'lte', 'range', 'exact'],
            'qy': ['gte',  'lte', 'range', 'exact'],
            'brightness': ['gte',  'lte', 'range', 'exact'],
            'pka': ['gte',  'lte', 'range', 'exact'],
            'bleach_measurements__rate': ['gte',  'lte', 'range', 'exact'],
            'spectral_brightness': ['gt', 'lt'],
            'ex_spectra': ['isnull'],
            'em_spectra': ['isnull'],
        }


class ProteinFilterForm(forms.Form):
    def clean_seq__cdna_contains(self):
        seq = self.cleaned_data['seq__cdna_contains'].replace(' ', '').replace('\n', '').upper()
        cdna_sequence_validator(seq)
        return seq

    def clean_seq__icontains(self):
        seq = self.cleaned_data['seq__icontains'].replace(' ', '').replace('\n', '').upper()
        return seq


class CharArrayFilter(filters.BaseCSVFilter, filters.CharFilter):
    pass


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
    pdb__contains = CharArrayFilter(name='pdb', lookup_expr='contains')
    # aliases__contains = CharArrayFilter(name='aliases', lookup_expr='icontains')

    class Meta:
        model = Protein
        form = ProteinFilterForm
        order_by = 'default_state__em_max'
        fields = {
            'name': ['icontains', 'iendswith', 'istartswith', 'iexact', ],
            # 'aliases': ['contains'],
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
            'genbank': ['iexact'],
            'pdb': ['contains'],
            'uniprot': ['iexact'],
            'status': ['exact'],
            'switch_type': ['exact', 'ne'],
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
            'ne': 'is not',
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
            'uniprot': 'UniProtKB ID',
            'genbank': 'GenBank ID',
            'pdb': 'PDB ID',
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
        print(qsALL)
        return [P for P in qsALL if P.default_state and P.default_state.local_brightness > value]

    def translate_cdna(self, queryset, name, value):
        coding_dna = Seq.Seq(value, Alphabet.IUPAC.unambiguous_dna)
        return queryset.filter(seq__icontains=coding_dna.translate())

