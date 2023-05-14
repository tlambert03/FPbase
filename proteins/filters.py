import django_filters
from Bio import Seq
from django import forms
from django_filters import rest_framework as filters

from .models import Organism, Protein, Spectrum, State
from .validators import cdna_sequence_validator


class SpectrumFilter(filters.FilterSet):
    class Meta:
        model = Spectrum
        fields = ("category", "subtype", "id", "owner_state")


class StateFilter(filters.FilterSet):
    ex_spectra = django_filters.BooleanFilter(field_name="ex_spectra", lookup_expr="isnull")
    em_spectra = django_filters.BooleanFilter(field_name="em_spectra", lookup_expr="isnull")
    spectral_brightness = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright",
        help_text="fold brightness relative to spectral neighbors",
    )
    spectral_brightness__gt = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright_gt",
        lookup_expr="gt",
        help_text="fold brightness relative to spectral neighbors",
    )
    spectral_brightness__lt = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright_lt",
        lookup_expr="lt",
        help_text="fold brightness relative to spectral neighbors",
    )

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
        order_by = "em_max"
        fields = {
            "name": ["icontains", "iendswith", "istartswith", "iexact"],
            "ex_max": ["around", "range", "lte", "gte", "exact"],
            "em_max": ["around", "range", "lte", "gte", "exact"],
            "lifetime": ["gte", "lte", "range", "exact"],
            "maturation": ["gte", "lte", "range", "exact"],
            "ext_coeff": ["gte", "lte", "range", "exact"],
            "qy": ["gte", "lte", "range", "exact"],
            "brightness": ["gte", "lte", "range", "exact"],
            "pka": ["gte", "lte", "range", "exact"],
            "bleach_measurements__rate": ["gte", "lte", "range", "exact"],
        }


class ProteinFilterForm(forms.Form):
    def clean_seq__cdna_contains(self):
        seq = self.cleaned_data["seq__cdna_contains"].replace(" ", "").replace("\n", "").upper()
        cdna_sequence_validator(seq)
        return seq

    def clean_seq__icontains(self):
        seq = self.cleaned_data["seq__icontains"].replace(" ", "").replace("\n", "").upper()
        return seq


class CharArrayFilter(filters.BaseCSVFilter, filters.CharFilter):
    pass


class ProteinFilter(filters.FilterSet):
    spectral_brightness = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright",
        help_text="fold brightness relative to spectral neighbors",
    )
    spectral_brightness__gt = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright_gt",
        lookup_expr="gt",
        help_text="fold brightness relative to spectral neighbors",
    )
    spectral_brightness__lt = django_filters.NumberFilter(
        field_name="spectral_brightness",
        method="get_specbright_lt",
        lookup_expr="lt",
        help_text="fold brightness relative to spectral neighbors",
    )
    seq__cdna_contains = django_filters.CharFilter(
        field_name="seq",
        method="translate_cdna",
        lookup_expr="cdna_contains",
        help_text="cDNA sequence (in frame)",
    )
    pdb__contains = CharArrayFilter(field_name="pdb", lookup_expr="contains")
    # aliases__contains = CharArrayFilter(field_name='aliases', lookup_expr='icontains')
    name__icontains = django_filters.CharFilter(
        field_name="name", method="name_or_alias_icontains", lookup_expr="icontains"
    )
    switch_type__ne = django_filters.ChoiceFilter(choices=Protein.SWITCHING_CHOICES, method="switch_type__notequal")
    cofactor__ne = django_filters.ChoiceFilter(choices=Protein.COFACTOR_CHOICES, method="cofactor__notequal")
    parent_organism__ne = django_filters.ModelChoiceFilter(
        queryset=Organism.objects.all(), method="parent_organism__notequal"
    )

    # name__iexact = django_filters.CharFilter(
    #     field_name='name',
    #     method='name_or_alias_iexact', lookup_expr='iexact')
    # name__iendswith = django_filters.CharFilter(
    #     field_name='name',
    #     method='name_or_alias_iendswith', lookup_expr='iendswith')
    # name__istartswith = django_filters.CharFilter(
    #     field_name='name',
    #     method='name_or_alias_istartswith', lookup_expr='istartswith')

    class Meta:
        model = Protein
        form = ProteinFilterForm
        order_by = "default_state__em_max"
        fields = {
            "name": ["icontains", "iendswith", "istartswith", "iexact"],
            # 'aliases': ['contains'],
            "seq": ["icontains", "iendswith", "istartswith", "cdna_contains"],
            "default_state__ex_max": ["around", "range", "lte", "gte", "exact"],
            "default_state__em_max": ["around", "range", "lte", "gte", "exact"],
            "default_state__lifetime": ["gte", "lte", "range", "exact"],
            "default_state__maturation": ["gte", "lte", "range", "exact"],
            "default_state__ext_coeff": ["gte", "lte", "range", "exact"],
            "default_state__qy": ["gte", "lte", "range", "exact"],
            "default_state__brightness": ["gte", "lte", "range", "exact"],
            "default_state__pka": ["gte", "lte", "range", "exact"],
            "default_state__bleach_measurements__rate": [
                "gte",
                "lte",
                "range",
                "exact",
            ],
            "agg": ["exact"],
            "uuid": ["iexact"],
            "genbank": ["iexact"],
            "pdb": ["contains"],
            "uniprot": ["iexact"],
            "status": ["exact"],
            "switch_type": ["exact", "ne"],
            "cofactor": ["exact", "ne"],
            "parent_organism": ["exact", "ne"],
            "primary_reference__year": ["gte", "gt", "lt", "lte", "range", "exact"],
            "primary_reference__author__family": ["icontains"],
            "slug": ["exact"],
            "id": ["exact"],
        }
        form_fields = dict(**fields, spectral_brightness=["gt", "lt"])
        operators = {
            "lt": "is less than",
            "gt": "is greater than",
            "lte": "is less than or equal to",
            "gte": "is greater than or equal to",
            "around": "is around",
            "exact": "is",
            "ne": "is not",
            "iexact": "is",
            "range": "is between",
            "contains": "contains (case sensitive)",
            "icontains": "contains",
            "iendswith": "ends with",
            "istartswith": "starts with",
            "cdna_contains": "cDNA could contain",
        }
        labels = {
            "default_state__ex_max": "Excitation Maximum",
            "default_state__em_max": "Emission Maximum",
            "default_state__lifetime": "Lifetime (ns)",
            "default_state__maturation": "Maturation (min)",
            "default_state__ext_coeff": "Extinction Coefficient",
            "default_state__qy": "Quantum Yield",
            "default_state__brightness": "Brightness",
            "default_state__pka": "pKa",
            "uniprot": "UniProtKB ID",
            "genbank": "GenBank ID",
            "pdb": "PDB ID",
            "uuid": "FPbase ID",
            "seq": "Sequence",
            "name": "Name or Alias",
            "agg": "Oligomerization",
            "primary_reference__year": "Year published",
            "default_state__bleach_measurements__rate": "Photostability (s)",
            "primary_reference__author__family": "Author",
        }

    def name_or_alias_icontains(self, queryset, name, value):
        return queryset.filter(name__icontains=value) | queryset.filter(aliases__icontains=value)

    def switch_type__notequal(self, queryset, name, value):
        return queryset.exclude(switch_type=value)

    def cofactor__notequal(self, queryset, name, value):
        return queryset.exclude(cofactor=value)

    def parent_organism__notequal(self, queryset, name, value):
        return queryset.exclude(parent_organism=value)

    # def name_or_alias_iexact(self, queryset, name, value):
    #     return queryset.filter(name__iexact=value) | queryset.filter(aliases__iexact=value)

    # def name_or_alias_iendswith(self, queryset, name, value):
    #     return queryset.filter(name__iendswith=value) | queryset.filter(aliases__iendswith=value)

    # def name_or_alias_istartswith(self, queryset, name, value):
    #     return queryset.filter(name__istartswith=value) | queryset.filter(aliases__istartswith=value)

    def get_specbright(self, queryset, name, value):
        qsALL = list(queryset.all())
        ids = [P.id for P in qsALL if P.default_state and P.default_state.local_brightness == value]
        return queryset.filter(id__in=ids)

    def get_specbright_lt(self, queryset, name, value):
        qsALL = list(queryset.all())
        ids = [P.id for P in qsALL if P.default_state and P.default_state.local_brightness < value]
        return queryset.filter(id__in=ids)

    def get_specbright_gt(self, queryset, name, value):
        qsALL = list(queryset.all())
        ids = [P.id for P in qsALL if P.default_state and P.default_state.local_brightness > value]
        return queryset.filter(id__in=ids)

    def translate_cdna(self, queryset, name, value):
        return queryset.filter(seq__icontains=Seq.translate(value))
