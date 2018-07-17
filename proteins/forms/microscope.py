import re
from django import forms
from ..util.importers import check_chroma_for_part, check_semrock_for_part, add_filter_to_database
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from ..models import (Light, Camera, Microscope, Filter, OpticalConfig, FilterPlacement, ProteinCollection)
from dal import autocomplete
from django.forms.models import inlineformset_factory
from collections import defaultdict


class FilterPromise(object):
    def __init__(self, part):
        self.part = part

    @property
    def is_valid(self):
        if hasattr(self, '_valid'):
            return self._valid
        if check_chroma_for_part(self.part):
            self.brand = 'chroma'
            self._valid = True
        elif check_semrock_for_part(self.part):
            self.brand = 'semrock'
            self._valid = True
        else:
            self.brand = None
            self._valid = False
        return self._valid

    def fetch(self, user=None):
        newObjects, _ = add_filter_to_database(self.brand, self.part, user)
        return newObjects[0].owner


class MicroscopeForm(forms.ModelForm):

    light_source = forms.ModelMultipleChoiceField(
        label='Light Source(s)',
        queryset=Light.objects.all(), required=False,
        help_text=('Specify laser lines individually in the optical configurations below. '
                   'If a single light source is selected here, it will be applied to '
                   'all configurations without lasers. If your light source is missing, '
                   'you can <a href="/spectra/submit"> submit a power density curve</a>'))

    detector = forms.ModelMultipleChoiceField(
        label='Detector(s)',
        queryset=Camera.objects.all(), required=False,
        help_text=(
            'If a single detector is chosen, it will be applied to all configurations. '
            'If your camera is missing, you can <a href="/spectra/submit">'
            'submit a QE curve</a>'
        )
    )

    optical_configs = forms.CharField(
        label="Optical Configurations",
        required=False,
        widget=forms.Textarea(attrs={'rows': 6, 'cols': 40, 'class': 'textarea form-control'}),
        help_text=('See extended help below')
    )

    collection = forms.ModelChoiceField(
        required=False, queryset=ProteinCollection.objects.exclude(private=True),
        label='Protein Collection',
        help_text='Subset of probes to show on microscope page. '
                  'Leave blank to enable all fluorophores in the database.')

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        super().__init__(*args, **kwargs)

    class Meta:
        model = Microscope
        fields = ('name', 'description', 'collection')
        help_texts = {
            'name': 'Name of this microscope or set of filter configurations',
            'description': 'This text will appear below the name on your microscope page'
        }
        widgets = {
            'name': forms.widgets.TextInput(attrs={'class': 'textinput textInput form-control'}),
            'description': forms.widgets.TextInput(attrs={'class': 'textinput textInput form-control'}),
            'detector': forms.widgets.Select(attrs={'class': 'selectmultiple form-control'}),
            'light_source': forms.widgets.Select(attrs={'class': 'selectmultiple form-control'}),
        }

    def create_oc(self, name, filters):
        if len(filters) == 4:
            bs_em_reflect = not bool(filters.pop())
        else:
            bs_em_reflect = False
        oc = OpticalConfig.objects.create(
            name=name,
            owner=self.user,
            microscope=self.instance)

        _paths = [FilterPlacement.EX, FilterPlacement.BS, FilterPlacement.EM]

        def _assign_filt(filt, i):
            if isinstance(filt, int) and i == 0:
                oc.laser = filt
                oc.save()
                return
            if isinstance(filt, FilterPromise):
                filt = filt.fetch(self.user)

            fp = FilterPlacement(filter=filt, config=oc, path=_paths[i],
                                 reflects=bs_em_reflect if i == 1 else False)
            fp.save()

        for i, fnames in enumerate(filters):
            if not fnames:
                continue
            if not isinstance(fnames, (tuple, list)):
                fnames = [fnames]
            for _f in fnames:
                _assign_filt(_f, i)

        return oc

    def save(self, commit=True):
        self.instance = super().save(commit=False)
        for row in self.cleaned_data['optical_configs']:
            newoc = self.create_oc(row[0], row[1:])
            if newoc:
                if self.cleaned_data['light_source'].count() == 1:
                    if not newoc.laser:
                        newoc.light = self.cleaned_data['light_source'].first()
                if self.cleaned_data['detector'].count() == 1:
                    newoc.camera = self.cleaned_data['detector'].first()
                newoc.save()
                self.instance.optical_configs.add(newoc)
        if self.cleaned_data['light_source'].count() > 1:
            [self.instance.lights.add(i) for i in
             self.cleaned_data['light_source'].all()]
        if self.cleaned_data['detector'].count() > 1:
            [self.instance.cameras.add(i) for i in
             self.cleaned_data['detector'].all()]
        self.instance.owner = self.user
        if commit:
            self.instance.save()
        return self.instance

    def clean_optical_configs(self):
        ocs = self.cleaned_data['optical_configs']
        cleaned = []
        namestore = []
        if self.instance:
            namestore.extend([oc.name for oc in self.instance.optical_configs.all()])
        # on update form allow for the same name (case insensitive)
        brackets = re.compile(r'[\[\]\(\)]')

        def _getpromise(fname):
            fp = FilterPromise(fname)
            if fp.is_valid:
                return fp
            else:
                self.add_error(
                    'optical_configs',
                    'Filter not found in database or at Chroma/Semrock: '
                    '{}'.format(fname))
                return None

        def lookup(fname):
            if not fname:
                return None
            if fname.isdigit():
                return int(fname)
            try:
                return Filter.objects.get(name__icontains=fname)
            except MultipleObjectsReturned:
                try:
                    return Filter.objects.get(part__iexact=fname)
                except ObjectDoesNotExist:
                    return _getpromise(fname)
            except ObjectDoesNotExist:
                return _getpromise(fname)
            return None

        for line in ocs.splitlines():
            if not line:
                continue
            _out = []
            if brackets.search(line):
                _splt = [i.strip() for i in re.split(r'(\([^)]*\))', line) if i.strip()]
                splt = []
                for item in _splt:
                    if brackets.search(item):
                        splt.append([n.strip() for n in brackets.sub('', item).split(',') if n.strip()])
                    else:
                        if item.endswith(','):
                            item = item[:-1]
                        if item.startswith(','):
                            item = item[1:]
                        splt.extend([n.strip() for n in item.split(',')])
            else:
                splt = [i.strip() for i in line.split(',')]
            if not len(splt) in (4, 5):
                self.add_error(
                    'optical_configs',
                    "Lines must have 4 or 5 comma-separated fields but this one "
                    "has {}: {}".format(len(splt), line))
            for n, f in enumerate(splt):
                if n == 0:
                    if f in namestore:
                        self.add_error(
                            'optical_configs',
                            'Optical config with the name %s already exists.' % f)
                    else:
                        namestore.append(f)
                        _out.append(f)
                elif n == 4:
                    if f.lower() in ('0', 'false', 'none'):
                        _out.append(False)
                    else:
                        _out.append(True)
                else:
                    if isinstance(f, list):
                        _out.append([lookup(x) for x in f])
                    else:
                        _out.append(lookup(f))
            cleaned.append(_out)
        return cleaned


class OpticalConfigForm(forms.ModelForm):

    ex_filters = forms.ModelMultipleChoiceField(
        label="Excitation Filter(s)",
        queryset=Filter.objects.all(), required=False,
        widget=autocomplete.ModelSelect2Multiple(
            url='proteins:filter-autocomplete',
            attrs={
                'class': 'custom-select',
                'data-theme': 'bootstrap',
                'data-width': "100%",
                'data-placeholder': '----------',
            }
        ),
    )
    em_filters = forms.ModelMultipleChoiceField(
        label="Emission Filter(s)",
        queryset=Filter.objects.all(), required=False,
        widget=autocomplete.ModelSelect2Multiple(
            url='proteins:filter-autocomplete',
            attrs={
                'data-theme': 'bootstrap',
                'data-width': "100%",
                'data-placeholder': '----------',
            }
        ),
    )
    bs_filters = forms.ModelMultipleChoiceField(
        label="Dichroic Filter(s)",
        queryset=Filter.objects.all(), required=False,
        widget=autocomplete.ModelSelect2Multiple(
            url='proteins:filter-autocomplete',
            attrs={
                'data-theme': 'bootstrap',
                'data-width': "100%",
                'data-placeholder': '----------',
            }
        ),
    )
    invert_bs = forms.BooleanField(
        required=False,
        widget=forms.widgets.CheckboxInput(attrs={'class': 'custom-control-input'})
    )

    def __init__(self, *args, **kwargs):
        instance = kwargs.get('instance', None)
        if instance:
            kwargs.update(
                initial={
                    'invert_bs': instance.inverted_bs.exists(),
                    'bs_filters': instance.bs_filters,
                    'em_filters': instance.em_filters,
                    'ex_filters': instance.ex_filters,
                }
            )
        super().__init__(*args, **kwargs)

    def save(self, commit=True):
        oc = super().save()

        for _p in ('em_filters', 'ex_filters', 'bs_filters'):
            for filt in self.initial.get(_p, []):
                if filt not in self.cleaned_data[_p].values_list('id', flat=True):
                    oc.filterplacement_set.filter(filter=filt).delete()

        for _p in ('em_filters', 'ex_filters', 'bs_filters'):
            for filt in self.cleaned_data.get(_p, []):
                if filt.id not in self.initial.get(_p, []):
                    path = _p.split('_')[0]
                    reflects = self.cleaned_data['invert_bs'] if path == 'bs' else False
                    oc.add_filter(filt, path, reflects)

        # for filt in self.cleaned_data['em_filters']:
        #     if filt.id not in self.initial.get('em_filters', []):
        #         oc.add_em_filter(filt)
        # for filt in self.cleaned_data['ex_filters']:
        #     if filt.id not in self.initial.get('ex_filters', []):
        #         oc.add_ex_filter(filt)
        # for filt in self.cleaned_data['bs_filters']:
        #     if filt.id not in self.initial.get('bs_filters', []):
        #         oc.add_bs_filter(filt, self.cleaned_data['invert_bs'])
        if self.initial:
            if self.cleaned_data['invert_bs'] != self.initial.get('invert_bs'):
                for filt in self.cleaned_data['bs_filters']:
                    for fp in oc.filterplacement_set.filter(
                            filter=filt, path=FilterPlacement.BS):
                        fp.reflects = self.cleaned_data['invert_bs']
                        fp.save()
        if commit:
            oc.save()
        return oc

    class Meta:
        model = OpticalConfig
        fields = ('name', 'laser', 'light', 'camera')
        help_texts = {
            'light': 'laser overrides light source',
            'name': 'name of this optical config'
        }
        widgets = {
            'name': forms.widgets.TextInput(attrs={'class': 'textinput textInput form-control'}),
            'camera': forms.widgets.Select(attrs={'class': 'form-control custom-select'}),
            'light': forms.widgets.Select(attrs={'class': 'form-control custom-select'}),
            'laser': forms.widgets.NumberInput(attrs={'class': 'numberinput form-control'}),
        }

    def is_valid(self):
        if hasattr(self, '_isvalid'):
            return self._isvalid
        else:
            self._isvalid = super().is_valid()
            if not self._isvalid:
                for field, error in self.errors.items():
                    if field in self.fields:
                        self.fields[field].widget.attrs['class'] += " is-invalid"
        return self._isvalid


class BaseOpticalConfigFormSet(forms.BaseInlineFormSet):
    pass


OpticalConfigFormSet = inlineformset_factory(
    Microscope, OpticalConfig, form=OpticalConfigForm,
    formset=BaseOpticalConfigFormSet, extra=1, can_delete=True)

