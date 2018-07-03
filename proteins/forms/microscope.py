import re
from django import forms
from ..util.importers import check_chroma_for_part, check_semrock_for_part, add_filter_to_database
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from ..models import (Light, Camera, Microscope, Filter, OpticalConfig, FilterPlacement)
from dal import autocomplete
from django.forms.models import inlineformset_factory


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

    optical_configurations = forms.CharField(
        required=False,
        widget=forms.Textarea(attrs={'rows': 6, 'cols': 40, 'class': 'textarea form-control'}),
        help_text=('Optical configurations represent a set of filters in your '
                   'microscope, (usually for a specific channel).  See help below')
    )

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        super().__init__(*args, **kwargs)

    class Meta:
        model = Microscope
        fields = ('name', 'description',)
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
            bs_ex_reflect = bool(filters.pop())
        else:
            bs_ex_reflect = True
        oc = OpticalConfig.objects.create(
            name=name,
            owner=self.user,
            microscope=self.instance)

        def _assign_filt(filt, i):
            if isinstance(filt, int) and i == 0:
                oc.laser = filt
                oc.save()
                return
            if isinstance(filt, FilterPromise):
                filt = filt.fetch(self.user)
            if i == 1:
                fpx = FilterPlacement(filter=filt, config=oc,
                                      reflects=bs_ex_reflect,
                                      path=FilterPlacement.EX)
                fpm = FilterPlacement(filter=filt, config=oc,
                                      reflects=not bs_ex_reflect,
                                      path=FilterPlacement.EM)
                fpx.save()
                fpm.save()
            else:
                _path = FilterPlacement.EX if i < 1 else FilterPlacement.EM
                fp = FilterPlacement(filter=filt, config=oc, path=_path)
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
        self.instance = super().save()
        for row in self.cleaned_data['optical_configurations']:
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
            [self.instance.lights.add(i) for i in
             self.cleaned_data['detector'].all()]
        self.instance.owner = self.user
        self.instance.save()
        return self.instance

    def clean_optical_configurations(self):
        ocs = self.cleaned_data['optical_configurations']
        cleaned = []
        # on update form allow for the same name (case insensitive)
        brackets = re.compile(r'[\[\]\(\)]')

        def _getpromise(fname):
            fp = FilterPromise(fname)
            if fp.is_valid:
                return fp
            else:
                self.add_error(
                    'optical_configurations',
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
                    'optical_configurations',
                    "Lines must have 4 or 5 comma-separated fields but this one "
                    "has {}: {}".format(len(splt), line))
            for n, f in enumerate(splt):
                if n == 0:
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

    class Meta:
        model = OpticalConfig
        fields = ('name', 'laser', 'light', 'camera',)
        help_texts = {'laser': 'overrides light source'}
        widgets = {
            'name': forms.widgets.TextInput(attrs={'class': 'textinput textInput form-control'}),
            'camera': forms.widgets.Select(attrs={'class': 'form-control custom-select'}),
            'light': forms.widgets.Select(attrs={'class': 'form-control custom-select'}),
            'laser': forms.widgets.NumberInput(attrs={'class': 'numberinput form-control'})
        }


class BaseOpticalConfigFormSet(forms.BaseInlineFormSet):
    def clean(self):
        # perform any cross-formset validation here
        super().clean()


OpticalConfigFormSet = inlineformset_factory(
    Microscope, OpticalConfig, form=OpticalConfigForm,
    formset=BaseOpticalConfigFormSet, extra=1, can_delete=True)

