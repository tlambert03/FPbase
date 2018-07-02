import re
from django import forms
from ..util.importers import check_chroma_for_part, check_semrock_for_part, add_filter_to_database
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from ..models import (Light, Camera, Microscope, Filter, OpticalConfig, FilterPlacement)
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Div, HTML


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
        widget=forms.Textarea(attrs={'rows': 6, 'cols': 40}),
        help_text=('Optical configurations represent a set of filters in your '
                   'microscope, (usually for a specific channel).  See help below')
    )

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        self.helper = FormHelper()
        self.helper.layout = Layout(
            Div('name'),
            Div('description'),
            Div(
                Div('light_source', css_class='col-sm-6 col-xs-12'),
                Div('detector', css_class='col-sm-6 col-xs-12'),
                css_class='row',
            ),
            Div('optical_configurations'),
            HTML('<input type="submit" class="btn btn-primary mb-4" value="Submit" />'
                 '<div class="text-primary ml-4 hidden" id="spinner" style="top: -12px;'
                 ' position: relative;"><i class="fas fa-cog fa-spin mr-1" '
                 'style="font-size: 1.7rem; top: 5px; position: relative;"></i> '
                 'Building your configuration...</div>')
        )
        super().__init__(*args, **kwargs)

    class Meta:
        model = Microscope
        fields = ('name', 'description',)
        help_texts = {
            'name': 'Name of this microscope or set of filter configurations',
            'description': 'This text will appear below the name on your microscope page'
        }

    def create_oc(self, name, filters):
        if len(filters) == 4:
            bs_ex_reflect = bool(filters.pop())
        else:
            bs_ex_reflect = True
        oc = OpticalConfig.objects.create(name=name, owner=self.user)

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
                self.instance.configs.add(newoc)
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
