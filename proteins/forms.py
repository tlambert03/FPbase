from django import forms
from django.utils.text import slugify
from django.utils.safestring import mark_safe
from django.forms.models import inlineformset_factory  # ,BaseInlineFormSet
from django.apps import apps
from django.urls import reverse
from django.core.exceptions import ObjectDoesNotExist
from proteins.models import (Protein, State, StateTransition, Spectrum,
                             ProteinCollection, BleachMeasurement)
from proteins.validators import validate_spectrum, validate_doi, protein_sequence_validator
from proteins.util.importers import text_to_spectra
from proteins.util.helpers import zip_wave_data
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Div, HTML, Submit


def popover_html(label, content, side='right'):
    return '<label data-toggle="tooltip" style="padding-' + side + ': 1rem;" data-placement="' + side + '" title="' + content + '">' + label + '</label>'


def check_existence(form, fieldname, value):
    # check whether another Protein already has this
    # value for this fieldname
    # return value if not
    if not value:
        return None

    # on update form allow for the same sequence (case insensitive)
    if hasattr(form, 'instance'):
        instanceVal = getattr(form.instance, fieldname)
        if isinstance(instanceVal, str) and instanceVal.upper() == value.upper():
            return value

    if fieldname == 'name':
        slug = slugify(value)
        query = Protein.objects.filter(slug=slug)
        query = query | Protein.objects.filter(name__iexact=value)
        query = query | Protein.objects.filter(name__iexact=value.replace(' ', ''))
        query = query | Protein.objects.filter(name__iexact=value.replace(' ', '').replace('monomeric', 'm'))
    else:
        query = Protein.objects.filter(**{fieldname: value})

    if query.exists():
        prot = query.first()
        raise forms.ValidationError(mark_safe(
            '<a href="{}" style="text-decoration: underline;">{}</a> already has this {}'.format(
                prot.get_absolute_url(), prot.name, Protein._meta.get_field(fieldname).verbose_name.lower())))
    return value


class DOIField(forms.CharField):
    max_length = 100
    required = True
    default_validators = [validate_doi]

    def to_python(self, value):
        if value and isinstance(value, str):
            value = value.lstrip('https://doi.org/')
        return super().to_python(value)


class SequenceField(forms.CharField):
    widget = forms.Textarea(attrs={
        'class': 'vLargeTextField',
        'rows': 3,
    })
    max_length = 1024
    label = 'AA Sequence'
    default_validators = [protein_sequence_validator]
    strip = True

    def to_python(self, value):
        if value and isinstance(value, str):
            value = value.replace(' ', '').upper()
        return super().to_python(value)


class SelectAddWidget(forms.widgets.Select):
    template_name = 'proteins/forms/widgets/select_add.html'


class ProteinForm(forms.ModelForm):
    '''Form class for user-facing protein creation/submission form '''
    reference_doi = DOIField(required=False, help_text='e.g. 10.1038/nmeth.2413', label='Reference DOI')
    seq = SequenceField(required=False, help_text='Amino acid sequence (IPG ID is preferred)',
        label=popover_html('Sequence', "If you enter an IPG ID, the sequence can be automatically fetched from NCBI"),)
    # reference_pmid = forms.CharField(max_length=24, label='Reference Pubmed ID',
    #     required=False, help_text='e.g. 23524392 (must provide either DOI or PMID)')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.form_tag = False
        self.helper.error_text_inline = True
        self.helper.layout = Layout(
            Div(
                Div('name', css_class='col-md-4 col-sm-12'),
                Div('aliases', css_class='col-md-4 col-sm-6'),
                Div('reference_doi', css_class='col-md-4 col-sm-6'),
                css_class='row',
            ),
            Div(
                Div('agg', css_class='col-sm-6'),
                Div('parent_organism', css_class='col-sm-6'),
                css_class='row',
            ),
            Div(
                Div('ipg_id', css_class='col-lg-3 col-sm-6'),
                Div('genbank', css_class='col-lg-3 col-sm-6'),
                Div('uniprot', css_class='col-lg-3 col-sm-6'),
                Div('pdb', css_class='col-lg-3 col-sm-6'),
                css_class='row',
            ),
            Div(
                Div('seq', css_class='col'),
                css_class='row',
            )
        )

    class Meta:
        model = Protein
        fields = ('name', 'ipg_id', 'seq', 'agg', 'parent_organism', 'pdb', 'reference_doi',
                  'aliases', 'chromophore', 'genbank', 'uniprot', 'blurb')
        labels = {
            "agg": "Oligomerization",
        }
        help_texts = {
            'aliases': 'Comma separated list of aliases',
            'pdb': 'Comma separated list of <a href="https://www.rcsb.org/pdb/staticHelp.do?p=help/advancedsearch/pdbIDs.html" target="_blank">PDB IDs</a>',
            'ipg_id': 'NCBI <a href="https://www.ncbi.nlm.nih.gov/ipg/docs/about/" target="_blank">Identical Protein Group ID</a>',
            'genbank': 'NCBI <a href="https://www.ncbi.nlm.nih.gov/genbank/sequenceids/" target="_blank">GenBank ID</a>',
            'uniprot': '<a href="https://www.uniprot.org/help/accession_numbers" target="_blank">UniProt accession number</a>'
        }
        widgets = {
            'parent_organism': SelectAddWidget(),
        }

    def clean_name(self):
        return check_existence(self, 'name', self.cleaned_data['name'])

    def clean_seq(self):
        return check_existence(self, 'seq', self.cleaned_data['seq'])

    def clean_ipg_id(self):
        return check_existence(self, 'ipg_id', self.cleaned_data['ipg_id'])

    def clean_genbank(self):
        return check_existence(self, 'genbank', self.cleaned_data['genbank'])

    def clean_uniprot(self):
        return check_existence(self, 'uniprot', self.cleaned_data['uniprot'])

    def save_new_only(self, commit=True):
        # check the current db Instance for all of the changed_data values
        # if there is currently a non-null value in the database,
        # don't overwrite it ...
        # repopulate self.instance with form data afterwards, so that calling
        # self.save() again in the future WILL overwrite the database values
        isupdate = bool(hasattr(self, 'instance') and self.instance.pk)
        backup = {}
        if not isupdate:
            return super().save(commit=commit)
        else:
            dbInstance = Protein.objects.get(pk=self.instance.pk)
            for field in self.changed_data:
                if not hasattr(dbInstance, field):
                    continue
                dbValue = getattr(dbInstance, field)
                formValue = getattr(self.instance, field)
                if field in ('pdb', 'aliases'):
                    continue
                    for val in formValue:
                        if val not in dbValue:
                            getattr(self.instance, field).append(val)
                else:
                    if dbValue and formValue != dbValue:
                        backup[field] = formValue
                        setattr(self.instance, field, dbValue)
        self.instance = super().save(commit=commit)
        for field, value in backup.items():
            setattr(self.instance, field, value)
        return self.instance


class SpectrumFormField(forms.CharField):
    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={
        'class': 'vLargeTextField',
        'rows': 2,
    })

    def __init__(self, *args, **kwargs):
        if 'help_text' not in kwargs:
            kwargs['help_text'] = 'List of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]. File data takes precedence.'
        super().__init__(*args, **kwargs)


class SpectrumForm(forms.ModelForm):
    lookup = {
        Spectrum.DYE: ('owner_dye', 'Dye'),
        Spectrum.PROTEIN: ('owner_state', 'State'),
        Spectrum.FILTER: ('owner_filter', 'Filter'),
        Spectrum.CAMERA: ('owner_camera', 'Camera'),
        Spectrum.LIGHT: ('owner_light', 'Light')
    }

    owner = forms.CharField(max_length=100, label='Owner Name', required=True, help_text='Name of protein, dye, filter, etc...')
    data = SpectrumFormField(required=False, label='Data')
    file = forms.FileField(required=False, label='File Upload',
                           help_text='2 column CSV/TSV file with wavelengths in first column and data in second column')

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        self.helper = FormHelper()
        self.helper.attrs = {
            "id": 'spectrum-submit-form',
            "data-validate-owner-url": reverse('proteins:validate_spectrumownername')
        }
        self.helper.add_input(Submit('submit', 'Submit'))
        self.helper.layout = Layout(
            Div('owner'),
            Div(
                Div('category', css_class='col-sm-6 col-xs-12'),
                Div('subtype', css_class='col-sm-6 col-xs-12'),
                css_class='row',
            ),
            Div('file', 'data',),
            Div(
                Div('ph', css_class='col-md-6 col-sm-12'),
                Div('solvent', css_class='col-md-6 col-sm-12'),
                css_class='row',
            ),
        )

        super().__init__(*args, **kwargs)

    class Meta:
        model = Spectrum
        fields = ('category', 'subtype', 'file', 'data', 'ph', 'solvent', 'owner')
        widgets = {'data': forms.Textarea(attrs={'class': 'vLargeTextField', 'rows': 2})}

    def clean(self):
        cleaned_data = super().clean()
        if not (cleaned_data.get('data') or self.files):
            self.add_error('data', 'Please either fill in the data field or '
                           'select a file to upload.')
            self.add_error('file', 'Please either fill in the data field or '
                           'select a file to upload.')

    def save(self, commit=True):
        owner_name = self.cleaned_data.get("owner")
        cat = self.cleaned_data.get('category')
        owner_model = apps.get_model('proteins', self.lookup[cat][1])
        if cat == Spectrum.PROTEIN:
            owner_obj = owner_model.objects.get(protein__name__iexact=owner_name)
        else:
            owner_obj, c = owner_model.objects.get_or_create(name=owner_name, defaults={'created_by': self.user})
            if not c:
                owner_obj.update_by = self.user
                owner_obj.save()
        setattr(self.instance, self.lookup[cat][0], owner_obj)
        self.instance.created_by = self.user
        return super().save(commit=commit)

    def clean_file(self):
        if self.files:
            filetext = ''
            try:
                for chunk in self.files['file'].chunks():
                    filetext += chunk.decode("utf-8")
                x, y, headers = text_to_spectra(filetext)
                if not len(y):
                    self.add_error('file', 'Did not find a data column in the provided file')
                if not len(x):
                    self.add_error('file', 'Could not parse wavelengths from first column')
            except Exception:
                self.add_error('file', 'Sorry, could not parse spectrum from this file. '
                               'Is it it two column csv with (wavelength, spectrum)?')
            if not self.errors:
                self.cleaned_data['data'] = zip_wave_data(x, y[0])
                self.data = self.data.copy()
                self.data['data'] = self.cleaned_data['data']

    def clean_owner(self):
        # make sure an owner with the same category and name doesn't already exist
        owner = self.cleaned_data.get("owner")
        cat = self.cleaned_data.get('category')
        stype = self.cleaned_data.get('subtype')
        try:
            mod = apps.get_model('proteins', self.lookup[cat][1])
            if cat == Spectrum.PROTEIN:
                obj = mod.objects.get(protein__slug=slugify(owner))
            else:
                obj = mod.objects.get(slug=slugify(owner))
        except ObjectDoesNotExist:
            if cat == Spectrum.PROTEIN:
                recs = Protein.objects.find_similar(owner, 0.5)
                tmplt = "<a href='{}'>{}</a>"
                sim_link = ", ".join(
                    [tmplt.format(p.get_absolute_url(), p.name) for p in recs])
                sim_link = "(Similar: {}) ".format(sim_link) if recs.count() else ''
                self.add_error('owner', forms.ValidationError(
                    mark_safe("Could not find {} in the database. "
                              .format(owner) + sim_link +
                              "<a href='{}' style='text-decoration: underline;'>"
                              .format(reverse('proteins:submit')) +
                              "Add it</a>?"),
                    params={'owner': owner}, code='no_protein_exists'))
            # new object will be made in save()
            return owner
        except KeyError:
            # this might be repetitive... since a missing category will already
            # throw an error prior to this point
            if not cat:
                raise forms.ValidationError("Category not provided")
            else:
                raise forms.ValidationError("Category not recognized")
        else:
            # object exists... check if it has this type of spectrum
            if obj.spectra.filter(subtype=stype).exists():
                if cat == Spectrum.PROTEIN:
                    self.add_error('owner', forms.ValidationError(
                        "Protein %(owner)s already has a %(stype)s spectrum",
                        params={
                            'owner': obj.protein.name,
                            'stype': obj.spectra.filter(subtype=stype).first().get_subtype_display()},
                        code='owner_exists'))

                    # raise forms.ValidationError(
                    #     "The default state for protein %(owner)s already has a spectrum of type %(stype)s.",
                    #     params={'owner': owner, 'stype': stype},
                    #     code='owner_exists')
                else:
                    self.add_error('owner', forms.ValidationError(
                        "A %(model)s with the name %(name)s already has a %(stype)s spectrum",
                        params={'model': self.lookup[cat][1].lower(),
                                'name': obj.name,
                                'stype': obj.spectra.filter(subtype=stype).first().get_subtype_display()},
                        code='owner_exists'))
                    # raise forms.ValidationError(
                    #     "A %(model)s with the slug %(slug)s already has a spectrum of type %(stype)s.",
                    #     params={'model': self.lookup[cat][1].lower(), 'slug': slugify(owner), 'stype': stype},
                    #     code='owner_exists')
            else:
                return owner


class StateForm(forms.ModelForm):
    # ex_spectra = SpectrumFormField(required=False, label='Excitation Spectrum')
    # em_spectra = SpectrumFormField(required=False, label='Emission Spectrum')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True
        self.helper.layout = Layout(
            Div(
                HTML('<h4>State</h4>'),
                Div(
                    Div('name', css_class='col-md-6'),
                    Div('is_dark', css_class='col-md-6 self-align-center dark_state_button'),
                    css_class='row',
                ),
                Div(
                    Div('ex_max', css_class='col-lg-3 col-md-6'),
                    Div('em_max', css_class='col-lg-3 col-md-6'),
                    Div('ext_coeff', css_class='col-lg-3 col-md-6'),
                    Div('qy', css_class='col-lg-3 col-md-6'),
                    css_class='row hide_if_dark',
                ),
                Div(
                    Div('pka', css_class='col-md-4'),
                    Div('maturation', css_class='col-md-4'),
                    Div('lifetime', css_class='col-md-4'),
                    css_class='row hide_if_dark',
                ),
                #Div(
                #    Div('ex_spectra', css_class='col-md-6 spectrum-field'),
                #    Div('em_spectra', css_class='col-md-6 spectrum-field'),
                #    css_class='row hide_if_dark',
                #),
                css_class='stateform_block'
            )
        )

    def clean_ex_max(self):
        ex_max = self.cleaned_data.get("ex_max")
        if not self.cleaned_data['is_dark'] and not ex_max:
            raise forms.ValidationError("Must provide Ex Max for non-dark state")
        return ex_max

    def clean_em_max(self):
        em_max = self.cleaned_data.get("em_max")
        if not self.cleaned_data['is_dark'] and not em_max:
            raise forms.ValidationError("Must provide Em Max for non-dark state")
        return em_max

    class Meta:
        model = State
        fields = ('name', 'is_dark', 'ex_max', 'em_max', 'ext_coeff', 'qy', 'protein',
                  'pka', 'maturation', 'lifetime',
                  # 'ex_spectra', 'em_spectra',
                  )
        # widgets = {
        #    'ex_spectra': forms.Textarea(attrs={'rows': 2}),
        #    'em_spectra': forms.Textarea(attrs={'rows': 2}),
        # }
        labels = {
            "ex_max": "Excitation Max (nm)",
            "em_max": "Emission Max (nm)",
            "pka": "pKa",
        }


class BaseStateFormSet(forms.BaseInlineFormSet):
    def clean(self):
        if any(self.errors):
            # Don't bother validating the formset unless each form is valid on its own
            return
        names = []
        darkcount = 0
        for form in self.forms:
            if not form.cleaned_data:
                continue
            name = form.cleaned_data.get('name')
            if name in names:
                raise forms.ValidationError("Different states must have distinct names.")
            names.append(name)
            darkcount += 1 if form.cleaned_data.get('is_dark') else 0
            if darkcount > 1:
                raise forms.ValidationError("A protein can only have a single dark state")


StateFormSet = inlineformset_factory(Protein, State, form=StateForm, formset=BaseStateFormSet, extra=1, can_delete=True)


class StateTransitionForm(forms.ModelForm):

    class Meta:
        model = StateTransition
        fields = ('trans_wave', 'protein', 'from_state', 'to_state',)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True
        self.helper.layout = Layout(
            Div(
                Div(
                    Div('from_state', css_class='col-md-4'),
                    Div('to_state', css_class='col-md-4'),
                    Div('trans_wave', css_class='col-md-4'),
                    css_class='row',
                ),
            )
        )


StateTransitionFormSet = inlineformset_factory(Protein, StateTransition, form=StateTransitionForm, extra=1)


class CollectionForm(forms.ModelForm):

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop('user', None)
        self.proteins = kwargs.pop('proteins', None)
        super().__init__(*args, **kwargs)

    class Meta:
        model = ProteinCollection
        fields = ('name', 'description', 'private')

    def clean_name(self):
        name = self.cleaned_data['name']
        # on update form allow for the same name (case insensitive)
        isupdate = bool(hasattr(self, 'instance') and self.instance.pk)
        if isupdate and self.instance.name.lower() == name.lower():
            return name
        try:
            col = ProteinCollection.objects.get(name__iexact=name, owner=self.user)
        except ProteinCollection.DoesNotExist:
            return name

        raise forms.ValidationError(mark_safe(
            'You already have a collection named <a href="{}" style="text-decoration: underline;">{}</a>'.format(
                col.get_absolute_url(), col.name)))


class BleachMeasurementForm(forms.ModelForm):
    reference_doi = DOIField(required=False, help_text='e.g. 10.1038/nmeth.2413', label='Reference DOI')

    class Meta:
        model = BleachMeasurement
        fields = ('rate', 'power', 'units', 'modality', 'reference_doi',
                  'state', 'light', 'temp', 'fusion', 'in_cell')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True
        self.helper.layout = Layout(
            Div(
                Div(
                    Div('state', css_class='col-md-4 col-xs-12'),
                    Div('rate', css_class='col-md-4 col-xs-12'),
                    Div('reference_doi', css_class='col-md-4 col-xs-12'),
                    Div('light', css_class='col-md-4 col-xs-12'),
                    Div('power', css_class='col-md-4 col-xs-12'),
                    Div('units', css_class='col-md-4 col-xs-12'),
                    Div('modality', css_class='col-md-3 col-xs-12'),
                    Div('temp', css_class='col-md-3 col-xs-12'),
                    Div('fusion', css_class='col-md-3 col-xs-12'),
                    Div('in_cell', css_class='col-md-3 col-xs-12'),
                    css_class='row',
                ),
            )
        )

