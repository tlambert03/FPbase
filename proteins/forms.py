from django import forms
from django.contrib.auth.models import User
from django.utils.text import slugify
from django.utils.safestring import mark_safe
from django.core.exceptions import ValidationError
from django.forms.models import inlineformset_factory  # ,BaseInlineFormSet
from proteins.models import Protein, State
from proteins.validators import validate_spectrum, validate_doi
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Div, HTML


def popover_html(label, content, side='right'):
    return '<label data-toggle="tooltip" style="padding-' + side + ': 1rem;" data-placement="' + side + '" title="' + content + '">' + label + '</label>'


class ProteinUpdateForm(forms.ModelForm):
    '''Form class for user-facing protein creation/submission form '''
    reference_doi = forms.CharField(max_length=100, label='Reference DOI',
        required=True, validators=[validate_doi],
        help_text='e.g. 10.1038/nmeth.2413')
    # reference_pmid = forms.CharField(max_length=24, label='Reference Pubmed ID',
    #     required=False, help_text='e.g. 23524392 (must provide either DOI or PMID)')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper(self)
        self.helper.form_tag = False
        self.helper.error_text_inline = True
        self.helper.layout = Layout(
            Div(
                Div('name', css_class='col-sm-6'),
                Div('reference_doi', css_class='col-sm-6'),
                css_class='row',
            ),
            Div(
                Div('ipg_id', css_class='col-sm-4'),
                Div('agg', css_class='col-sm-4'),
                Div('parent_organism', css_class='col-sm-4'),
                css_class='row',
            ),
            Div(
                Div('seq', css_class='col'),
                css_class='row',
            )
        )

    class Meta:
        model = Protein
        fields = ('name', 'ipg_id', 'seq', 'agg', 'parent_organism', 'reference_doi')
        widgets = {
            'seq': forms.Textarea(attrs={'rows': 3}),
        }
        labels = {
            "seq": popover_html('Sequence', "If you enter an IPG ID, the sequence can be automatically fetched from NCBI"),
            "agg": "Oligomerization",
        }
        help_texts = {
            'ipg_id': 'NCBI <a href="https://www.ncbi.nlm.nih.gov/ipg/docs/about/">Identical Protein Group ID</a>',
        }


class ProteinSubmitForm(ProteinUpdateForm):
    def clean_name(self):
        name = self.cleaned_data['name']
        slug = slugify(name)
        query = Protein.objects.filter(slug=slug)
        query = query | Protein.objects.filter(name__iexact=name)
        query = query | Protein.objects.filter(name__iexact=name.replace(' ', ''))
        query = query | Protein.objects.filter(name__iexact=name.replace(' ', '').replace('monomeric', 'm'))
        if query.exists():
            prot = query.first()
            raise ValidationError(mark_safe(
                '<a href="{}" style="text-decoration: underline;">{}</a> already exists in the database'.format(
                    prot.get_absolute_url(), prot.name)))
        return name

    def clean_seq(self):
        seq = self.cleaned_data['seq']
        if not seq:
            return None
        query = Protein.objects.filter(seq=seq)
        if query.exists():
            prot = query.first()
            raise ValidationError(mark_safe(
                '<a href="{}" style="text-decoration: underline;">{}</a> already has this sequence'.format(
                    prot.get_absolute_url(), prot.name)))
        return seq

    def clean_ipg_id(self):
        ipg_id = self.cleaned_data['ipg_id']
        if not ipg_id:
            return None
        query = Protein.objects.filter(ipg_id=ipg_id)
        if query.exists():
            prot = query.first()
            raise ValidationError(mark_safe(
                '<a href="{}" style="text-decoration: underline;">{}</a> already has this ID'.format(
                    prot.get_absolute_url(), prot.name)))
        return ipg_id


class SpectrumFormField(forms.CharField):
    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={
        'class': 'vLargeTextField',
        'rows': 2,
    })

    def __init__(self, *args, **kwargs):
        if 'help_text' not in kwargs:
            kwargs['help_text'] = 'List of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]'
        super().__init__(*args, **kwargs)


class StateSubmitForm(forms.ModelForm):
    ex_spectra = SpectrumFormField(required=False, label='Excitation Spectrum')
    em_spectra = SpectrumFormField(required=False, label='Emission Spectrum')

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
                Div(
                    Div('ex_spectra', css_class='col-md-6 spectrum-field'),
                    Div('em_spectra', css_class='col-md-6 spectrum-field'),
                    css_class='row hide_if_dark',
                ),
                css_class='stateform_block'
            )
        )

    def clean(self):
        super().clean()
        is_dark = self.cleaned_data.get("is_dark")
        ex_max = self.cleaned_data.get("ex_max")
        em_max = self.cleaned_data.get("em_max")
        if not is_dark and not (ex_max and em_max):
            raise forms.ValidationError("Must provide both ex & em maxima for non-dark state")

    class Meta:
        model = State
        fields = ('name', 'is_dark', 'ex_max', 'em_max', 'ext_coeff', 'qy',
                  'pka', 'maturation', 'lifetime', 'ex_spectra', 'em_spectra')
        widgets = {
            'ex_spectra': forms.Textarea(attrs={'rows': 2}),
            'em_spectra': forms.Textarea(attrs={'rows': 2}),
        }
        labels = {
            "ext_coeff": "Extinction Coefficient",
            "ex_max": "Excitation Max (nm)",
            "em_max": "Emission Max (nm)",
            "qy": "Quantum Yield",
            "pka": "pKa",
            "ex_spectra": popover_html('Excitation Spectrum', "If you enter an IPG ID, the sequence can be automatically fetched from NCBI"),
            "em_spectra": "Emission Spectrum",
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


StateFormSet = inlineformset_factory(Protein, State, form=StateSubmitForm, formset=BaseStateFormSet, extra=1)
StateUpdateFormSet = inlineformset_factory(Protein, State, form=StateSubmitForm, formset=BaseStateFormSet, extra=1, can_delete=True)


class ProteinForm(forms.ModelForm):
    class Meta:
        model = Protein
        # fields = ['name', 'slug', 'gb_prot', 'gb_nuc', 'seq', 'parent_organism', 'primary_reference', 'references', 'FRET_partner', 'created_by', 'updated_by']
        fields = ['name', 'slug', 'gb_prot', 'gb_nuc', 'seq', 'parent_organism', 'FRET_partner', 'created_by', 'updated_by']

    def clean_created_by(self):
        if not self.cleaned_data['created_by']:
            return User()
        return self.cleaned_data['created_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']


class StateForm(forms.ModelForm):
    ex_spectra = SpectrumFormField(required=False)
    em_spectra = SpectrumFormField(required=False)

    class Meta:
        model = State
        fields = ['protein', 'is_dark', 'name', 'ex_max', 'em_max', 'ex_spectra', 'em_spectra', 'ext_coeff', 'qy', 'pka', 'maturation', 'lifetime', 'created_by', 'updated_by']

    def __init__(self, *args, **kwargs):
        super(StateForm, self).__init__(*args, **kwargs)  # populates the post
        try:
            self.fields['to_state'].queryset = State.objects.filter(protein=self.instance.protein).exclude(slug=self.instance.slug)
            if self.instance.protein.switch_type == '1':
                pass
                # would like to remove fields from basic type proteins
                # del self.fields['to_state']
        except Exception:
            #FIXME: update state business
            pass
            # self.fields['to_state'].queryset = State.objects.filter(name='')

    def clean_created_by(self):
        if not self.cleaned_data['created_by']:
            return User()
        return self.cleaned_data['created_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']
