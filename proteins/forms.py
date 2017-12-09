from django.contrib.auth.models import User
from django import forms
from django.forms.models import BaseInlineFormSet
from django.utils.translation import ugettext_lazy as _
from proteins.models import Protein, State
from proteins.validators import validate_spectrum
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, ButtonHolder, Submit, Div, Fieldset
from crispy_forms.bootstrap import Accordion, AccordionGroup


class ProteinSubmitForm(forms.ModelForm):
    class Meta:
        model = Protein
        exclude = ['slug', 'base_name', 'FRET_partner', 'mw']


class ProteinSearchForm(forms.ModelForm):
    ex_max = forms.IntegerField(required=False, help_text='Search for protein with Ex max around this wavelength')
    ex_range = forms.IntegerField(required=False, help_text='Bandwidth range for Ex max search')
    em_max = forms.IntegerField(required=False, help_text='Search for protein with Em max around this wavelength')
    em_range = forms.IntegerField(required=False, help_text='Bandwidth range for Em max search')
    name = forms.CharField(required=False)

    class Meta:
        model = Protein
        fields = ['name', 'gb_prot', 'gb_nuc', 'agg', 'switch_type', 'parent_organism', 'FRET_partner']
        help_texts = {
            'name': _('Query string (will search within names)'),
            'gb_prot': _('GenBank protein Accession number (e.g. AFR60231)'),
            'gb_nuc': _('GenBank nucleotide Accession number'),
            'switch_type': _('Photoswitching type'),
        }
        widgets = {
            'name': forms.TextInput(attrs={'class': 'myfieldclass'}),
        }

    def clean(self):
        cleaned_data = super(ProteinSearchForm, self).clean()
        form_empty = True
        for field_value in cleaned_data.values():
            # Check for None or '', so IntegerFields with 0 or similar things don't seem empty.
            if not (field_value is None or field_value == ''):
                if not field_value:
                    continue
                form_empty = False
                break
        if form_empty:
            raise forms.ValidationError(_("You must fill at least one field!"))
        return cleaned_data   # Important that clean should return cleaned_data!

    helper = FormHelper()
    helper.form_method = 'post'
    helper.form_action = '/search/'
    helper.layout = Layout(
        Fieldset(
            'Search',
            'name',
        ),
        Accordion(
            AccordionGroup(
                'Advanced Search Parameters...',
                Div(
                    Div('gb_prot',
                        css_class='col-sm-6',
                        ),
                    Div('gb_nuc',
                        css_class='col-sm-6',
                        ),
                    css_class='row',
                ),
                Div(
                    Div('switch_type',
                        css_class='col-sm-4',
                        ),
                    Div('agg',
                        css_class='col-sm-4',
                        ),
                    Div('parent_organism',
                        css_class='col-sm-4',
                        ),
                    css_class='row',
                ),
                Div(
                    Div('ex_max',
                        css_class='col-sm-4',
                        ),
                    Div('ex_range',
                        css_class='col-sm-2',
                        ),
                    Div('em_max',
                        css_class='col-sm-4',
                        ),
                    Div('em_range',
                        css_class='col-sm-2',
                        ),
                    css_class='row',
                ),
                'FRET_partner',
                active=False,
            ),
            ),
        ButtonHolder(
            Submit('submit', 'Submit', css_class='button white')
        ),
    )


class ProteinForm(forms.ModelForm):
    class Meta:
        model = Protein
        # fields = ['name', 'slug', 'gb_prot', 'gb_nuc', 'seq', 'parent_organism', 'primary_reference', 'references', 'FRET_partner', 'added_by', 'updated_by']
        fields = ['name', 'slug', 'gb_prot', 'gb_nuc', 'seq', 'parent_organism', 'FRET_partner', 'added_by', 'updated_by']

    def clean_added_by(self):
        if not self.cleaned_data['added_by']:
            return User()
        return self.cleaned_data['added_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']


class SpectrumFormField(forms.CharField):
    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={'class': 'vLargeTextField'})


class StateForm(forms.ModelForm):
    ex_spectra = SpectrumFormField(required=False)
    em_spectra = SpectrumFormField(required=False)

    class Meta:
        model = State
        fields = ['protein', 'is_dark', 'name', 'ex_max', 'em_max', 'ex_spectra', 'em_spectra', 'ext_coeff', 'qy', 'pka', 'maturation', 'lifetime', 'added_by', 'updated_by']

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

    def clean_added_by(self):
        if not self.cleaned_data['added_by']:
            return User()
        return self.cleaned_data['added_by']

    def clean_updated_by(self):
        if not self.cleaned_data['updated_by']:
            return User()
        return self.cleaned_data['updated_by']


# class StateFormSet(BaseInlineFormSet):
#     def clean(self):
#         super(StateFormSet, self).clean()
#         if any(self.errors):
#             # Don't bother validating the formset unless each form is valid on its own
#             return
#         defaults = 0
#         for form in self.forms:
#             if form.cleaned_data['default']:
#                 defaults += 1
#             if defaults > 1:
#                 raise forms.ValidationError("Only one default state is allowed.")


# notes:

# it is strongly recommended that you explicitly set all fields that should be edited in the form using the fields attribute
# Failure to do so can easily lead to security problems when a form unexpectedly allows a user to set certain fields
# https://docs.djangoproject.com/en/1.7/topics/forms/modelforms/#model-formsets
