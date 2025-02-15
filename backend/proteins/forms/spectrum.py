from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Field, Layout, Submit
from dal import autocomplete
from django import forms
from django.apps import apps
from django.core.exceptions import ObjectDoesNotExist
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.utils.text import slugify

from proteins.models import Fluorophore, Spectrum, State
from proteins.util.helpers import zip_wave_data
from proteins.util.importers import text_to_spectra
from proteins.validators import validate_spectrum


class SpectrumFormField(forms.CharField):
    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={"class": "vLargeTextField", "rows": 2})

    def __init__(self, *args, **kwargs):
        if "help_text" not in kwargs:
            kwargs["help_text"] = (
                "List of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]. File data takes precedence."
            )
        super().__init__(*args, **kwargs)


class SpectrumForm(forms.ModelForm):
    lookup = {
        Spectrum.DYE: ("owner_dye", "Dye"),
        Spectrum.PROTEIN: ("owner_state", "State"),
        Spectrum.FILTER: ("owner_filter", "Filter"),
        Spectrum.CAMERA: ("owner_camera", "Camera"),
        Spectrum.LIGHT: ("owner_light", "Light"),
    }

    owner_state = forms.ModelChoiceField(
        required=False,
        label=mark_safe('Protein<span class="asteriskField">*</span>'),
        queryset=State.objects.select_related("protein"),
        widget=autocomplete.ModelSelect2(
            url="proteins:state-autocomplete",
            attrs={"data-theme": "bootstrap", "data-width": "100%"},
        ),
    )
    owner = forms.CharField(
        max_length=100,
        label=mark_safe('<span class="owner-type">Owner</span> Name<span class="asteriskField">*</span>'),
        required=False,
        help_text="Name of protein, dye, filter, etc...",
    )
    data = SpectrumFormField(required=False, label="Data")
    file = forms.FileField(
        required=False,
        label="File Upload",
        help_text="2 column CSV/TSV file with wavelengths in first column and data in second column",
    )
    confirmation = forms.BooleanField(
        required=True,
        label=mark_safe(
            "<span class='small'>I understand that I am adding a spectrum to the <em>public</em> "
            "FPbase spectra database, and confirm that I have verified the validity of the data</span>"
        ),
    )

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user", None)
        self.helper = FormHelper()
        self.helper.attrs = {
            "id": "spectrum-submit-form",
            "data-validate-owner-url": reverse("proteins:validate_spectrumownername"),
        }
        self.helper.add_input(Submit("submit", "Submit"))
        self.helper.layout = Layout(
            Div("category"),
            Div(
                Div("owner_state", css_class="col-sm-6 col-xs-12 protein-owner hidden"),
                Div("owner", css_class="col-sm-6 col-xs-12 non-protein-owner"),
                Div("subtype", css_class="col-sm-6 col-xs-12"),
                css_class="row",
            ),
            Div("file", "data"),
            Div(
                Div("ph", css_class="col-md-6 col-sm-12"),
                Div("solvent", css_class="col-md-6 col-sm-12"),
                css_class="row",
            ),
            Field("confirmation", css_class="custom-checkbox"),
        )

        super().__init__(*args, **kwargs)

    class Meta:
        model = Spectrum
        fields = (
            "category",
            "subtype",
            "file",
            "data",
            "ph",
            "source",
            "solvent",
            "owner_state",
            "owner",
        )
        widgets = {"data": forms.Textarea(attrs={"class": "vLargeTextField", "rows": 2})}

    def clean(self):
        cleaned_data = super().clean()
        if not (cleaned_data.get("data") or self.files):
            self.add_error(
                "data",
                "Please either fill in the data field or select a file to upload.",
            )
            self.add_error(
                "file",
                "Please either fill in the data field or select a file to upload.",
            )

    def save(self, commit=True):
        cat = self.cleaned_data.get("category")
        if cat != Spectrum.PROTEIN:
            owner_model = apps.get_model("proteins", self.lookup[cat][1])
            owner_name = self.cleaned_data.get("owner")
            owner_obj, c = owner_model.objects.get_or_create(name=owner_name, defaults={"created_by": self.user})
            if not c:
                owner_obj.update_by = self.user
                owner_obj.save()
            setattr(self.instance, self.lookup[cat][0], owner_obj)
        self.instance.created_by = self.user
        return super().save(commit=commit)

    def clean_file(self):
        if self.files:
            filetext = ""
            try:
                for chunk in self.files["file"].chunks():
                    try:
                        filetext += chunk.decode("utf-8")
                    except AttributeError:
                        filetext += chunk
                x, y, headers = text_to_spectra(filetext)
                if not len(y):
                    self.add_error("file", "Did not find a data column in the provided file")
                if not len(x):
                    self.add_error("file", "Could not parse wavelengths from first column")
            except Exception:
                self.add_error(
                    "file",
                    "Sorry, could not parse spectrum from this file. "
                    "Is it it two column csv with (wavelength, spectrum)?",
                )
            if not self.errors:
                self.cleaned_data["data"] = zip_wave_data(x, y[0])
                self.data = self.data.copy()
                self.data["data"] = self.cleaned_data["data"]

    def clean_owner_state(self):
        owner_state = self.cleaned_data.get("owner_state")
        stype = self.cleaned_data.get("subtype")
        if self.cleaned_data.get("category") == Spectrum.PROTEIN:
            spectra = Spectrum.objects.all_objects().filter(owner_state=owner_state, subtype=stype)
            if spectra.exists():
                first = spectra.first()
                self.add_error(
                    "owner_state",
                    forms.ValidationError(
                        "%(owner)s already has a%(n)s %(stype)s spectrum %(status)s",
                        params={
                            "owner": owner_state,
                            "stype": first.get_subtype_display().lower(),
                            "n": "n" if stype != Spectrum.TWOP else "",
                            "status": " (pending)" if first.status == Spectrum.STATUS.pending else "",
                        },
                        code="owner_exists",
                    ),
                )
        return owner_state

    def clean_owner(self):
        # make sure an owner with the same category and name doesn't already exist
        owner = self.cleaned_data.get("owner")
        cat = self.cleaned_data.get("category")
        stype = self.cleaned_data.get("subtype")
        if cat == Spectrum.PROTEIN:
            return owner

        try:
            mod = apps.get_model("proteins", self.lookup[cat][1])
            obj = mod.objects.get(slug=slugify(owner))
        except ObjectDoesNotExist:
            return owner
        except KeyError as e:
            # this might be repetitive... since a missing category will already
            # throw an error prior to this point
            if not cat:
                raise forms.ValidationError("Category not provided") from e
            else:
                raise forms.ValidationError("Category not recognized") from e
        else:
            # object exists... check if it has this type of spectrum
            exists = False
            if isinstance(obj, Fluorophore):
                if obj.spectra.filter(subtype=stype).exists():
                    exists = True
                    stype = obj.spectra.filter(subtype=stype).first().get_subtype_display()
            elif hasattr(obj, "spectrum") and obj.spectrum:
                exists = True
                stype = obj.spectrum.get_subtype_display()
            if exists:
                self.add_error(
                    "owner",
                    forms.ValidationError(
                        "A %(model)s with the name %(name)s already has a %(stype)s spectrum",
                        params={
                            "model": self.lookup[cat][1].lower(),
                            "name": obj.name,
                            "stype": stype,
                        },
                        code="owner_exists",
                    ),
                )
                # raise forms.ValidationError(
                #     "A %(model)s with the slug %(slug)s already has a spectrum of type %(stype)s.",
                #     params={'model': self.lookup[cat][1].lower(), 'slug': slugify(owner), 'stype': stype},
                #     code='owner_exists')
            else:
                return owner
