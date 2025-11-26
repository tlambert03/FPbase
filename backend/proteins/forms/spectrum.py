from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Field, Layout, Submit
from dal import autocomplete
from django import forms
from django.apps import apps
from django.core.exceptions import ObjectDoesNotExist
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.utils.text import slugify

from proteins.models import Dye, DyeState, FluorState, Spectrum, State
from proteins.util.helpers import zip_wave_data
from proteins.util.importers import text_to_spectra
from proteins.validators import validate_spectrum


class SpectrumFormField(forms.CharField):
    default_validators = [validate_spectrum]
    widget = forms.Textarea(attrs={"class": "vLargeTextField", "rows": 2})

    def __init__(self, *args, **kwargs):
        if "help_text" not in kwargs:
            kwargs["help_text"] = (
                "List of [wavelength, value] pairs, e.g. [[300, 0.5], [301, 0.6],... ]. "
                "File data takes precedence."
            )
        super().__init__(*args, **kwargs)


class SpectrumForm(forms.ModelForm):
    # Lookup for non-protein, non-dye categories (filter/camera/light)
    # Proteins use owner_fluor Select2 field directly
    # Dyes are handled specially in save() and clean_owner()
    lookup = {
        Spectrum.FILTER: ("owner_filter", "Filter"),
        Spectrum.CAMERA: ("owner_camera", "Camera"),
        Spectrum.LIGHT: ("owner_light", "Light"),
    }

    owner_fluor = forms.ModelChoiceField(
        required=False,
        label=mark_safe('Protein<span class="asteriskField">*</span>'),
        queryset=State.objects.select_related("protein"),
        widget=autocomplete.ModelSelect2(
            url="proteins:state-autocomplete",
            attrs={"data-theme": "bootstrap-5", "data-width": "100%"},
        ),
    )
    owner = forms.CharField(
        max_length=100,
        label=mark_safe(
            '<span class="owner-type">Owner</span> Name<span class="asteriskField">*</span>'
        ),
        required=False,
        help_text="Name of protein, dye, filter, etc...",
    )
    # FIXME!!
    # this is an actual problem.  this field needs to be renamed to avoid name with the Form
    data = SpectrumFormField(required=False, label="Data")  # pyright: ignore[reportAssignmentType]
    file = forms.FileField(
        required=False,
        label="File Upload",
        help_text=(
            "2 column CSV/TSV file with wavelengths in first column and data in second column"
        ),
    )
    confirmation = forms.BooleanField(
        required=True,
        label=mark_safe(
            "<span class='small'>I understand that I am adding a spectrum to the <em>public</em> "
            "FPbase spectra database, and confirm that I have verified the validity of the "
            "data</span>"
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
                Div("owner_fluor", css_class="col-sm-6 col-xs-12 protein-owner hidden"),
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
            "owner_fluor",
            "owner",
        )
        widgets = {"data": forms.Textarea(attrs={"class": "vLargeTextField", "rows": 2})}

    def clean(self):
        cleaned_data = super().clean()

        # Check which data source was selected based on form submission
        data_source = self.data.get("data_source", "file") if self.data else "file"

        # Validate based on the selected data source
        if data_source == "manual":
            # Manual data tab: require manual data, file is optional
            if not cleaned_data.get("data"):
                raise forms.ValidationError("Please enter valid spectrum data.")
        else:
            # File tab: require file upload, manual data is optional
            if not self.files.get("file"):
                raise forms.ValidationError("Please select a file to upload.")

        return cleaned_data

    def save(self, commit=True):
        # Set spectrum data from form - data is now a property, not a model field,
        # so ModelForm won't automatically set it on the instance
        if self.cleaned_data.get("data"):
            self.instance.data = self.cleaned_data["data"]

        cat = self.cleaned_data.get("category")
        if cat == Spectrum.DYE:
            # Dyes require special handling: create Dye first, then DyeState
            # This mirrors protein behavior where State belongs to Protein
            # TODO: unify dye/protein behavior - currently users can create dyes
            # but not proteins. Consider using Select2 for dyes too.
            owner_name = self.cleaned_data.get("owner")
            dye, created = Dye.objects.get_or_create(
                slug=slugify(owner_name),
                defaults={"name": owner_name, "created_by": self.user},
            )
            if not created:
                dye.updated_by = self.user
                dye.save()
            # Get or create the default DyeState for this dye
            dye_state, _ = DyeState.objects.get_or_create(
                dye=dye,
                name=FluorState.DEFAULT_NAME,
                defaults={"created_by": self.user},
            )
            self.instance.owner_fluor = dye_state
        elif cat != Spectrum.PROTEIN:
            # Filter, Camera, Light - create owner directly
            owner_model = apps.get_model("proteins", self.lookup[cat][1])
            owner_name = self.cleaned_data.get("owner")
            owner_obj, created = owner_model.objects.get_or_create(
                name=owner_name, defaults={"created_by": self.user}
            )
            if not created:
                owner_obj.updated_by = self.user
                owner_obj.save()
            setattr(self.instance, self.lookup[cat][0], owner_obj)
        self.instance.created_by = self.user
        return super().save(commit=commit)

    def clean_file(self):
        if self.files:
            filetext = ""
            x = []
            y = []
            try:
                for chunk in self.files["file"].chunks():
                    try:
                        filetext += chunk.decode("utf-8")
                    except AttributeError:
                        filetext += chunk
                x, y, _headers = text_to_spectra(filetext)
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
                self.data: dict = self.data.copy()
                self.data["data"] = self.cleaned_data["data"]

    def clean_owner_fluor(self):
        owner_fluor = self.cleaned_data.get("owner_fluor")
        stype = self.cleaned_data.get("subtype")
        # Only proteins use the owner_fluor Select2 field
        # Dyes use the owner text field and are validated in clean_owner()
        if self.cleaned_data.get("category") == Spectrum.PROTEIN and owner_fluor:
            spectra = Spectrum.objects.all_objects().filter(owner_fluor=owner_fluor, subtype=stype)
            if spectra.exists():
                first = spectra.first()
                self.add_error(
                    "owner_fluor",
                    forms.ValidationError(
                        "%(owner)s already has a%(n)s %(stype)s spectrum %(status)s",
                        params={
                            "owner": owner_fluor,
                            "stype": first.get_subtype_display().lower(),
                            "n": "n" if stype != Spectrum.TWOP else "",
                            "status": " (pending)"
                            if first.status == Spectrum.STATUS.pending
                            else "",
                        },
                        code="owner_exists",
                    ),
                )
        return owner_fluor

    def clean_owner(self):
        # make sure an owner with the same category and name doesn't already exist
        owner = self.cleaned_data.get("owner")
        cat = self.cleaned_data.get("category")
        stype = self.cleaned_data.get("subtype")
        if cat == Spectrum.PROTEIN:
            return owner

        # Dyes need special handling: look up Dye, then check its default DyeState
        if cat == Spectrum.DYE:
            try:
                dye = Dye.objects.get(slug=slugify(owner))
            except Dye.DoesNotExist:
                return owner
            dye_state = dye.states.filter(name=FluorState.DEFAULT_NAME).first()
            if dye_state and dye_state.spectra.filter(subtype=stype).exists():
                stype_display = (
                    dye_state.spectra.filter(subtype=stype).first().get_subtype_display()
                )
                self.add_error(
                    "owner",
                    forms.ValidationError(
                        "A dye with the name %(name)s already has a %(stype)s spectrum",
                        params={"name": dye.name, "stype": stype_display},
                        code="owner_exists",
                    ),
                )
            return owner

        # Filter, Camera, Light - look up by slug directly
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
            raise forms.ValidationError("Category not recognized") from e

        # object exists... check if it has this type of spectrum
        if hasattr(obj, "spectrum") and obj.spectrum:
            self.add_error(
                "owner",
                forms.ValidationError(
                    "A %(model)s with the name %(name)s already has a %(stype)s spectrum",
                    params={
                        "model": self.lookup[cat][1].lower(),
                        "name": obj.name,
                        "stype": obj.spectrum.get_subtype_display(),
                    },
                    code="owner_exists",
                ),
            )
        return owner
