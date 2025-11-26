import io

from django import forms
from django.contrib import admin
from django.db.models import Count, Prefetch
from django.forms import TextInput
from django.urls import reverse
from django.utils.safestring import mark_safe
from mptt.admin import MPTTModelAdmin
from reversion.admin import VersionAdmin
from reversion_compare.admin import CompareVersionAdmin

from fpbase.util import uncache_protein_page
from proteins.models import (
    BleachMeasurement,
    Camera,
    Dye,
    DyeState,
    Excerpt,
    Filter,
    FilterPlacement,
    FluorState,
    Light,
    Lineage,
    Microscope,
    OpticalConfig,
    Organism,
    OSERMeasurement,
    Protein,
    ProteinCollection,
    SnapGenePlasmid,
    Spectrum,
    State,
    StateTransition,
)
from proteins.models.lineage import MutationSetField
from proteins.util.maintain import validate_node

# ############ INLINES ###############


class SpectrumOwner:
    list_display = ("__str__", "spectra", "created_by", "created")
    list_select_related = ("created_by",)
    list_filter = ("created",)
    search_fields = ("name",)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.readonly_fields = (*self.readonly_fields, "spectra")

    @admin.display(description="spectra")
    def spectra(self, obj):
        def _makelink(sp):
            url = reverse("admin:proteins_spectrum_change", args=(sp.pk,))
            pending = " (pending)" if sp.status == Spectrum.STATUS.pending else ""
            return f'<a href="{url}">{sp.get_subtype_display()}{pending}</a>'

        links = []
        if isinstance(obj, FluorState):
            [links.append(_makelink(sp)) for sp in obj.spectra.all()]
        else:
            links.append(_makelink(obj.spectrum))
        return mark_safe(", ".join(links))

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related(Prefetch("spectrum", queryset=Spectrum.objects.all_objects()))


class MultipleSpectraOwner(SpectrumOwner):
    def get_queryset(self, request):
        qs = super(SpectrumOwner, self).get_queryset(request)
        return qs.prefetch_related(Prefetch("spectra", queryset=Spectrum.objects.all_objects()))


class BleachInline(admin.TabularInline):
    model = BleachMeasurement
    autocomplete_fields = ("reference",)
    extra = 1


class OSERInline(admin.StackedInline):
    model = OSERMeasurement
    autocomplete_fields = ("reference",)
    extra = 1
    readonly_fields = ("created", "created_by", "modified", "updated_by")
    fieldsets = [
        (None, {"fields": (("reference", "celltype"),)}),
        (None, {"fields": (("percent", "percent_stddev", "percent_ncells"),)}),
        (None, {"fields": (("oserne", "oserne_stddev", "oserne_ncells"),)}),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": (("created", "created_by"), ("modified", "updated_by")),
            },
        ),
    ]


class FluorStateInline(MultipleSpectraOwner):
    model = FluorState
    extra = 0
    can_delete = True
    show_change_link = True
    fieldsets = [
        (None, {"fields": (("name", "slug", "is_dark"),)}),
        (
            None,
            {
                "fields": (
                    ("ex_max", "em_max"),
                    ("ext_coeff", "qy"),
                    ("twop_ex_max", "twop_peak_gm", "twop_qy"),
                    "pka",
                    "lifetime",
                    "spectra",
                )
            },
        ),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": (("created", "created_by"), ("modified", "updated_by")),
            },
        ),
    ]
    readonly_fields = (
        "slug",
        "created",
        "created_by",
        "modified",
        "updated_by",
    )


class StateInline(FluorStateInline, admin.StackedInline):
    # form = StateForm
    # formset = StateFormSet
    model = State
    fieldsets = [
        *FluorStateInline.fieldsets,
        (None, {"fields": ("bleach_links",)}),
        (None, {"fields": ("maturation",)}),
    ]

    readonly_fields = (*FluorStateInline.readonly_fields, "bleach_links")

    @admin.display(description="BleachMeasurements")
    def bleach_links(self, obj):
        links = []
        for bm in obj.bleach_measurements.all():
            url = reverse("admin:proteins_bleachmeasurement_change", args=(bm.pk,))
            link = f'<a href="{url}">{bm}</a>'
            links.append(link)
        return mark_safe(", ".join(links))


class DyeStateInline(FluorStateInline, admin.StackedInline):
    model = DyeState


class LineageInline(admin.TabularInline):
    model = Lineage
    autocomplete_fields = ("parent", "root_node")
    readonly_fields = ("root_node", "created_by", "created", "updated_by")
    fields = ("parent", "mutation", "created_by", "created")


# ############# MODELS ###############


@admin.register(Light)
class LightAdmin(SpectrumOwner, admin.ModelAdmin):
    model = Light
    ordering = ("-created",)


@admin.register(Dye)
class DyeAdmin(admin.ModelAdmin):
    model = Dye
    list_display = ("__str__", "created_by", "created")
    ordering = ("-created",)
    list_filter = ("created", "manufacturer")
    fields = (
        "name",
        "slug",
        "manufacturer",
        "default_state",
        "primary_reference",
        "created",
        "modified",
        "created_by",
        "updated_by",
    )
    readonly_fields = ("created", "modified")
    inlines = (DyeStateInline,)


@admin.register(Filter)
class FilterAdmin(SpectrumOwner, VersionAdmin):
    model = Filter
    ordering = ("-created",)
    readonly_fields = ("configs",)

    @admin.display(description="OC Memberships")
    def configs(self, obj):
        def _makelink(oc):
            url = reverse("admin:proteins_opticalconfig_change", args=(oc.pk,))
            return f'<a href="{url}">{oc}</a>'

        links = []
        [links.append(_makelink(oc)) for oc in obj.optical_configs.all()]
        return mark_safe(", ".join(links))


@admin.register(Camera)
class CameraAdmin(SpectrumOwner, VersionAdmin):
    model = Camera
    ordering = ("-created",)


@admin.register(Spectrum)
class SpectrumAdmin(VersionAdmin):
    model = Spectrum
    autocomplete_fields = ["reference"]
    list_select_related = (
        "owner_fluor",
        "owner_filter",
        "owner_camera",
        "owner_light",
        "created_by",
    )
    list_display = ("__str__", "category", "subtype", "owner", "created_by")
    list_filter = ("status", "created", "category", "subtype")
    readonly_fields = ("owner", "name", "created", "modified", "spectrum_preview")
    search_fields = (
        "owner_fluor__name",
        "owner_filter__name",
        "owner_camera__name",
        "owner_light__name",
    )

    def get_fields(self, request, obj=None):
        fields = []
        if not obj or not obj.category:
            # If no category yet, allow selecting any owner type
            own = [
                "owner_fluor",
                "owner_filter",
                "owner_camera",
                "owner_light",
            ]
        elif obj.category in (Spectrum.PROTEIN, Spectrum.DYE):
            # Protein and Dye both use owner_fluor (Fluorophore)
            own = ["owner_fluor"]
        else:
            # Filter, Camera, Light
            own = ["owner_" + obj.get_category_display().split(" ")[0].lower()]
        fields.extend(own)
        # Add clickable link to owner's admin page (for existing objects)
        if obj and obj.owner:
            fields.append("owner")
        self.autocomplete_fields.extend(own)
        fields += [
            "category",
            "subtype",
            "data",
            "spectrum_preview",
            "ph",
            "solvent",
            "source",
            "reference",
            "status",
            ("created", "created_by"),
            ("modified", "updated_by"),
        ]
        return fields

    @admin.display(description="Owner")
    def owner(self, obj):
        owner = obj.owner
        # FluorState is a base class - resolve to the actual subclass admin
        if isinstance(owner, FluorState):
            model_name = "state" if owner.entity_type == FluorState.EntityTypes.PROTEIN else "dyestate"
        else:
            model_name = owner._meta.model.__name__.lower()

        url = reverse(f"admin:proteins_{model_name}_change", args=(owner.pk,))
        link = f'<a href="{url}">{owner}</a>'
        return mark_safe(link)

    @admin.display(description="Spectrum Preview")
    def spectrum_preview(self, obj: Spectrum) -> str:
        """Show a visual preview of the spectrum for pending spectra."""
        if not obj.data:
            return "N/A"

        try:
            # Generate the spectrum image as SVG
            output = io.BytesIO()
            obj.spectrum_img(fmt="svg", output=output, xlabels=True, ylabels=False, figsize=(9, 2))
            output.seek(0)
            svg_data = output.read().decode("utf-8")

            # Return the SVG directly embedded in the admin
            return mark_safe(f'<div style="max-width: 1200px; overflow: auto;">{svg_data}</div>')
        except Exception as e:
            return f"Error generating preview: {e!s}"

    def get_queryset(self, request):
        """
        Return a QuerySet of all model instances that can be edited by the
        admin site. This is used by changelist_view.
        """
        qs = Spectrum.objects.all_objects()
        ordering = self.get_ordering(request)
        if ordering:
            qs = qs.order_by(*ordering)
        return qs


@admin.register(OSERMeasurement)
class OSERMeasurementAdmin(VersionAdmin):
    model = OSERMeasurement
    autocomplete_fields = ("protein", "reference")
    list_select_related = ("protein", "reference")
    list_display = ("protein", "percent", "reference", "celltype")
    list_filter = ("celltype",)
    search_fields = ("protein__name",)
    readonly_fields = ("created", "created_by", "modified", "updated_by")
    orders = ("reference", "protein")
    fieldsets = [
        (None, {"fields": (("protein", "reference", "celltype"),)}),
        (None, {"fields": (("percent", "percent_stddev", "percent_ncells"),)}),
        (None, {"fields": (("oserne", "oserne_stddev", "oserne_ncells"),)}),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": (("created", "created_by"), ("modified", "updated_by")),
            },
        ),
    ]


@admin.register(BleachMeasurement)
class BleachMeasurementAdmin(VersionAdmin):
    model = BleachMeasurement
    autocomplete_fields = ("state", "reference")
    list_select_related = ("state", "state__protein")


@admin.register(State)
class StateAdmin(CompareVersionAdmin):
    # form = StateForm
    model = State
    list_select_related = ("protein", "created_by", "updated_by")
    search_fields = ("protein__name",)
    list_display = (
        "__str__",
        "protein_link",
        "ex_max",
        "em_max",
        "created_by",
        "updated_by",
        "modified",
    )
    list_filter = ("created", "modified")
    inlines = (BleachInline,)
    fieldsets = [
        (None, {"fields": (("name", "slug", "is_dark"),)}),
        (
            None,
            {
                "fields": (
                    ("ex_max", "em_max"),
                    ("ext_coeff", "qy"),
                    ("twop_ex_max", "twop_peak_gm", "twop_qy"),
                    ("pka", "maturation"),
                    "lifetime",
                )
            },
        ),
    ]

    @admin.display(description="Protein")
    def protein_link(self, obj):
        url = reverse("admin:proteins_protein_change", args=([obj.protein.pk]))
        return mark_safe(f'<a href="{url}">{obj.protein}</a>')


@admin.register(DyeState)
class DyeStateAdmin(MultipleSpectraOwner, CompareVersionAdmin):
    # form = StateForm
    model = State
    list_select_related = ("dye", "created_by", "updated_by")
    search_fields = ("dye__name",)
    list_display = (
        "__str__",
        "dye_link",
        "ex_max",
        "em_max",
        "created_by",
        "updated_by",
        "modified",
    )
    list_filter = ("created", "modified")
    fieldsets = [
        (None, {"fields": (("name", "slug", "is_dark"),)}),
        (
            None,
            {
                "fields": (
                    ("ex_max", "em_max"),
                    ("ext_coeff", "qy"),
                    ("twop_ex_max", "twop_peak_gm", "twop_qy"),
                    ("pka", "lifetime"),
                )
            },
        ),
        (
            None,
            {
                "fields": ("spectra",),
            },
        ),
    ]

    @admin.display(description="Dye")
    def dye_link(self, obj):
        url = reverse("admin:proteins_dye_change", args=([obj.dye.pk]))
        return mark_safe(f'<a href="{url}">{obj.dye}</a>')


class StateTransitionAdmin(VersionAdmin):
    model = StateTransition
    list_select_related = ("protein", "from_state", "to_state")
    autocomplete_fields = ("protein", "from_state", "to_state")


class StateTransitionInline(admin.TabularInline):
    model = StateTransition
    extra = 1

    def formfield_for_foreignkey(self, db_field, request=None, **kwargs):
        field = super().formfield_for_foreignkey(db_field, request, **kwargs)

        # restrict
        if db_field.related_model == State:
            obj = getattr(request, "object", None)
            if request is not None:
                field.queryset = field.queryset.filter(protein=obj)
            else:
                field.queryset = field.queryset.none()

        return field


@admin.register(Organism)
class OrganismAdmin(CompareVersionAdmin):
    list_select_related = ("created_by", "updated_by")
    list_display = (
        "scientific_name",
        "id",
        "created",
        "created_by",
        "modified",
        "updated_by",
    )
    list_filter = ("created", "modified")
    search_fields = (
        "scientific_name",
        "common_name",
        "id",
        "created_by__username",
        "created_by__first_name",
        "created_by__last_name",
    )

    fieldsets = [
        (
            "Organism",
            {
                "fields": (
                    "id",
                    "scientific_name",
                    "common_name",
                    "rank",
                    "division",
                    "genus",
                    "species",
                )
            },
        ),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": ("created", "created_by", "modified", "updated_by"),
            },
        ),
    ]
    readonly_fields = (
        "scientific_name",
        "common_name",
        "division",
        "rank",
        "genus",
        "species",
        "created",
        "created_by",
        "modified",
        "updated_by",
    )

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()


@admin.action(description="Mark selected proteins as approved")
def approve_protein(modeladmin, request, queryset):
    # note, this will fail if the list is ordered by numproteins
    queryset.update(status=Protein.STATUS.approved)


@admin.register(Protein)
class ProteinAdmin(CompareVersionAdmin):
    autocomplete_fields = ("parent_organism", "references", "primary_reference")
    filter_horizontal = ("snapgene_plasmids",)
    list_display = (
        "__str__",
        "created",
        "modified",
        "created_by",
        "updated_by",
        "nstates",
    )
    list_filter = ("status", "created", "modified", "switch_type")
    list_select_related = ("created_by", "updated_by")
    search_fields = (
        "name",
        "aliases",
        "slug",
        "ipg_id",
        "created_by__username",
        "created_by__first_name",
        "created_by__last_name",
    )
    prepopulated_fields = {"slug": ("name",)}
    inlines = (StateInline, StateTransitionInline, OSERInline, LineageInline)
    actions = [approve_protein]
    fieldsets = [
        (
            None,
            {
                "fields": (
                    ("uuid", "name", "slug"),
                    ("aliases", "status", "chromophore"),
                    ("seq", "seq_validated", "seq_comment"),
                    ("ipg_id", "genbank", "uniprot", "pdb"),
                    ("parent_organism", "switch_type"),
                    ("agg", "mw"),
                    "blurb",
                )
            },
        ),
        (None, {"fields": (("primary_reference", "references"),)}),
        (None, {"fields": ("snapgene_plasmids",)}),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": (("created", "created_by"), ("modified", "updated_by")),
            },
        ),
    ]
    readonly_fields = (
        "created",
        "created_by",
        "modified",
        "updated_by",
        "switch_type",
        "uuid",
    )

    def nstates(self, obj):
        if obj:
            return obj.states.all().count()
        else:
            return 0

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related("states")

    def get_form(self, request, obj=None, **kwargs):
        request.object = obj
        return super().get_form(request, obj, **kwargs)

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        # obj.status = 'approved'
        obj.save()
        uncache_protein_page(obj.slug, request)

    def save_formset(self, request, form, formset, change):
        instances = formset.save(commit=False)
        for instance in instances:
            if not instance.created_by:
                instance.created_by = request.user
            instance.updated_by = request.user
            instance.save()
        formset.save()


class FilterPlacementInline(admin.TabularInline):
    model = FilterPlacement
    extra = 1  # how many rows to show
    autocomplete_fields = ("filter",)


@admin.register(OpticalConfig)
class OpticalConfigAdmin(admin.ModelAdmin):
    model = OpticalConfig
    inlines = (FilterPlacementInline,)
    list_display = ("__str__", "microscope", "owner_link", "created")
    fields = (("name", "microscope"), ("laser", "light"), ("camera", "owner"))

    @admin.display(description="Owner")
    def owner_link(self, obj):
        if obj.microscope and obj.microscope.owner:
            url = reverse("admin:users_user_change", args=([obj.microscope.owner.pk]))
            return mark_safe(f'<a href="{url}">{obj.microscope.owner}</a>')

    def save_model(self, request, obj, form, change):
        obj.save()
        obj.microscope.save()

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related("microscope__owner")


@admin.register(Microscope)
class MicroscopeAdmin(admin.ModelAdmin):
    model = Microscope
    autocomplete_fields = ("extra_lights", "extra_cameras")
    readonly_fields = ("configs",)
    list_display = ("__str__", "owner_link", "created", "modified", "OCs")
    list_select_related = ("owner",)
    list_filter = ("created",)
    ordering = ("-modified",)

    @admin.display(description="Owner")
    def owner_link(self, obj):
        if obj.owner:
            url = reverse("admin:users_user_change", args=([obj.owner.pk]))
            return mark_safe(f'<a href="{url}">{obj.owner}</a>')

    @admin.display(ordering="oc_count")
    def OCs(self, obj):
        return obj.optical_configs.count()

    @admin.display(description="Optical Configs")
    def configs(self, obj):
        def _makelink(oc):
            url = reverse("admin:proteins_opticalconfig_change", args=(oc.pk,))
            return f'<a href="{url}">{oc}</a>'

        links = []
        [links.append(_makelink(oc)) for oc in obj.optical_configs.all()]
        return mark_safe(", ".join(links))

    def get_queryset(self, request):
        qs = super().get_queryset(request).annotate(oc_count=Count("optical_configs"))
        return qs.prefetch_related("optical_configs")


@admin.action(description="Mark selected collections as private")
def make_private(modeladmin, request, queryset):
    # note, this will fail if the list is ordered by numproteins
    queryset.update(private=True)


@admin.register(ProteinCollection)
class ProteinCollectionAdmin(admin.ModelAdmin):
    model = ProteinCollection
    list_display = ("__str__", "owner_link", "created", "private", "numproteins")
    list_filter = ("created", "private")
    readonly_fields = ("numproteins", "owner_link")
    list_select_related = ("owner",)
    autocomplete_fields = ("proteins",)
    search_fields = ("name", "proteins__name")
    actions = [make_private]

    @admin.display(description="Owner")
    def owner_link(self, obj):
        url = reverse("admin:users_user_change", args=([obj.owner.pk]))
        return mark_safe(f'<a href="{url}">{obj.owner}</a>')

    @admin.display(ordering="proteins_count")
    def numproteins(self, obj):
        return obj.proteins.count()

    def get_queryset(self, request):
        qs = super().get_queryset(request).annotate(proteins_count=Count("proteins"))
        return qs.prefetch_related("proteins")


class LineageAdminForm(forms.ModelForm):
    class Meta:
        model = Lineage
        fields = ("protein", "parent", "mutation", "root_node", "rootmut")
        readonly_fields = "root_node"

    parent = forms.ModelChoiceField(required=False, queryset=Lineage.objects.prefetch_related("protein").all())
    root_node = forms.ModelChoiceField(queryset=Lineage.objects.prefetch_related("protein").all())


@admin.register(Excerpt)
class ExcerptAdmin(VersionAdmin):
    model = Excerpt
    list_display = ("id", "reference", "content", "created_by", "created")
    list_filter = ("status",)
    list_select_related = ("reference", "created_by")
    autocomplete_fields = ("proteins", "reference", "created_by", "updated_by")


@admin.register(SnapGenePlasmid)
class SnapGenePlasmidAdmin(admin.ModelAdmin):
    model = SnapGenePlasmid
    list_display = ("name", "plasmid_id", "author", "size", "protein_count")
    search_fields = ("name", "plasmid_id", "description")
    list_filter = ("topology", "author")
    readonly_fields = ("protein_count",)

    @admin.display(description="# Proteins")
    def protein_count(self, obj):
        return obj.proteins.count()

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related("proteins")


@admin.register(Lineage)
class LineageAdmin(MPTTModelAdmin, CompareVersionAdmin):
    form = LineageAdminForm
    mptt_level_indent = 10
    mptt_indent_field = "protein"
    list_display = (
        "protein",
        "mutation_ellipsis",
        "rootmut_ellipsis",
        "created",
        "status",
    )
    autocomplete_fields = ("protein", "parent", "reference")
    search_fields = ("protein__name",)
    readonly_fields = (
        "created",
        "modified",
        "created_by",
        "updated_by",
        "rootmut",
        "mutation_ellipsis",
        "status",
        "errors",
        "root_node",
    )

    fieldsets = [
        (
            None,
            {
                "fields": (
                    "protein",
                    "parent",
                    "mutation",
                    "reference",
                    "root_node",
                    "rootmut",
                    "errors",
                )
            },
        ),
        (
            "Change History",
            {
                "classes": ("collapse",),
                "fields": (("created", "created_by", "modified", "updated_by")),
            },
        ),
    ]

    formfield_overrides = {MutationSetField: {"widget": TextInput(attrs={"size": "100%"})}}

    max_length = 25

    @admin.display(description="Mutation from parent")
    def mutation_ellipsis(self, obj):
        if len(str(obj.mutation)) > self.max_length:
            return f"{str(obj.mutation)[: self.max_length]}..."
        return str(obj.mutation)

    @admin.display(description="Mutation from root")
    def rootmut_ellipsis(self, obj):
        if len(str(obj.rootmut)) > self.max_length:
            return f"{str(obj.rootmut)[: self.max_length]}..."
        return str(obj.rootmut)

    @admin.display(description="OK?")
    def status(self, obj):
        if obj.parent and obj.parent.protein.seq and obj.protein.seq:
            try:
                newseq = obj.parent.protein.seq.mutate(obj.mutation)
                if newseq != obj.protein.seq:
                    return mark_safe("⚠️")
            except Exception:
                return mark_safe("❌")
        return ""

    def errors(self, obj):
        errors = validate_node(obj)
        return "\n".join(errors)

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related("protein", "parent__protein")

    def get_form(self, request, obj=None, **kwargs):
        form = super().get_form(request, **kwargs)
        qs = Lineage.objects.prefetch_related("protein").all().order_by("protein__name")
        if obj and hasattr(obj, "id"):
            qs = qs.exclude(id=obj.id)
        form.base_fields["parent"].queryset = qs
        return form
