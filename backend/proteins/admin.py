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
    Excerpt,
    Filter,
    FilterPlacement,
    Fluorophore,
    Light,
    Lineage,
    Microscope,
    OpticalConfig,
    Organism,
    OSERMeasurement,
    Protein,
    ProteinCollection,
    Spectrum,
    State,
    StateTransition,
)
from proteins.models.lineage import MutationSetField
from proteins.util.maintain import validate_node

# from reversion.models import Version

# ############ INLINES ###############

# placeholder... not yet used
# @admin.register(Mutation)
# class MutationAdmin(VersionAdmin):
#     model = Mutation


class SpectrumOwner:
    list_display = ("__str__", "spectra", "created_by", "created")
    list_select_related = ("created_by",)
    list_filter = ("created", "manufacturer")
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
        if isinstance(obj, Fluorophore):
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


class StateInline(MultipleSpectraOwner, admin.StackedInline):
    # form = StateForm
    # formset = StateFormSet
    model = State
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
                    ("twop_ex_max", "twop_peakGM", "twop_qy"),
                    ("pka", "maturation"),
                    "lifetime",
                    "bleach_links",
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
        "bleach_links",
        "created",
        "created_by",
        "modified",
        "updated_by",
    )

    @admin.display(description="BleachMeasurements")
    def bleach_links(self, obj):
        links = []
        for bm in obj.bleach_measurements.all():
            url = reverse("admin:proteins_bleachmeasurement_change", args=(bm.pk,))
            link = f'<a href="{url}">{bm}</a>'
            links.append(link)
        return mark_safe(", ".join(links))


class LineageInline(admin.TabularInline):
    model = Lineage
    autocomplete_fields = ("parent", "root_node")
    readonly_fields = ("root_node", "created_by", "created", "updated_by")
    fields = ("parent", "mutation", "created_by", "created")


# ############# MODELS ###############


# @admin.register(Version)
# class myVersionAdmin(admin.ModelAdmin):
#     model = Version
#     # autocomplete_fields = ('revision',)
#
#     def get_queryset(self, request):
#         qs = super().get_queryset(request)
#         return qs.prefetch_related('revision')


@admin.register(Light)
class LightAdmin(SpectrumOwner, admin.ModelAdmin):
    model = Light
    ordering = ("-created",)


@admin.register(Dye)
class DyeAdmin(MultipleSpectraOwner, VersionAdmin):
    model = Dye
    ordering = ("-created",)


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
        "owner_state__protein",
        "owner_filter",
        "owner_camera",
        "owner_light",
        "owner_dye",
        "created_by",
    )
    list_display = ("__str__", "category", "subtype", "owner", "created_by")
    list_filter = ("status", "created", "category", "subtype")
    readonly_fields = ("owner", "name", "created", "modified")
    search_fields = (
        "owner_state__protein__name",
        "owner_filter__name",
        "owner_camera__name",
        "owner_light__name",
        "owner_dye__name",
    )

    def get_fields(self, request, obj=None):
        fields = []
        if not obj or not obj.category:
            own = [
                "owner_state",
                "owner_filter",
                "owner_camera",
                "owner_light",
                "owner_dye",
            ]
        elif obj.category == Spectrum.PROTEIN:
            own = ["owner_state"]
        else:
            own = ["owner_" + obj.get_category_display().split(" ")[0].lower()]
        fields.extend(own)
        self.autocomplete_fields.extend(own)
        fields += [
            "category",
            "subtype",
            "data",
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
        url = reverse(
            f"admin:proteins_{obj.owner._meta.model.__name__.lower()}_change",
            args=(obj.owner.pk,),
        )
        link = f'<a href="{url}">{obj.owner}</a>'
        return mark_safe(link)

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
                    ("twop_ex_max", "twop_peakGM", "twop_qy"),
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


# class FRETpairAdmin(CompareVersionAdmin):
#     list_display = ('__str__', 'donor', 'acceptor', 'radius', 'created_by', 'created', 'modified')
#     list_filter = ('created', 'modified')
#     search_fields = ('name', 'donor', 'acceptor', 'created_by__username',
#                      'created_by__first_name', 'created_by__last_name')
#     fieldsets = [
#         ('FRET Pair', {
#             'fields': ('name', 'donor', 'acceptor', 'radius')
#         }),
#         ('References', {
#             'classes': ('collapse',),
#             'fields': ('pair_references',)
#         }),
#         ('Change History', {
#             'classes': ('collapse',),
#             'fields': ('created', 'created_by', 'modified', 'updated_by')
#         })
#     ]
#     readonly_fields = ('name', 'created', 'created_by', 'modified', 'updated_by')
#
#     def save_model(self, request, obj, form, change):
#         if not obj.created_by:
#             obj.created_by = request.user
#         obj.updated_by = request.user
#         obj.save()


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


# class OpticalConfigInline(admin.TabularInline):
#     model = OpticalConfig
#     fields = ('name', 'light', 'laser', 'camera')
#     extra = 1  # how many rows to show


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
