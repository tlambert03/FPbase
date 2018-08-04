from django.contrib import admin
from django.urls import reverse
from django.utils.safestring import mark_safe
from django.db.models import Count
from proteins.models import (Protein, State, StateTransition, Organism,
                             BleachMeasurement, Spectrum, Dye,
                             Light, Filter, Camera, Mutation, Microscope,
                             OpticalConfig, FilterPlacement, Fluorophore,
                             ProteinCollection)
from reversion_compare.admin import CompareVersionAdmin
from reversion.admin import VersionAdmin
# from reversion.models import Version


# ############ INLINES ###############

# placeholder... not yet used
# @admin.register(Mutation)
# class MutationAdmin(VersionAdmin):
#     model = Mutation


class SpectrumOwner(object):
    list_display = ('__str__', 'spectra', 'created_by', 'created')
    list_select_related = ('created_by',)
    list_filter = ('created', 'manufacturer')
    search_fields = ('name',)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.readonly_fields = self.readonly_fields + ('spectra',)

    def spectra(self, obj):
        def _makelink(sp):
            url = reverse("admin:proteins_spectrum_change", args=(sp.pk,))
            return '<a href="{}">{}</a>'.format(url, sp.get_subtype_display())
        links = []
        if isinstance(obj, Fluorophore):
            [links.append(_makelink(sp)) for sp in obj.spectra.all()]
        else:
            links.append(_makelink(obj.spectrum))
        return mark_safe(", ".join(links))
    spectra.short_description = 'spectra'

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related('spectrum')


class MultipleSpectraOwner(SpectrumOwner):
    def get_queryset(self, request):
        qs = super(SpectrumOwner, self).get_queryset(request)
        return qs.prefetch_related('spectra')


class BleachInline(admin.TabularInline):
    model = BleachMeasurement
    autocomplete_fields = ("reference",)
    extra = 1


class StateInline(MultipleSpectraOwner, admin.StackedInline):
    # form = StateForm
    # formset = StateFormSet
    model = State
    extra = 0
    can_delete = True
    show_change_link = True
    fieldsets = [
        (None, {
            'fields': (('name', 'slug', 'is_dark',),)
        }),
        (None, {
            'fields': (('ex_max', 'em_max'), ('ext_coeff', 'qy'), ('pka', 'maturation'), 'lifetime', 'bleach_links', 'spectra')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': (('created', 'created_by'), ('modified', 'updated_by'))
        })
    ]
    readonly_fields = ('slug', 'bleach_links', 'created', 'created_by', 'modified', 'updated_by')

    def bleach_links(self, obj):
        links = []
        for bm in obj.bleach_measurements.all():
            url = reverse("admin:proteins_bleachmeasurement_change", args=(bm.pk,))
            link = '<a href="{}">{}</a>'.format(url, bm)
            links.append(link)
        return mark_safe(", ".join(links))
    bleach_links.short_description = 'BleachMeasurements'


# class FRETpairInline(admin.TabularInline):
#     model = FRETpair
#     autocomplete_fields = ("pair_references", 'acceptor')
#     extra = 1
#     can_delete = True
#     show_change_link = True
#     fk_name = 'donor'  # or 'acceptor'
#     fields = ('acceptor', 'radius', 'pair_references')


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


@admin.register(Dye)
class DyeAdmin(MultipleSpectraOwner, VersionAdmin):
    model = Dye


@admin.register(Filter)
class FilterAdmin(SpectrumOwner, VersionAdmin):
    model = Filter


@admin.register(Camera)
class CameraAdmin(SpectrumOwner, VersionAdmin):
    model = Camera


@admin.register(Mutation)
class MutationAdmin(VersionAdmin):
    model = Mutation


@admin.register(Spectrum)
class SpectrumAdmin(admin.ModelAdmin):
    model = Spectrum
    autocomplete_fields = []
    list_select_related = ('owner_state__protein', 'owner_filter', 'owner_camera', 'owner_light', 'owner_dye', 'created_by')
    list_display = ('__str__', 'category', 'subtype', 'owner', 'created_by')
    list_filter = ('created', 'category', 'subtype')
    readonly_fields = ('owner', 'name', 'created')
    search_fields = ('owner_state__protein__name', 'owner_filter__name', 'owner_camera__name', 'owner_light__name', 'owner_dye__name')

    def get_fields(self, request, obj=None):
        fields = []
        if not obj or not obj.category:
            own = ['owner_state', 'owner_filter', 'owner_camera', 'owner_light', 'owner_dye']
        elif obj.category == Spectrum.PROTEIN:
            own = ['owner_state']
        else:
            own = ['owner_' + obj.get_category_display().split(' ')[0].lower()]
        fields.extend(own)
        self.autocomplete_fields.extend(own)
        fields += ['category', 'subtype', 'data', 'ph', 'solvent',
                   ('created', 'created_by', 'updated_by')]
        return fields

    def owner(self, obj):
        url = reverse("admin:proteins_{}_change".format(obj.owner._meta.model.__name__.lower()), args=(obj.owner.pk,))
        link = '<a href="{}">{}</a>'.format(url, obj.owner)
        return mark_safe(link)
    owner.short_description = 'Owner'

    # def get_queryset(self, request):
    #     qs = super().get_queryset(request)
    #     return qs.prefetch_related('owner_state__protein')


@admin.register(BleachMeasurement)
class BleachMeasurementAdmin(VersionAdmin):
    model = BleachMeasurement
    autocomplete_fields = ('state', 'reference')
    list_select_related = ('state', 'state__protein',)


@admin.register(State)
class StateAdmin(CompareVersionAdmin):
    # form = StateForm
    model = State
    list_select_related = ('protein',)
    search_fields = ('protein__name',)
    list_display = ('__str__', 'protein_link', 'ex_max', 'em_max', 'created_by')
    list_filter = ('created', 'modified', 'created_by__username')
    inlines = (BleachInline,)
    fieldsets = [
        (None, {
            'fields': (('name', 'slug', 'is_dark',),)
        }),
        (None, {
            'fields': (('ex_max', 'em_max'), ('ext_coeff', 'qy'), ('pka', 'maturation'), 'lifetime')
        }),
    ]

    def protein_link(self, obj):
        url = reverse("admin:proteins_protein_change", args=([obj.protein.pk]))
        return mark_safe('<a href="{}">{}</a>'.format(url, obj.protein))
    protein_link.short_description = 'Protein'


class StateTransitionAdmin(VersionAdmin):
    model = StateTransition
    list_select_related = ('protein', 'from_state', 'to_state')
    autocomplete_fields = ('protein', 'from_state', 'to_state')


class StateTransitionInline(admin.TabularInline):
    model = StateTransition
    extra = 1

    def formfield_for_foreignkey(self, db_field, request=None, **kwargs):
        field = super().formfield_for_foreignkey(db_field, request, **kwargs)

        # restrict
        if db_field.related_model == State:
            obj = getattr(request, 'object', None)
            if request is not None:
                field.queryset = field.queryset.filter(protein=obj)
            else:
                field.queryset = field.queryset.none()

        return field


# class FRETpairAdmin(CompareVersionAdmin):
#     list_display = ('__str__', 'donor', 'acceptor', 'radius', 'created_by', 'created', 'modified')
#     list_filter = ('created', 'modified')
#     search_fields = ('name', 'donor', 'acceptor', 'created_by__username', 'created_by__first_name', 'created_by__last_name')
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
    list_display = ('scientific_name', 'id', 'created_by', 'created', 'modified')
    list_filter = ('created', 'modified')
    search_fields = ('scientific_name', 'common_name', 'id', 'created_by__username', 'created_by__first_name', 'created_by__last_name')

    fieldsets = [
        ('Organism', {
            'fields': ('id', 'scientific_name', 'common_name', 'rank', 'division', 'genus', 'species')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('scientific_name', 'common_name', 'division', 'rank', 'genus', 'species', 'created', 'created_by', 'modified', 'updated_by')

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()


@admin.register(Protein)
class ProteinAdmin(CompareVersionAdmin):
    autocomplete_fields = ('parent_organism', 'references', 'primary_reference')
    list_display = ('__str__', 'ipg_id', 'switch_type', 'created', 'modified', 'states_all_count')
    list_filter = ('created', 'modified', 'switch_type', 'created_by__username', 'status')
    search_fields = ('name', 'aliases', 'slug', 'ipg_id', 'created_by__username', 'created_by__first_name', 'created_by__last_name')
    prepopulated_fields = {'slug': ('name',)}
    inlines = [
        StateInline, StateTransitionInline
    ]
    fieldsets = [
        (None, {
            'fields': (('name', 'slug',), ('aliases', 'chromophore'),
                       ('seq', 'seq_validated'), ('ipg_id', 'genbank', 'uniprot', 'pdb'),
                       ('parent_organism', 'switch_type'), ('agg', 'mw'))
        }),
        (None, {
            'fields': (('primary_reference', 'references'),)
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': (('created', 'created_by'), ('modified', 'updated_by'))
        })
    ]
    readonly_fields = ('created', 'created_by', 'modified', 'updated_by', 'switch_type')

    def states_all_count(self, obj):
        if obj:
            return obj.states.all().count()
        else:
            return 0

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related('states')

    def get_form(self, request, obj=None, **kwargs):
        request.object = obj
        return super().get_form(request, obj, **kwargs)

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.status = 'approved'
        obj.save()

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
    inlines = (FilterPlacementInline, )
    list_display = ('__str__', 'microscope', 'owner_link',)
    fields = (('name', 'microscope'), ('laser', 'light'), ('camera', 'owner'))

    def owner_link(self, obj):
        if obj.microscope and obj.microscope.owner:
            url = reverse("admin:users_user_change", args=([obj.microscope.owner.pk]))
            return mark_safe('<a href="{}">{}</a>'.format(url, obj.microscope.owner))
    owner_link.short_description = 'Owner'

    def save_model(self, request, obj, form, change):
        obj.save()
        obj.microscope.save()

    def get_queryset(self, request):
        qs = super().get_queryset(request)
        return qs.prefetch_related('microscope__owner')


# class OpticalConfigInline(admin.TabularInline):
#     model = OpticalConfig
#     fields = ('name', 'light', 'laser', 'camera')
#     extra = 1  # how many rows to show


@admin.register(Microscope)
class MicroscopeAdmin(admin.ModelAdmin):
    model = Microscope
    autocomplete_fields = ('extra_lights', 'extra_cameras')
    readonly_fields = ('configs', )
    list_display = ('__str__', 'owner_link', 'created', 'numocs',)
    list_select_related = ('owner', )
    list_filter = ('created', )

    def owner_link(self, obj):
        if obj.owner:
            url = reverse("admin:users_user_change", args=([obj.owner.pk]))
            return mark_safe('<a href="{}">{}</a>'.format(url, obj.owner))
    owner_link.short_description = 'Owner'

    def numocs(self, obj):
        return obj.optical_configs.count()
    numocs.admin_order_field = 'oc_count'

    def configs(self, obj):
        def _makelink(oc):
            url = reverse("admin:proteins_opticalconfig_change", args=(oc.pk,))
            return '<a href="{}">{}</a>'.format(url, oc)
        links = []
        [links.append(_makelink(oc)) for oc in obj.optical_configs.all()]
        return mark_safe(", ".join(links))
    configs.short_description = 'Optical Configs'

    def get_queryset(self, request):
        qs = super().get_queryset(request).annotate(oc_count=Count('optical_configs'))
        return qs.prefetch_related('optical_configs')


def make_private(modeladmin, request, queryset):
    # note, this will fail if the list is ordered by numproteins
    queryset.update(private=True)


make_private.short_description = "Mark selected collections as private"


@admin.register(ProteinCollection)
class ProteinCollectionAdmin(admin.ModelAdmin):
    model = ProteinCollection
    list_display = ('__str__', 'owner_link', 'created', 'private', 'numproteins',)
    list_filter = ('created', 'private')
    readonly_fields = ('numproteins', 'owner_link')
    list_select_related = ('owner', )
    autocomplete_fields = ('proteins', )
    search_fields = ('name', 'proteins__name')
    actions = [make_private]

    def owner_link(self, obj):
        url = reverse("admin:users_user_change", args=([obj.owner.pk]))
        return mark_safe('<a href="{}">{}</a>'.format(url, obj.owner))
    owner_link.short_description = 'Owner'

    def numproteins(self, obj):
        return obj.proteins.count()
    numproteins.admin_order_field = 'proteins_count'

    def get_queryset(self, request):
        qs = super().get_queryset(request).annotate(proteins_count=Count('proteins'))
        return qs.prefetch_related('proteins')
