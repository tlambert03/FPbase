from django.contrib import admin
from proteins.models import Protein, State, Organism, FRETpair
from proteins.forms import ProteinForm, StateForm, StateFormSet


class StateInline(admin.StackedInline):
    form = StateForm
    formset = StateFormSet
    model = State
    extra = 0
    can_delete = True
    fieldsets = [
        ('State', {
            'fields': ('default', 'state_name', 'state_id')
        }),
        ('Attributes', {
            'classes': ('collapse',),
            'fields': ('ex_max', 'em_max', 'ext_coeff', 'qy', 'pka', 'maturation', 'lifetime',)
        }),
        # ('Switching', {
        #     'fields': ('trans_wave', 'to_state')
        # }),
        # ('Bleach Measurements', {
        #     'fields': ()
        # }),
        ('Spectra', {
            'classes': ('collapse',),
            'fields': ('ex_spectra', 'em_spectra',)
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created_at', 'added_by', 'updated_at', 'updated_by')
        })
    ]
    readonly_fields = ('state_id', 'created_at', 'added_by', 'updated_at', 'updated_by')


@admin.register(Protein)
class ProteinAdmin(admin.ModelAdmin):
    form = ProteinForm
    list_display = ('__str__', 'gb_prot', 'gb_nuc', 'switch_type', 'created_at', 'updated_at')
    list_filter = ('created_at', 'updated_at', 'switch_type')
    search_fields = ('name', 'slug', 'gb_prot', 'gb_nuc', 'added_by__username', 'added_by__first_name', 'added_by__last_name')
    prepopulated_fields = {'slug': ('name',)}
    inlines = [
        StateInline,
    ]
    fieldsets = [
        ('Protein', {
            'fields': ('name', 'slug', 'ipg_id', 'seq', 'gb_prot', 'gb_nuc', 'parent_organism', 'switch_type', 'agg', 'mw')
        }),
        ('References', {
            'classes': ('collapse',),
            'fields': ('primary_reference', 'references')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created_at', 'added_by', 'updated_at', 'updated_by')
        })
    ]
    readonly_fields = ('created_at', 'added_by', 'updated_at', 'updated_by')
    filter_horizontal = ('references',)

    def save_model(self, request, obj, form, change):
        if not obj.added_by:
            obj.added_by = request.user
        obj.updated_by = request.user
        obj.save()

    def save_formset(self, request, form, formset, change):
        instances = formset.save(commit=False)
        for instance in instances:
            if not instance.added_by:
                instance.added_by = request.user
            instance.updated_by = request.user
            instance.save()
        formset.save()


@admin.register(FRETpair)
class FRETpairAdmin(admin.ModelAdmin):
    list_display = ('__str__', 'donor', 'acceptor', 'radius', 'added_by', 'created_at', 'updated_at')
    list_filter = ('created_at', 'updated_at')
    search_fields = ('name', 'donor', 'acceptor', 'added_by__username', 'added_by__first_name', 'added_by__last_name')
    fieldsets = [
        ('FRET Pair', {
            'fields': ('name', 'donor', 'acceptor', 'radius')
        }),
        ('References', {
            'classes': ('collapse',),
            'fields': ('pair_references',)
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created_at', 'added_by', 'updated_at', 'updated_by')
        })
    ]
    readonly_fields = ('name', 'created_at', 'added_by', 'updated_at', 'updated_by')

    def save_model(self, request, obj, form, change):
        if not obj.added_by:
            obj.added_by = request.user
        obj.updated_by = request.user
        obj.save()


@admin.register(Organism)
class OrganismAdmin(admin.ModelAdmin):
    list_display = ('scientific_name', 'tax_id', 'added_by', 'created_at', 'updated_at')
    list_filter = ('created_at', 'updated_at')
    search_fields = ('scientific_name', 'common_name', 'tax_id', 'added_by__username', 'added_by__first_name', 'added_by__last_name')

    fieldsets = [
        ('Organism', {
            'fields': ('tax_id', 'scientific_name', 'common_name', 'division', 'genus', 'species')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created_at', 'added_by', 'updated_at', 'updated_by')
        })
    ]
    readonly_fields = ('scientific_name', 'common_name', 'division', 'genus', 'species', 'created_at', 'added_by', 'updated_at', 'updated_by')

    def save_model(self, request, obj, form, change):
        if not obj.added_by:
            obj.added_by = request.user
        obj.updated_by = request.user
        obj.save()
