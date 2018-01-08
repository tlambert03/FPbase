from django.contrib import admin
from proteins.models import Protein, State, StateTransition, Organism, FRETpair, BleachMeasurement
from proteins.forms import ProteinForm, StateForm  # StateFormSet
from reversion_compare.admin import CompareVersionAdmin


class BleachInline(admin.StackedInline):
    model = BleachMeasurement

@admin.register(State)
class StateAdmin(CompareVersionAdmin):
    form = StateForm
    model = State
    search_fields = ('protein__name',)
    inlines = (BleachInline,)

class StateInline(admin.StackedInline):
    form = StateForm
#    formset = StateFormSet
    model = State
    extra = 0
    can_delete = True
    fieldsets = [
        ('State', {
            'fields': ('is_dark', 'name', 'slug')
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
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('slug', 'created', 'created_by', 'modified', 'updated_by')


class StateTransitionInline(admin.StackedInline):
    model = StateTransition


@admin.register(Protein)
class ProteinAdmin(CompareVersionAdmin):
    form = ProteinForm
    list_display = ('__str__', 'ipg_id', 'switch_type', 'created', 'modified')
    list_filter = ('created', 'modified', 'switch_type', 'created_by__username')
    search_fields = ('name', 'slug', 'ipg_id', 'created_by__username', 'created_by__first_name', 'created_by__last_name')
    prepopulated_fields = {'slug': ('name',)}
    inlines = [
        StateInline, StateTransitionInline,
    ]
    fieldsets = [
        ('Protein', {
            'fields': ('name', 'newtest', 'slug', 'ipg_id', 'seq', 'gb_prot', 'gb_nuc', 'parent_organism', 'switch_type', 'agg', 'mw')
        }),
        ('References', {
            'classes': ('collapse',),
            'fields': ('primary_reference', 'references')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('created', 'created_by', 'modified', 'updated_by', 'switch_type')
    filter_horizontal = ('references',)

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()

    def save_formset(self, request, form, formset, change):
        instances = formset.save(commit=False)
        for instance in instances:
            if not instance.created_by:
                instance.created_by = request.user
            instance.updated_by = request.user
            instance.save()
        formset.save()


@admin.register(FRETpair)
class FRETpairAdmin(CompareVersionAdmin):
    list_display = ('__str__', 'donor', 'acceptor', 'radius', 'created_by', 'created', 'modified')
    list_filter = ('created', 'modified')
    search_fields = ('name', 'donor', 'acceptor', 'created_by__username', 'created_by__first_name', 'created_by__last_name')
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
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('name', 'created', 'created_by', 'modified', 'updated_by')

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()


@admin.register(Organism)
class OrganismAdmin(CompareVersionAdmin):
    list_display = ('scientific_name', 'tax_id', 'created_by', 'created', 'modified')
    list_filter = ('created', 'modified')
    search_fields = ('scientific_name', 'common_name', 'tax_id', 'created_by__username', 'created_by__first_name', 'created_by__last_name')

    fieldsets = [
        ('Organism', {
            'fields': ('tax_id', 'scientific_name', 'common_name', 'division', 'genus', 'species')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('scientific_name', 'common_name', 'division', 'genus', 'species', 'created', 'created_by', 'modified', 'updated_by')

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()
