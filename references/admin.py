from django.contrib import admin
from django.urls import reverse
from django.utils.safestring import mark_safe

from references.models import Reference, Author
from references.forms import ReferenceForm
from reversion_compare.admin import CompareVersionAdmin

from proteins.models import Protein, Excerpt


@admin.register(Author)
class AuthorAdmin(admin.ModelAdmin):
    list_display = ('full_name', 'num_refs')
    fields = (('family', 'given'), 'ref_links', 'protein_links')
    readonly_fields = ('ref_links', 'protein_links')
    search_fields = ('family',)
    ordering = ('family',)

    def num_refs(self, obj):
        return mark_safe(obj.publications.all().count())

    def ref_links(self, obj):
        refs = obj.publications.all()
        links = []
        for ref in refs:
            url = reverse("admin:references_reference_change", args=(ref.pk,))
            link = '<a href="{}">{}</a>'.format(url, ref)
            links.append(link)
        return mark_safe(", ".join(links))

    def protein_links(self, obj):
        proteins = obj.protein_contributions
        links = []
        for prot in proteins:
            url = reverse("admin:proteins_protein_change", args=(prot.pk,))
            link = '<a href="{}">{}</a>'.format(url, prot)
            links.append(link)
        return mark_safe(", ".join(links))

    ref_links.short_description = 'References'
    protein_links.short_description = 'Proteins'

    def get_queryset(self, request):
        queryset = super().get_queryset(request).prefetch_related('publications')
        return queryset


class ReferenceInline(admin.StackedInline):
    # form = StateForm
    # formset = StateFormSet
    model = Reference
    extra = 0
    can_delete = True
    show_change_link = True


class PrimaryProteinInline(admin.TabularInline):
    model = Protein.references.through
    extra = 0
    can_delete = True
    show_change_link = True


class ExcerptInline(admin.StackedInline):
    # form = StateForm
    # formset = StateFormSet
    model = Excerpt
    extra = 0
    can_delete = True
    show_change_link = True
    fields = (('created', 'created_by'), 'content', ('proteins', 'status'))
    autocomplete_fields = ('proteins', 'created_by')
    readonly_fields = ('created', )
    list_select_related = ('proteins')


@admin.register(Reference)
class ReferenceAdmin(CompareVersionAdmin):
    form = ReferenceForm
    model = Reference
    ordering = ('-date', 'citation', 'created',)
    list_display = ('id', 'citation', 'protein_links', 'title', 'date', 'doi', 'created')
    list_filter = ('created', 'modified')
    search_fields = ('pmid', 'doi', 'title', 'citation', 'firstauthor')
    inlines = (PrimaryProteinInline, ExcerptInline)
    fieldsets = [
        ('Reference', {
            'fields': ('pmid', 'doi', 'title', 'author_links', 'protein_links',
                       'secondary_proteins', 'journal', 'volume', 'pages',
                       'issue', 'date', 'refetch_info_on_save')
        }),
        ('Bleach Measurements', {
            'fields': ('bleach_links',)
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': (('created', 'created_by'), ('modified', 'updated_by'))
        })
    ]
    readonly_fields = ('author_links', 'protein_links', 'secondary_proteins', 'bleach_links',
                       'created', 'created_by', 'modified', 'updated_by')

    def author_links(self, obj):
        authors = obj.authors.all()
        links = []
        for author in authors:
            url = reverse("admin:references_author_change", args=(author.pk,))
            link = '<a href="{}">{}</a>'.format(url, author)
            links.append(link)
        return mark_safe(", ".join(links))

    def protein_links(self, obj):
        proteins = obj.primary_proteins.all()
        links = []
        for prot in proteins:
            url = reverse("admin:proteins_protein_change", args=(prot.pk,))
            link = '<a href="{}">{}</a>'.format(url, prot)
            links.append(link)
        return mark_safe(", ".join(links))

    def secondary_proteins(self, obj):
        primary = obj.primary_proteins.all()
        proteins = obj.proteins.exclude(id__in=primary)
        links = []
        for prot in proteins:
            url = reverse("admin:proteins_protein_change", args=(prot.pk,))
            link = '<a href="{}">{}</a>'.format(url, prot)
            links.append(link)
        return mark_safe(", ".join(links))

    def bleach_links(self, obj):
        links = []
        for bm in obj.bleach_measurements.all():
            url = reverse("admin:proteins_bleachmeasurement_change", args=(bm.pk,))
            link = '<a href="{}">{}</a>'.format(url, bm)
            links.append(link)
        return mark_safe(", ".join(links))
    bleach_links.short_description = 'BleachMeasurements'

    author_links.short_description = 'Authors'
    protein_links.short_description = 'Primary Proteins'
    secondary_proteins.short_description = 'Secondary Proteins'

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save(skipdoi=not form.cleaned_data['refetch_info_on_save'])

    def get_queryset(self, request):
        queryset = super().get_queryset(request).prefetch_related('authors', 'primary_proteins')
        return queryset

