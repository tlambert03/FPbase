from django.contrib import admin
from django.urls import reverse
from django.utils.safestring import mark_safe

from references.models import Reference, Author
from references.forms import ReferenceForm, AuthorForm


@admin.register(Author)
class AuthorAdmin(admin.ModelAdmin):
    list_display = ('family', 'given')
    fields = (('family', 'given'), 'ref_links', 'protein_links')
    readonly_fields = ('ref_links', 'protein_links')
    search_fields = ('family',)

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


@admin.register(Reference)
class ReferenceAdmin(admin.ModelAdmin):
    form = ReferenceForm

    list_display = ('__str__',  'protein_links', 'title', 'pmid', 'doi', 'created')
    list_filter = ('created', 'modified')
    search_fields = ('pmid', 'doi', 'created_by__username', 'created_by__first_name', 'created_by__last_name', 'title')

    fieldsets = [
        ('Reference', {
            'fields': ('pmid', 'doi', 'title', 'author_links', 'protein_links', 'journal', 'volume', 'pages', 'issue', 'year')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created', 'created_by', 'modified', 'updated_by')
        })
    ]
    readonly_fields = ('title', 'author_links', 'protein_links', 'journal', 'pages', 'volume', 'issue', 'year', 'created', 'created_by', 'modified', 'updated_by')

    def author_links(self, obj):
        authors = obj.authors.all()
        links = []
        for author in authors:
            url = reverse("admin:references_author_change", args=(author.pk,))
            link = '<a href="{}">{}</a>'.format(url, author)
            links.append(link)
        return mark_safe(", ".join(links))

    def protein_links(self, obj):
        proteins = obj.protein_primary_reference.all()
        links = []
        for prot in proteins:
            url = reverse("admin:proteins_protein_change", args=(prot.pk,))
            link = '<a href="{}">{}</a>'.format(url, prot)
            links.append(link)
        return mark_safe(", ".join(links))

    author_links.short_description = 'Authors'
    protein_links.short_description = 'Proteins'

    def save_model(self, request, obj, form, change):
        if not obj.created_by:
            obj.created_by = request.user
        obj.updated_by = request.user
        obj.save()
