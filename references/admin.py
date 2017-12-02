from django.contrib import admin
from django.urls import reverse
from django.utils.safestring import mark_safe

from references.models import Reference, Author
from references.forms import ReferenceForm, AuthorForm


@admin.register(Author)
class AuthorAdmin(admin.ModelAdmin):
    list_display = ('last_name', 'initials')
    fields = (('last_name', 'initials'), 'ref_links', 'protein_links')
    readonly_fields = ('ref_links', 'protein_links')
    search_fields = ('last_name',)

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

    list_display = ('__str__',  'protein_links', 'title', 'pmid', 'doi', 'created_at')
    list_filter = ('created_at', 'updated_at')
    search_fields = ('pmid', 'doi', 'added_by__username', 'added_by__first_name', 'added_by__last_name', 'title')

    fieldsets = [
        ('Reference', {
            'fields': ('pmid', 'doi', 'title', 'author_links', 'protein_links', 'journal', 'volume', 'pages', 'pubdate', 'so', 'ref')
        }),
        ('Change History', {
            'classes': ('collapse',),
            'fields': ('created_at', 'added_by', 'updated_at', 'updated_by')
        })
    ]
    readonly_fields = ('title', 'author_links', 'protein_links', 'journal', 'pages', 'volume', 'pubdate', 'so', 'ref', 'created_at', 'added_by', 'updated_at', 'updated_by')

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
        if not obj.added_by:
            obj.added_by = request.user
        obj.updated_by = request.user
        obj.save()
