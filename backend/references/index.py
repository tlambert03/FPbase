from algoliasearch_django import AlgoliaIndex
from algoliasearch_django.decorators import register

from .models import Reference


@register(Reference)
class ReferenceIndex(AlgoliaIndex):
    fields = (
        "doi",
        "journal",
        "pmid",
        "year",
        "first_author",
        "title",
        "citation",
        "date",
        "prot_primary",
        "prot_secondary",
        "_excerpts",
        "url",
    )
