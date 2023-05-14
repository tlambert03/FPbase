from algoliasearch_django import AlgoliaIndex
from algoliasearch_django.decorators import register

from .models import Organism, Protein


@register(Protein)
class ProteinIndex(AlgoliaIndex):
    fields = (
        "name",
        "uuid",
        "aliases",
        "pdb",
        "genbank",
        "uniprot",
        "ipg_id",
        "_agg",
        "img_url",
        "switchType",
        "url",
        "date_published",
        "created",
        "rank",
        "ga_views",
        "n_faves",
        "n_cols",
        "ex",
        "em",
        "pka",
        "ec",
        "qy",
        "em_css",
        "local_brightness",
        "seq",
        "first_author",
        "cofactor",
        "color",
    )
    should_index = "is_visible"
    tags = "tags"


@register(Organism)
class OrganismIndex(AlgoliaIndex):
    fields = ("scientific_name", "division", "url")
