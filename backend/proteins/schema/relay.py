import graphene

from proteins import models
from proteins.schema.types import Protein


class ProteinNode(Protein):
    class Meta:
        model = models.Protein
        interfaces = (graphene.relay.Node,)
        fields = "__all__"

    @classmethod
    def get_queryset(cls, queryset, info):
        # (every connection of proteins comes through here: allProteins, a reference's, ...)
        return super().get_queryset(queryset.exclude(status=models.Protein.STATUS.hidden), info)
