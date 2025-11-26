import graphene

from proteins import models
from proteins.schema.types import Protein


class ProteinNode(Protein):
    class Meta:
        model = models.Protein
        interfaces = (graphene.relay.Node,)
        fields = "__all__"
