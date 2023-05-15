import graphene

from .. import models
from .types import Protein


class ProteinNode(Protein):
    class Meta:
        model = models.Protein
        interfaces = (graphene.relay.Node,)
        fields = "__all__"
