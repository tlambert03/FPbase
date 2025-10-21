import graphene
from query_optimizer import DjangoListField, DjangoObjectType, optimize

from . import models


class Author(DjangoObjectType):
    publications = DjangoListField(lambda: Reference)

    class Meta:
        model = models.Author
        exclude = ("reference_set",)


class Reference(DjangoObjectType):
    authors = DjangoListField(Author)

    class Meta:
        model = models.Reference
        exclude = ("author_set",)


class Query(graphene.ObjectType):
    references = graphene.List(Reference)
    reference = graphene.Field(Reference, doi=graphene.String())

    def resolve_references(self, info, **kwargs):
        return optimize(models.Reference.objects.all(), info)

    def resolve_reference(self, info, **kwargs):
        doi = kwargs.get("doi")
        if doi is not None:
            try:
                return optimize(models.Reference.objects.filter(doi=doi), info).get()
            except models.Reference.DoesNotExist:
                return None
        return None
