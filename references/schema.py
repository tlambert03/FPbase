import graphene
import graphene_django_optimizer as gdo

from . import models


class Author(gdo.OptimizedDjangoObjectType):
    publications = graphene.List(lambda: Reference)

    class Meta:
        model = models.Author
        exclude = ("reference_set",)

    @gdo.resolver_hints(select_related=("publications"), only=("publications"))
    def resolve_publications(self, info):
        return self.publications.all()


class Reference(gdo.OptimizedDjangoObjectType):
    authors = graphene.List(Author)

    class Meta:
        model = models.Reference
        exclude = ("author_set",)

    @gdo.resolver_hints(select_related=("authors"), only=("authors"))
    def resolve_authors(self, info):
        return self.authors.all()


class Query(graphene.ObjectType):
    references = graphene.List(Reference)
    reference = graphene.Field(Reference, doi=graphene.String())

    def resolve_references(self, info, **kwargs):
        return gdo.query(models.Reference.objects.all(), info)

    def resolve_reference(self, info, **kwargs):
        doi = kwargs.get("doi")
        if doi is not None:
            try:
                return gdo.query(models.Reference.objects.filter(doi=doi), info).get()
            except models.Reference.DoesNotExist:
                return None
        return None
