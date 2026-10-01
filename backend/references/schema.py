import graphene
import graphene_django_optimizer as gdo

from proteins.visibility import hidden_ids
from references import models


class Author(gdo.OptimizedDjangoObjectType):
    publications = graphene.List(lambda: Reference)

    class Meta:
        model = models.Author
        exclude = ("reference_set",)

    @gdo.resolver_hints(model_field="publications")
    def resolve_publications(self, info):
        return self.publications.all()


class Reference(gdo.OptimizedDjangoObjectType):
    authors = graphene.List(Author)

    class Meta:
        model = models.Reference
        exclude = ("author_set",)

    @gdo.resolver_hints(model_field="authors")
    def resolve_authors(self, info):
        return self.authors.all()

    # (not what belongs to a hidden protein: see `hidden_ids`)
    @gdo.resolver_hints(model_field="spectra")
    def resolve_spectra(self, info):
        hidden = hidden_ids(info.context).spectra
        return [spectrum for spectrum in self.spectra.all() if spectrum.id not in hidden]

    @gdo.resolver_hints(model_field="oser_measurements")
    def resolve_oser_measurements(self, info):
        hidden = hidden_ids(info.context).oser_measurements
        return [m for m in self.oser_measurements.all() if m.id not in hidden]


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
