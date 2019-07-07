import graphene

import proteins.schema
import references.schema


class Query(proteins.schema.Query, references.schema.Query, graphene.ObjectType):
    pass


schema = graphene.Schema(query=Query)
