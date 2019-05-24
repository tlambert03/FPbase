import graphene

import proteins.schema


class Query(proteins.schema.Query, graphene.ObjectType):
    pass


schema = graphene.Schema(query=Query)
