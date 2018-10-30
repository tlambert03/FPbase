from ..models import Lineage


def validate_lineage():
    ALL = list(Lineage.objects.all().select_related('protein'))
    for lin in ALL:
        try:
            MutationSet(lin.mutation)
        except Exception as e:
            print("Mutation {} in {}->{} failed: {}"
                .format(lin.mutation, lin.parent.protein, lin.protein, e)

    roots = Lineage.objects.filter(parent=None).select_related('protein')
    for root in roots:
        print('Validating root: ', root.protein.name)
