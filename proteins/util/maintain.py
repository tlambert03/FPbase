import re
from fpseq.mutations import MutationSet, Mutation, NEUTRAL_MUTATIONS
from django.core.exceptions import ValidationError

from ..validators import validate_mutation
from proteins.util.helpers import getprot


def check_offset(protname):
    prot = getprot(protname)
    ms = MutationSet(prot.lineage.mutation)
    offset = ms.detect_offset(prot.lineage.parent.protein.seq)
    if offset:
        print("{}:\t{} -> {}".format(protname, ms, ms.shift(offset)))


def add_missing_seqs():
    from ..models import Lineage
    for node in Lineage.objects.all():
        if not node.protein.seq and (node.parent and node.parent.protein.seq):
            seq, _ = node.parent.protein.seq.mutate(node.mutation, correct_offset=True)
            node.protein.seq = seq
            node.protein.save()
            print("saved seq for {}".format(node.protein))


def check_node_sequence_mutation_consistent(node, correct_offset=False):
    parent = node.parent.protein.seq
    if correct_offset:
        seq, _ = parent.mutate(node.mutation, correct_offset=True)
    else:
        seq = parent.mutate(node.mutation)

    if node.protein.seq and seq != node.protein.seq:
        ms = seq.mutations_to(node.protein.seq)
        #ms -= "/".join(NEUTRAL_MUTATIONS)
        return ms


def validate_node(node):
    if not node.mutation:
        return []
    errors = []
    # check that the mutation aligns with parent and actually yields the sequence
    try:
        if node.parent and node.parent.protein.seq:
            ms = check_node_sequence_mutation_consistent(node)
            if ms:
                errors.append(f'{node.parent.protein} + {node.mutation} does not match the current {node.protein} sequence (Δ: {ms})')  # noqa
    except Mutation.SequenceMismatch as e:
        errors.append(str(e).replace('parent', node.parent.protein.name))
    return errors


def check_lineages(qs=None, correct_offset=False):
    errors = dict()
    good = set()

    if not qs:
        from ..models import Lineage
        qs = Lineage.objects.all()

    for node in list(qs.prefetch_related('protein', 'parent__protein', 'root_node__protein')):
        err = validate_node(node)
        if err:
            errors[node] = err
        else:
            good.add(node)

    return errors, good
