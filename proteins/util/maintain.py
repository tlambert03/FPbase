import re
from fpseq.mutations import MutationSet, Mutation, NEUTRAL_MUTATIONS
from collections import defaultdict
from django.core.exceptions import ValidationError

from ..models import Lineage
from ..validators import validate_mutation
from proteins.util.helpers import getprot


def check_offset(protname):
    prot = getprot(protname)
    ms = MutationSet(prot.lineage.mutation)
    offset = ms.detect_offset(prot.lineage.parent.protein.seq)
    if offset:
        print("{}:\t{} -> {}".format(protname, ms, ms.shift(offset)))


def add_missing_seqs():
    for node in Lineage.objects.all():
        if not node.protein.seq and (node.parent and node.parent.protein.seq):
            seq, _ = node.parent.protein.seq.mutate(node.mutation, correct_offset=True)
            node.protein.seq = seq
            node.protein.save()
            print("saved seq for {}".format(node.protein))


def check_lineages(correct_offset=False):
    errors = defaultdict(set)  # each is tuple of (Node, error)
    good = set()

    def recursive_check_children(node):
        for mut in re.split(' |,|/|\\\\', str(node.mutation)):
            if mut:
                try:
                    validate_mutation(mut)
                except ValidationError as e:
                    errors[node.protein.name].add('Bad Mutation String: {}'.format(e))

        try:
            if node.parent and node.parent.protein.seq:
                root = node.get_root()
                parent = node.parent.protein.seq
                if correct_offset:
                    seq, _ = parent.mutate(node.mutation, correct_offset=True)
                else:
                    seq = parent.mutate(node.mutation)

                if node.protein.seq:
                    if seq != node.protein.seq:
                        ms = seq.mutations_to(node.protein.seq, root.protein.seq)
                        ms -= "/".join(NEUTRAL_MUTATIONS)
                        if ms:
                            errors[node.protein.name].add(
                                'Sequence does not match Parent sequence + mutation.  diff: {}'
                                .format(ms))
                    else:
                        good.add(node.protein.name)
                else:
                    seq, _ = parent.mutate(node.mutation, correct_offset=True)
                    # node.protein.seq = seq
                    # node.protein.save()
        except Mutation.SequenceMismatch as e:
            errors[node.protein.name].add('SequenceMismatch: {}'.format(e))
            return
        except ValueError as e:
            errors[node.protein.name].add('ValueError: {}'.format(e))
            return
        for c in node.get_children():
            recursive_check_children(c)

    roots = Lineage.objects.filter(parent=None).select_related('protein')
    [recursive_check_children(root) for root in roots]
    return errors, good


def validate_lineage():
    ALL = list(Lineage.objects.all().select_related('protein'))
    for lin in ALL:
        try:
            MutationSet(lin.mutation)
        except Exception as e:
            print("Mutation {} in {}->{} failed: {}"
                .format(lin.mutation, lin.parent.protein, lin.protein, e))

    errors = defaultdict(set)  # each is tuple of (Node, error)
    good = set()
    roots = Lineage.objects.filter(parent=None).select_related('protein')
    leaves = Lineage.objects.filter(children=None).exclude(parent=None).select_related('protein')
    for leaf in leaves:
        ancestors = list(leaf.get_ancestors()) + [leaf]
        root = ancestors.pop(0)
        seq = root.protein.seq
        for ancestor in ancestors:
            for mut in re.split(' |,|/|\\\\', ancestor.mutation):
                if mut:
                    try:
                        validate_mutation(mut)
                    except ValidationError as e:
                        errors[ancestor.protein.name].add('Bad Mutation String: {}'.format(e))

            try:
                seq = seq.mutate(ancestor.mutation)
                # if not ancestor.protein.seq:
                #    ancestor.protein.seq = seq
                #    ancestor.protein.save()
                if ancestor.protein.seq:
                    if seq != ancestor.protein.seq:
                        ms = seq.mutations_to(ancestor.protein.seq, root.protein.seq)
                        ms -= 'M1_S2insV/H231L/V2del/Q80R/L232H'
                        if ms:
                            errors[ancestor.protein.name].add(
                                'Sequence does not match Parent sequence + mutation.  diff: {}'
                                .format(ms))
                    else:
                        good.add(ancestor.protein.name)
                else:
                    ancestor.protein.seq = seq
                    ancestor.protein.save()
            except Mutation.SequenceMismatch as e:
                errors[ancestor.protein.name].add('SequenceMismatch: {}'.format(e))
                break
            except ValueError as e:
                errors[ancestor.protein.name].add('ValueError: {}'.format(e))
                break
    print(errors)

    for root in roots:
        print('Validating root: ', root.protein.name)
