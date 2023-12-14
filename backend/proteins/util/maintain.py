from typing import TYPE_CHECKING, cast

from fpseq.mutations import Mutation, MutationSet
from proteins.util.helpers import getprot

if TYPE_CHECKING:
    from fpseq import FPSeq
    from proteins.models import Lineage


def check_offset(protname):
    prot = getprot(protname)
    ms = MutationSet(prot.lineage.mutation)
    offset = ms.detect_offset(prot.lineage.parent.protein.seq)
    if offset:
        print(f"{protname}:\t{ms} -> {ms.shift(offset)}")


def add_missing_seqs():
    from ..models import Lineage

    for node in Lineage.objects.all():
        if not node.protein.seq and (node.parent and node.parent.protein.seq):
            seq, _ = node.parent.protein.seq.mutate(node.mutation, correct_offset=True)
            node.protein.seq = seq
            node.protein.save()
            print(f"saved seq for {node.protein}")


def check_node_sequence_mutation_consistent(node, correct_offset=False):
    parent = cast("FPSeq", node.parent.protein.seq)
    if correct_offset:
        seq, _ = parent.mutate(node.mutation, correct_offset=True)
    else:
        seq = parent.mutate(node.mutation)

    if node.protein.seq and seq != node.protein.seq:
        ms = seq.mutations_to(node.protein.seq)
        # ms -= "/".join(NEUTRAL_MUTATIONS)
        return ms


def suggested_switch_type(protein):
    """return the "apparent" switch type based on states and transitions
    for best performance, pre-annotate the protein with ndark and nfrom:

        .annotate(ndark=Count('states', filter=Q(states__is_dark=True)))
        .annotate(nfrom=Count('transitions__from_state', distinct=True))
    """
    nstates = protein.states.count()
    if not nstates:
        return None
    if nstates == 1:
        return protein.BASIC
    # 2 or more states...
    n_transitions = protein.transitions.count()
    if hasattr(protein, "ndark"):
        darkstates = protein.ndark
    else:
        darkstates = protein.states.filter(is_dark=True).count()
    if not n_transitions:
        return protein.OTHER
    elif nstates == 2:
        # 2 transitions with unique from_states
        if hasattr(protein, "nfrom"):
            nfrom = protein.nfrom
        else:
            nfrom = len(set(protein.transitions.values_list("from_state", flat=True)))
        if nfrom >= 2:
            return protein.PHOTOSWITCHABLE
        if darkstates == 0:
            return protein.PHOTOCONVERTIBLE
        if darkstates == 1:
            return protein.PHOTOACTIVATABLE
        if darkstates > 1:
            return None
    elif nstates > 2:
        return protein.MULTIPHOTOCHROMIC


def validate_switch_type(protein):
    """returns False if the protein has an unusual switch type
    for its states & transitions.
    """
    return protein.switch_type == suggested_switch_type(protein)


def validate_node(node: "Lineage") -> list[str]:
    if not node.mutation:
        return []
    errors = []
    # check that the mutation aligns with parent and actually yields the sequence
    try:
        if node.parent and node.parent.protein.seq:
            ms = check_node_sequence_mutation_consistent(node)
            if ms:
                errors.append(
                    f"{node.parent.protein} + {node.mutation} does not match "
                    + f"the current {node.protein} sequence (Δ: {ms})"
                )
    except Mutation.SequenceMismatch as e:
        errors.append(str(e).replace("parent", node.parent.protein.name))
    # except Exception as e:
    # errors.append(str(e))
    return errors


def check_lineages(qs=None, correct_offset=False):
    errors = {}
    good = set()

    if not qs:
        from ..models import Lineage

        qs = Lineage.objects.all()

    for node in list(qs.prefetch_related("protein", "parent__protein", "root_node__protein")):
        err = validate_node(node)
        if err:
            errors[node] = err
        else:
            good.add(node)

    return errors, good
