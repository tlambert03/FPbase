try:
    from parasail import blosum62, nw_banded, nw_trace_scan_sat
except ImportError:
    import warnings

    warnings.warn("ERROR!!! could not import parasail... will not be able to align", stacklevel=2)
from .util import chunked_lines


def align_seqs(query, target, gop=5, gep=1, band_size=0):
    """basic parasail global alignment of two sequences
    result is wrapped in ParasailAlignment Class"""
    query = str(query)
    target = str(target)
    if band_size:
        result = nw_banded(target, query, gop, gep, band_size, blosum62)
    else:
        result = nw_trace_scan_sat(target, query, gop, gep, blosum62)
    return ParasailAlignment(result)


def parental_numbering(aseq1, aseq2):
    """given two ALIGNED sequences, return a 'position list' for the second
    sequence based on the parental sequence"""
    idx = 1
    numlist = []
    insertchars = "abcdefghijklmnopqrstuvwxyz"
    insertidx = 0
    for s1, s2 in zip(aseq1, aseq2):
        if s2 == "-":
            idx += 1
            continue
        if s1 == "-":
            numlist.append(str(idx - 1) + insertchars[insertidx % len(insertchars)])
            insertidx += 1
            continue
        insertidx = 0
        numlist.append(str(idx))
        idx += 1
    return numlist


class ParasailAlignment:
    """Convenience class to wrap the results of a parasail alignment"""

    def __init__(self, result):
        self.cigar = result.cigar.decode
        if isinstance(self.cigar, bytes):
            self.cigar = self.cigar.decode()
        # confusing nomenclature is for consistency with scikit-bio
        # where "query" is the initial sequence
        self.target = result.query
        self.query = result.ref
        self.score = result.score
        self._mutations = None

    def _tuples_from_cigar(self):
        tuples = []
        length_stack = []
        for character in self.cigar:
            if character.isdigit():
                length_stack.append(character)
            else:
                tuples.append((int("".join(length_stack)), character))
                length_stack = []
        return tuples

    @property
    def cigar_tuple(self):
        if hasattr(self, "_cigar_tuple"):
            return self._cigar_tuple
        self._cigar_tuple = self._tuples_from_cigar()
        return self._cigar_tuple

    def __repr__(self):
        return str(self)

    def __str__(self):
        a = chunked_lines(self.aligned_target_sequence(), spacer="")
        b = chunked_lines(self.aligned_query_sequence(), spacer="")
        out = []
        for t, q in zip(a, b):
            out.append(q)
            out.append("".join(["*" if x != y else (" " if x == " " else "|") for x, y in zip(t, q)]))
            out.append(t + "\n")
        return "\n".join(out)

    def __iter__(self):
        yield self.aligned_query_sequence()
        yield self.aligned_target_sequence()

    def as_mutations(self, reference=None):
        from .mutations import MutationSet, _get_aligned_muts

        seq1, seq2 = self
        mutstring = "/".join(_get_aligned_muts(seq1, seq2))
        if reference is not None:
            return MutationSet(mutstring, parental_numbering(*align_seqs(reference, seq1)))
        return MutationSet(mutstring)

    def print_alignment(self, max_length=80):
        print(self.aligned_query_sequence() + "\n" + self.aligned_target_sequence())

    def aligned_query_sequence(self):
        return self._get_aligned_sequence(self.query, "I")

    def aligned_target_sequence(self):
        return self._get_aligned_sequence(self.target, "D")

    def _get_aligned_sequence(self, seq, gap_type, gap_char="-", eq_char="="):
        # assume zero based
        # gap_type is 'D' when returning aligned query sequence
        # gap_type is 'I' when returning aligned target sequence
        aligned_sequence = ""
        index = 0
        for length, symbol in self.cigar_tuple:
            if symbol in (eq_char, "X"):
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
            elif symbol == gap_type:
                aligned_sequence += gap_char * length
            elif symbol in ("D", "I"):
                end = length + index
                aligned_sequence += seq[index:end]
                index += length
        return aligned_sequence

    @classmethod
    def from_seqs(cls, query, target, **kwargs):
        return align_seqs(query, target, **kwargs)
