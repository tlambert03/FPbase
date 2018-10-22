try:
    import parasail
except ImportError:
    print('ERROR!!! could not import parasail... will not be able to align')

import textwrap


def nw_align(query, target, gop=5, gep=1, band_size=0):
    """ basic parasail global alignment of two sequences
    result is wrapped in ParasailAlignment Class
    """
    if band_size:
        result = parasail.nw_banded(target, query, gop, gep,
                                    band_size, parasail.blosum62)
    else:
        result = parasail.nw_trace_scan_sat(target, query, gop, gep, parasail.blosum62)
    return ParasailAlignment(result)


def align_seqs(seq1, seq2):
    """take two strings and return two strings aligned to each other
    with gap chars if necessary"""
    algn = nw_align(str(seq1), str(seq2))
    return algn.aligned_query_sequence(), algn.aligned_target_sequence()


# @property
# def mutations(self):
#     if self._mutations is not None:
#         return self._mutations
#     off = 1
#     AQS, AQT = self.aligned_query_sequence(), self.aligned_target_sequence()
#     muts = set((x, off + i, y) for x, (i, y) in zip(AQS, enumerate(AQT)) if x != y)
#     self._mutations = Mutations(muts)
#     return self._mutations


class ParasailAlignment:

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
        if hasattr(self, '_cigar_tuple'):
            return self._cigar_tuple
        self._cigar_tuple = self._tuples_from_cigar()
        return self._cigar_tuple

    def __str__(self, width=70):
        a = textwrap.wrap(self.aligned_target_sequence(), width)
        b = textwrap.wrap(self.aligned_query_sequence(), width)
        out = []
        for t, q in zip(a, b):
            out.append(q)
            out.append("".join(['|' if x == y else '*' for x, y in zip(t, q)]))
            out.append(t + '\n')
        return "\n".join(out)

    def print_alignment(self, max_length=80):
        print(self.aligned_query_sequence() + '\n' + self.aligned_target_sequence())

    def aligned_query_sequence(self):
        return self._get_aligned_sequence(self.query, 'I')

    def aligned_target_sequence(self):
        return self._get_aligned_sequence(self.target, 'D')

    def _get_aligned_sequence(self, seq, gap_type, gap_char='-', eq_char='='):
        # assume zero based
        # gap_type is 'D' when returning aligned query sequence
        # gap_type is 'I' when returning aligned target sequence
        aligned_sequence = ''
        index = 0
        for length, symbol in self.cigar_tuple:
            if symbol in (eq_char, 'X'):
                aligned_sequence += seq[index:length + index]
                index += length
            elif symbol == gap_type:
                aligned_sequence += gap_char * length
            elif symbol in ('D', 'I'):
                aligned_sequence += seq[index:length + index]
                index += length
        return aligned_sequence

    @classmethod
    def from_seqs(cls, query, target, **kwargs):
        return nw_align(query, target, **kwargs)
