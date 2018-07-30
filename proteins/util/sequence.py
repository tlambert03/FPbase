import csv
import re
from skbio.sequence import GrammaredSequence
from skbio.sequence import Protein as skProtein
from skbio.util import classproperty
from skbio.alignment import StripedSmithWaterman
try:
    import parasail as _parasail
except ImportError:
    print('ERROR!!! could not import parasail... will not be able to align')

# from skbio.alignment import  make_identity_substitution_matrix
# EYE_MTX = make_identity_substitution_matrix(1, -1, skProtein.alphabet)


def nw_align(query, target, gop=2, gep=1, band_size=0):
    query = str(query)
    target = str(target)
    if band_size:
        result = _parasail.nw_banded(target, query, gop, gep,
                                     band_size, _parasail.blosum62)
    else:
        result = _parasail.nw_trace_scan_sat(target, query, gop, gep, _parasail.blosum62)
    return ParasailAlignment(result)


def get_mutations(seq1, seq2, **kwargs):
    if not isinstance(seq1, FPSeq):
        seq1 = FPSeq(seq1)
    aligned = seq1.align_to(seq2, **kwargs)
    off = aligned.query_begin + 1
    AQS, AQT = aligned.aligned_query_sequence, aligned.aligned_target_sequence
    muts = set((x, off + i, y) for x, (i, y) in zip(AQS, enumerate(AQT)) if x != y)
    return Mutations(muts)


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

    @property
    def mutations(self):
        if self._mutations is not None:
            return self._mutations
        off = 0
        AQS, AQT = self.aligned_query_sequence(), self.aligned_target_sequence()
        muts = set((x, off + i, y) for x, (i, y) in zip(AQS, enumerate(AQT)) if x != y)
        self._mutations = Mutations(muts)
        return self._mutations

    def print_alignment(self, max_length=80):
        print(a.aligned_query_sequence() + '\n' + a.aligned_target_sequence())

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
            if symbol == eq_char:
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


def get_mutations(seq1, seq2, tuple_cigar):
    muts = []
    index = 0
    for length, mid in tuple_cigar:
        if mid == 'X':
            aligned_sequence += sequence[index:length + index]
            index += length
        elif mid == gap_type:
            aligned_sequence += gap_char * length
        else:
            pass
    return aligned_sequence


def mustring_to_list(mutstring):
    aa_alph = "[%s]" % "".join(skProtein.definite_chars)
    return re.findall(r'(?P<pre>{0}+)(?P<pos>\d+)(?P<post>{0}+)'
                      .format(aa_alph), mutstring)


class FPSeq(GrammaredSequence):

    @classproperty
    def degenerate_map(cls):
        return skProtein.degenerate_map

    @classproperty
    def definite_chars(cls):
        return skProtein.definite_chars

    @classproperty
    def default_gap_char(cls):
        return '-'

    @classproperty
    def gap_chars(cls):
        return set('-.')

    def same_as(self, other):
        return str(self) == str(other)

    def align_to(self, other, gop=2, gep=1, **kwargs):
        ssw = StripedSmithWaterman(str(self), gop, gep, **kwargs)
        return ssw(str(other))

    def mutations_to(self, other, **kwargs):
        return get_mutations(self, other, **kwargs)

    def mutations_from(self, other, **kwargs):
        return get_mutations(other, self, **kwargs)


class Mutations(object):

    def __init__(self, muts=None):
        if isinstance(muts, str):
            muts = mustring_to_list(muts)
        elif not isinstance(muts, (list, set, tuple)):
            raise ValueError('Mutations argument must be str, list, set, or tuple')
        if any([len(i) != 3 for i in muts]):
            raise ValueError('All mutations items must have 3 elements')
        self.muts = set([(a, int(b), c) for a, b, c in muts]) or set()

    def __eq__(self, other):
        """ should be rather robust way to compare mutations to some string
        or other Mutations instance, even if there's a little offset between
        the positions"""
        otherm = False
        if isinstance(other, Mutations):
            otherm = other.muts
        elif isinstance(other, str):
            otherm = Mutations.from_str(other).muts
        elif isinstance(other, (set, list, tuple)):
            try:
                otherm = Mutations(other).muts
            except Exception:
                raise ValueError(
                    'Could not compare Mutations object with other: {}'.format(other))
        if not otherm:
            raise ValueError(
                'operation not valid between type Mutations and {}'.format(type(other)))
        else:
            if self.muts == otherm:
                return True
            else:
                if len(self.muts) == len(otherm):
                    if len(self.detect_offset(other)) == 1:
                        return True
        return False

    def __len__(self):
        return len(self.muts)

    def __repr__(self):
        return '<Mutations: {}>'.format(self)

    def __str__(self):
        joiner = '/'
        out = []
        dels = 0
        tmpi = 0
        for n, (before, idx, after) in enumerate(sorted(self.muts, key=lambda x: int(x[1]))):
            after = 'del' if after == '-' else after
            out.append("".join([str(n) for n in [before, idx, after]]))
        return joiner.join(out)

    @classmethod
    def from_str(cls, mutstring, sep='/'):
        return cls(mutstring_to_list(mutstring))

    @property
    def deletions(self):
        return [i for i in self.muts if i[2] == '-']

    @property
    def insertions(self):
        return [i for i in self.muts if i[0] == '-']

    def detect_offset(self, other):
        """ looks for a probable equality with frame shift between
        two mutation sets.

        returns a dict where the keys are the detected offsets, and the
        values are the number of times that offset was observed in the mutations
        keys with higher values are more likely to represent true frame shifts
        """

        from collections import Counter
        mymuts = {'%s%s' % (a, c): b for a, b, c in self.muts}
        theirmuts = {'%s%s' % (a, c): b for a, b, c in other.muts}
        matches = []
        for key, pos in mymuts.items():
            if key in theirmuts:
                matches.append((pos, theirmuts[key]))
        return dict(Counter([a - b for a, b in matches]))


def show_align(seq1, seq2, match=' ', mismatch='*', gap='-', **kwargs):
    algn, score, startend = align_prot(seq1, seq2, **kwargs)
    lendiff = startend[0][0] - startend[1][0]
    if lendiff < 0:
        print(' ' * abs(lendiff), end='')
    else:
        print(seq1[:abs(lendiff)], end='')
    print(str(algn[0]))
    print(' ' * abs(lendiff), end='')
    for i in range(algn.shape[1]):
        if '-' in (str(algn[0][i]), str(algn[1][i])):
            char = gap
        else:
            char = match if algn[0][i] == algn[1][i] else mismatch
        print(char, end='')
    print()
    if lendiff > 0:
        print(' ' * abs(lendiff), end='')
    else:
        print(seq1[:abs(lendiff)], end='')
    print(str(algn[1]))
    return algn


# re.sub(r' \([A-Z][0-9_]+[A-Z]\)', '', name)


def osfp_import():
    from collections import defaultdict
    with open('_data/osfp-full-data-set.csv') as f:
        csvrows = csv.reader(f)
        D = defaultdict(dict)
        for name, agg, seq, doi in csvrows:
            D[name]['seq'] = seq.replace('\n', '')
            D[name]['agg'] = agg
            D[name]['doi'] = doi
        return D
        """
        for name, agg, seq, doi in csvrows:
            seq = seq.replace('\n', '')
            p = getname(name)
            if not p:
                continue
            if not p.seq:
                print('ADD SEQ: ', name)
            elif p.seq == seq:
                pass
            else:
                print('seq mismatch: ', name)

            # DOI
            if p.primary_reference:
                if not doi == p.primary_reference.doi:
                    print('DOI mismatch on {} ({}->{})'.format(name, p.primary_reference.doi, doi))
            elif doi:
                print('Add DOI to {}: {}'.format(p, doi))

            # AGG
            if not agg == p.get_agg_display():
                if not p.agg == Protein.WEAK_DIMER:
                    print('change {} agg from {} to {}'.format(name, p.get_agg_display(), agg))
        """

def snapgene_import():
    with open('snapgene.csv') as f:
        csvrows = csv.reader(f)
        for row in csvrows:
            name = row[4]
            seq = row[6]
            try:
                p = Protein.objects.get(name=name)
            except Protein.DoesNotExist:
                continue
            if not p:
                continue
            if p.seq:
                if not ParasailAlignment.from_seqs(seq, p.seq).mutations:
                    continue
                print(name, ParasailAlignment.from_seqs(seq, p.seq).mutations)
