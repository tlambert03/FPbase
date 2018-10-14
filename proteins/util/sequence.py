import csv
import re
import textwrap
try:
    from skbio.sequence import GrammaredSequence
    from skbio.sequence import Protein as skProtein
    from skbio.util import classproperty
    from skbio.alignment import StripedSmithWaterman
    SKB = True
except ImportError:
    print('ERROR!! NO scikit-bio, things will break')
    SKB = False
try:
    import parasail as _parasail
except ImportError:
    print('ERROR!!! could not import parasail... will not be able to align')

# from skbio.alignment import  make_identity_substitution_matrix
# EYE_MTX = make_identity_substitution_matrix(1, -1, skProtein.alphabet)

AA_WEIGHTS = {
    'A': 89.0932,
    'C': 121.1582,
    'D': 133.1027,
    'E': 147.1293,
    'F': 165.1891,
    'G': 75.0666,
    'H': 155.1546,
    'I': 131.1729,
    'K': 146.1876,
    'L': 131.1729,
    'M': 149.2113,
    'N': 132.1179,
    'O': 255.3134,
    'P': 115.1305,
    'Q': 146.1445,
    'R': 174.201,
    'S': 105.0926,
    'T': 119.1192,
    'U': 168.0532,
    'V': 117.1463,
    'W': 204.2252,
    'Y': 181.1885
}
WATER = 18.0153


def protein_weight(seq):
    try:
        return sum(AA_WEIGHTS[x] for x in seq) - (len(seq) - 1) * WATER
    except KeyError as e:
        raise ValueError('%s is not a valid unambiguous amino acid letter' % e)


def nw_align(query, target, gop=2, gep=1, band_size=0):
    if isinstance(query, (str, FPSeq)):
        query = str(query)
    else:
        raise ValueError('query sequence must be str or FPseq')
    if isinstance(target, (str, FPSeq)):
        target = str(target)
    else:
        raise ValueError('query sequence must be str or FPseq')
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
        off = 1
        AQS, AQT = self.aligned_query_sequence(), self.aligned_target_sequence()
        muts = set((x, off + i, y) for x, (i, y) in zip(AQS, enumerate(AQT)) if x != y)
        self._mutations = Mutations(muts)
        return self._mutations

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
    # FIX ME: doesn't deal with deletions and insertions
    aa_alph = "[%s]" % "".join(skProtein.definite_chars)
    return re.findall(r'(?P<pre>{0}+)(?P<pos>\d+)(?P<post>{0}+)'
                      .format(aa_alph), mutstring)


if not SKB:
    class FPSeq(str):
        pass
else:
    class FPSeq(GrammaredSequence):

        def __eq__(self, other):
            return str(self) == str(other)

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

        @property
        def weight(self):
            try:
                return protein_weight(str(self)) / 1000
            except ValueError:
                pass

        def same_as(self, other):
            return str(self) == str(other)

        def align_to(self, other, **kwargs):
            return nw_align(self, other, **kwargs)

        def mutations_to(self, other, **kwargs):
            return self.align_to(other, **kwargs).mutations

        def mutate(self, mutations, first_index=1):
            out = list(str(self))
            if isinstance(mutations, str):
                mutations = Mutations(mutations)
            for before, index, after in mutations:
                if before and out[index - first_index] != before:
                    raise ValueError('Requested mutation {0}{1}{2} inconsistent with current sequence: {3}{1}'
                                     .format(before, index, after, out[index - first_index]))
                else:
                    out[index - first_index] = after
            return "".join(out)


class Mutations(object):

    def __init__(self, muts=None):
        if isinstance(muts, str):
            muts = mustring_to_list(muts)
        elif not isinstance(muts, (list, set, tuple)):
            raise ValueError('Mutations argument must be str, list, set, or tuple')
        if any([len(i) != 3 for i in muts]):
            raise ValueError('All mutations items must have 3 elements')
        self.muts = set([(a, int(b), c) for a, b, c in muts]) or set()

    def __add__(self, offset):
        assert isinstance(offset, int), 'offset must be an integer'
        return {(a, b + offset, d) for a, b, d in self.muts}

    def __sub__(self, offset):
        assert isinstance(offset, int), 'offset must be an integer'
        return {(a, b - offset, d) for a, b, d in self.muts}

    def __iter__(self):
        return self.muts.__iter__()

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

def get_gb_data(file):
    with open(file, 'r') as handle:
        text = handle.read()
    q = ('DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE')
    pat = ''
    for i in range(len(q) - 1):
        pat += r'(?=.*{} (?P<{}>.+){})?'.format(q[i], q[i].lower().replace(' ', ''), q[i + 1])
    pat += r'(?=.*COMMENT (?P<comment>.+)FEAT)?'
    pat += r'(?=.*PUBMED\s+(?P<pub>.+)REF)?'
    pat += r'(?=.*/translation="(?P<seq>.+)")?'
    D = re.search(pat, text, re.DOTALL).groupdict()
    for k, v in D.items():
        if v:
            D[k] = v.replace('\n', '').replace('  ', '').strip()
    D['seq'] = D.get('seq', '').replace(' ', '')
    D['name'] = os.path.basename(file).strip('.gb')
    return D

