import csv
import re
import textwrap
import json
import unicodedata
try:
    import requests
except ImportError:
    print('Could not import requests. Cannot pull sequences from fpbase')

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

letters_3to1 = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Gln': 'Q',
    'Glu': 'E',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Pyl': 'O',
    'Sec': 'U',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Xaa': 'X',
    'Xle': 'J'
}

letters_1to3 = {
    'A': 'Ala',
    'B': 'Asx',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'J': 'Xle',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'O': 'Pyl',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'U': 'Sec',
    'V': 'Val',
    'W': 'Trp',
    'X': 'Xaa',
    'Y': 'Tyr',
    'Z': 'Glx'
}


def seq1(seq):
    """Convert protein sequence from three-letter to one-letter code. """
    return ''.join(letters_3to1.get(aa) for aa in seq)


def seq3(seq):
    """Convert protein sequence from one-letter to three-letter code. """
    return ''.join(letters_1to3.get(aa) for aa in seq)


def protein_weight(seq):
    try:
        return sum(AA_WEIGHTS[x] for x in seq) - (len(seq) - 1) * WATER
    except KeyError as e:
        raise ValueError('%s is not a valid unambiguous amino acid letter' % e)


def slugify(value, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces to hyphens.
    Remove characters that aren't alphanumerics, underscores, or hyphens.
    Convert to lowercase. Also strip leading and trailing whitespace.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value).strip().lower()
    return re.sub(r'[-\s]+', '-', value)


def nw_align(query, target, gop=5, gep=1, band_size=0):
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


def mustring_to_list(mutstring):
    # FIX ME: doesn't deal with deletions and insertions
    aa_alph = "[%s]" % "".join(skProtein.definite_chars)
    return re.findall(r'(?P<pre>{0}+)(?P<pos>\d+)(?P<post>{0}+)'
                      .format(aa_alph), mutstring)


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
    # mock to allow heroku installation even if SKB fails
    # so that we can first install numpy ...
    class FPSeq(str):
        pass
else:
    class FPSeq(GrammaredSequence):

        def __init__(self, sequence, position_lables=None, **kwargs):
            if isinstance(sequence, str):
                sequence = sequence.replace(' ', '').replace('\n', '')
            if isinstance(position_lables, dict):
                pos_labels = list()
                i = 1
                for n in range(i, len(sequence) + i):
                    if n in position_lables:
                        pos_labels.append(str(position_lables[n]))
                    else:
                        pos_labels.append(str(i))
                        i += 1
                kwargs['positional_metadata'] = {'position_lables': pos_labels}
            super().__init__(sequence, **kwargs)

        def __eq__(self, other):
            if isinstance(other, str):
                other = other.replace(' ', '').replace('\n', '')
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

        def mutate(self, mutations, zeroindex=1, err_on_shift=False):
            if isinstance(mutations, str):
                mutations = Mutations(mutations)

            offset = detect_offset(self, mutations)
            if offset != 0:
                if offset is None or err_on_shift:
                    raise ValueError('Requested mutation numbers inconsistent '
                                     'with current sequence frame')
                else:
                    print('WARNING: mutations shifted {} position{} relative to sequence frame'
                          .format(offset, 's' if abs(offset) > 1 else ''))
            out = list(str(self))
            for before, index, after in mutations:
                    out[index - zeroindex + offset] = after
            return FPSeq("".join(out))

        @classmethod
        def from_fpbase(cls, slug):
            url = 'https://www.fpbase.org/api/{}/?format=json'.format(slugify(slug))
            response = requests.get(url)
            return FPSeq(json.loads(response.content).get('seq'))


def detect_offset(refseq, muts, maxshift=20, zeroindex=1):
    """ looks for a probable equality with frame shift between
    a sequence and a mutation set

    muts is a set of mutation 3-tuples
    returns offset if there is a match, otherwise None
    """
    if isinstance(muts, str):
        muts = mustring_to_list(muts)
    offsets = [0]
    for i in range(1, maxshift + 1):
        offsets.extend([i, -i])
    seqlist = list(str(refseq))
    mutD = {int(i): b for b, i, a in muts}
    for offset in offsets:
        try:
            if all([seqlist[k - zeroindex + offset] == v for k, v in mutD.items()]):
                return offset
        except Exception as e:
            print(e)
            continue
    return None


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

