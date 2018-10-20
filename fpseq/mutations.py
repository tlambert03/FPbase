import re
from .skbio_protein import SkbSequence


def detect_offset(refseq, muts, maxshift=20, zeroindex=1):
    """ looks for a probable equality with frame shift between
    a sequence and a mutation set

    muts is a set of mutation 3-tuples
    returns offset if there is a match, otherwise None
    """
    if isinstance(muts, str):
        muts = mutstring_to_list(muts)
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


def get_mutations(seq1, seq2, **kwargs):
    from .align import nw_align
    aligned = nw_align(seq1, seq2, **kwargs)
    off = aligned.query_begin + 1
    AQS, AQT = aligned.aligned_query_sequence, aligned.aligned_target_sequence
    muts = set((x, off + i, y) for x, (i, y) in zip(AQS, enumerate(AQT)) if x != y)
    return Mutations(muts)


def mutstring_to_list(mutstring):
    # FIX ME: doesn't deal with deletions and insertions
    aa_alph = "[%s]" % "".join(SkbSequence.definite_chars)
    return re.findall(r'(?P<pre>{0}+)(?P<pos>\d+)(?P<post>{0}+)'
                      .format(aa_alph), mutstring)


class Mutations(object):

    def __init__(self, muts=None):
        if isinstance(muts, str):
            muts = mutstring_to_list(muts)
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
