import csv
import re
from proteins.models import Protein
from skbio import sequence as skbseq
from skbio.util import classproperty
from skbio.alignment import local_pairwise_align_ssw, make_identity_substitution_matrix

import skbio


def mustring_to_list(mutstring):
    aa_alph = "[%s]" % "".join(skbseq.Protein.definite_chars)
    return re.findall(r'(?P<pre>{0}+)(?P<pos>\d+)(?P<post>{0}+)'
                      .format(aa_alph), mutstring)


class Mutations(object):

    def __init__(self, muts=None):
        if isinstance(muts, str):
            muts = mustring_to_list(muts)
        assert all([len(i) == 3 for i in muts]), 'All mutations items must have 3 elements'
        self.muts = set([(a, int(b), c) for a, b, c in muts]) or set()

    def __eq__(self, other):
        if isinstance(other, Mutations):
            return self.muts == other.muts
        elif isinstance(other, str):
            return self == Mutations.from_str(other)
        elif isinstance(other, (set, list, tuple)):
            try:
                other = Mutations(other)
                return self.muts == other.muts
            except Exception:
                raise ValueError(
                    'Could not compare Mutations object with other: {}'.format(other))
        raise ValueError(
            'operation not valid between type Mutations and {}'.format(type(other)))

    def __len__(self):
        return len(self.muts)

    def __repr__(self):
        return '<Mutations({})>'.format(repr(self.muts))

    def __str__(self):
        joiner = '/'
        out = []
        dels = 0
        tmpi = 0
        for n, (before, idx, after) in enumerate(sorted(self.muts, key=lambda x: int(x[1]))):
            if after == '-':  # deletion
                if dels and idx != tmpi:
                    out.append(temp + 'del')
                    tmpi = 0
                    dels = 0
                if not dels:
                    temp = str(before) + str(idx)
                    tmpi = idx
                dels += 1
                tmpi += 1
            else:
                if dels:
                    if dels > 1:
                        temp += '_' + self.muts[n - 1][0] + str(self.muts[n - 1][1])
                    out.append(temp + 'del')
                    dels = 0
                    tmpi = 0
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


def get_mutations(seq1, seq2, **kwargs):
    if not isinstance(seq1, FPSeq):
        if isinstance(seq1, Protein):
            seq1 = seq1.seq
        seq1 = FPSeq(seq1)
    if not isinstance(seq2, skbseq.Protein):
        if isinstance(seq2, Protein):
            if not seq2.seq:
                raise ValueError('{} has no sequence data'.format(seq2))
            seq2 = seq2.seq
        seq2 = skbseq.Protein(seq2)
    algn, score, startend = seq1.align_to(seq2, **kwargs)
    offset = startend[0][0] + 1
    muts = set()
    for i in range(algn.shape[1]):
        if algn[0][i] != algn[1][i]:
            muts.add((str(algn[0][i]), i + offset, str(algn[1][i])))
    return Mutations(muts)


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


class FPSeq(skbseq.GrammaredSequence):

    def __init__(self, seq, *args, **kwargs):
        if isinstance(seq, Protein):
            seq = seq.seq
        super().__init__(seq, *args, **kwargs)

    @classproperty
    def degenerate_map(cls):
        return {
            "B": set("DN"), "Z": set("EQ"),
            "X": set("ACDEFGHIKLMNPQRSTVWY")
        }

    @classproperty
    def definite_chars(cls):
        return set("ACDEFGHIKLMNPQRSTVWY")

    @classproperty
    def default_gap_char(cls):
        return '-'

    @classproperty
    def gap_chars(cls):
        return set('-.')

    def same_as(self, other):
        return str(self) == str(other)

    def align_to(self, other, gop=4, gep=1):
        mtx = make_identity_substitution_matrix(1, -1, skbio.Protein.alphabet)
        if isinstance(other, FPSeq):
            other = skbseq.Protein(other.values)
        elif isinstance(other, str):
            other = skbseq.Protein(other)
        elif not isinstance(other, skbseq.Protein):
            raise ValueError('other must be either a str or subclass of skbio.Protein')
        this = skbseq.Protein(self.values)
        return local_pairwise_align_ssw(this, other, protein=True,
                                        substitution_matrix=mtx,
                                        gap_open_penalty=gop,
                                        gap_extend_penalty=gep)

    def mutations_to(self, other, **kwargs):
        return get_mutations(self, other, **kwargs)

    def mutations_from(self, other, **kwargs):
        return get_mutations(other, self, **kwargs)


# re.sub(r' \([A-Z][0-9_]+[A-Z]\)', '', name)


def getname(name):
    queries = [
        {'name__iexact': name},
        {'aliases__icontains': name},
        # {'name__icontains': name},
        {'name__iexact': re.sub(r' \((Before|Planar|wild).*', '', name)},
        {'aliases__icontains': re.sub(r' \((Before|Planar|wild).*', '', name)},
        # {'name__icontains': re.sub(r' \((Before|Planar).*', '', name)},
        {'name__iexact': name.strip('1')},
        # {'name__icontains': name.strip('1')},
    ]
    for query in queries:
        try:
            return Protein.objects.get(**query)
        except Exception:
            pass
    return None


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
