"""
mutation strings attempt to follow HGVS-nomenclature
http://varnomen.hgvs.org/recommendations/protein/
but for now, the assumption is generally for single-letter AA codes
whereas the HGVS prefers three letter codes.

examples (converted to single letter codes):
    (prefix below specifies the reference protein)

    Substitution: [prefix][amino_acid][position][new_amino_acid]
        S65T
    Deletion: [prefix][start_amino_acid[_stop_amino_acid]]del
        C76del      -> single amino acid deletion
        C76_G79del  -> (inclusive) deletion of 4 amino acids
    Insertion: [prefix][previous-amino-acid_next-amino-acid]ins[inserted_seq]
        K23_L24insRSG   -> insert RSG in between K23 and L 24
        note: insertions must specify the full range.  "K23insRSG" is not allowed
    Deletion-insertion: [prefix][start_amino_acid[_stop_amino_acid]]delins[inserted_seq]
        C76delinsRRGY     -> delete one amino acid C76 and insert RRGY
        C76_G78delinsRRGY -> delete three amino acids (C76-G78) and insert RRGY
    Extension:
        for N Terminal: *[stop_codon_position][new_amino_acid]ext[]
            *315TextAKGT  -> add TAKGT to the end of the amino acid sequence
        for C Terminal:
            mostly this will be treated as an insertion...
            eg MLIK -> MVKSGEELIK = M1_L2insVKSGEE

full string:
    'S65T/C76del/C76_G79del/K23_L24insRSG/C76delinsRRGY/C76_G78delinsRRGY/*315TextAKGT/M1_L2insVKSGEE'
"""

import re
from .skbio_protein import SkbSequence
from .align import align_seqs


# optional prefix could be added
# (?:(?P<prefix>[A-Za-z0-9]+)\.)?
mutpattern = re.compile(r"""
    (?P<start_char>[{0}*]{1})
    (?P<start_idx>\d+)
    (?:_(?P<stop_char>[{0}]{1})(?P<stop_idx>\d+))?
    (?P<operation>delins|del|ins|)?
    (?P<new_chars>[A-Z*]+)?
    (?:ext(?P<ext>[{0}]+))?
    (?:,|$|/|\s)""".format("".join(SkbSequence.definite_chars), '{1}'), re.X)


def parse_mutstring(string):
    return [Mutation(*mut) for mut in mutpattern.findall(string)]


class Mutation(object):

    def __init__(self, start_char, start_idx, stop_char, stop_idx,
                 operation, new_chars, ext, idx0=1):
        self.start_char = start_char
        try:
            self.start_idx = int(start_idx)
        except ValueError:
            raise ValueError('Mutation must have integer start index')
        self.stop_char = stop_char
        try:
            self.stop_idx = int(stop_idx)
        except ValueError:
            self.stop_idx = None
        self.operation = 'ext' if ext else (operation or 'sub')
        assert self.operation in ('sub', 'del', 'ins', 'delins', 'ext'), 'Unrecognized operation: {}'.format(self.operation)
        if self.operation == 'sub' and (stop_char or self.stop_idx):
            raise ValueError('Substitution mutations cannot specify a range (or a stop_char/idx)')
        if self.operation == 'ins' and not (stop_char and self.stop_idx):
            raise ValueError('Insertion mutations must specify a range (with stop_char/idx')
        if self.operation == 'del' and new_chars:
            raise ValueError('Deletion mutations cannot specify new_chars (use delins instead)')
        if self.operation == 'sub' and len(new_chars) != 1:
            raise ValueError('A substitution must have a single new character')
        self.new_chars = new_chars
        self.ext = ext
        self.idx0 = idx0

    def __str__(self):
        out = '{}{}'.format(self.start_char, self.start_idx)
        if self.stop_char and self.stop_idx:
            out += '_{}{}'.format(self.stop_char, self.stop_idx)
        if self.operation == 'ext':
            out += self.new_chars + self.operation + self.ext
        elif self.operation == 'sub':
            out += self.new_chars
        else:
            out += self.operation + self.new_chars
        return out

    def __repr__(self):
        return '<Mutation: {}>'.format(self)

    def __eq__(self, other):
        for attr in ('start_idx', 'start_char', 'stop_idx',
                     'stop_char', 'operation', 'ext', 'new_chars'):
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __hash__(self):
        return hash((self.start_idx, self.start_char, self.stop_idx,
                     self.stop_char, self.operation, self.ext, self.new_chars))

    class SequenceMismatch(Exception):
        pass

    def _assert_position_consistency(self, seq, idx0=1):
        startpos = self.start_idx - idx0
        # leaving start_char out should prevent this check
        if self.operation == 'ext':
            if self.start_idx != len(seq) - idx0 + 2:
                raise ValueError('Extension start char not consistent')
            return
        if self.start_char and seq[startpos] != self.start_char:
            raise self.SequenceMismatch(
                'Mutation {} starting at {}{} does not match the sequence provided: {} (with idx0={})'
                .format(self, self.start_char, self.start_idx,
                        '{}>{}<{}'.format(seq[startpos - 3:startpos],
                                          seq[startpos],
                                          seq[startpos + 1:startpos + 4]),
                        idx0))
        if self.stop_idx and self.stop_char:
            stoppos = self.stop_idx - idx0
            if not seq[stoppos] == self.stop_char:
                raise self.SequenceMismatch(
                    'Mutation {} stopping at {}{} does not match the sequence provided: {} (with idx0={})'
                    .format(self, self.stop_char, self.stop_idx,
                            '{}>{}<{}'.format(seq[stoppos - 3:stoppos],
                                              seq[stoppos],
                                              seq[stoppos + 1:stoppos + 4]),
                            idx0))

    def __call__(self, seq, idx0=1):
        # not calling this since the start-CHAR may have changed from a
        # previous mutation when stringing...
        # self._assert_position_consistency(seq, idx0)
        startpos = self.start_idx - idx0
        if self.operation == 'sub':
            return (seq[:startpos] + self.new_chars + seq[startpos + 1:], 0)
        if self.operation == 'ins':
            return (seq[:startpos + 1] + self.new_chars + seq[startpos + 1:], len(self.new_chars))
        if self.operation in ('del', 'delins'):
            stoppos = startpos
            if self.stop_idx:
                stoppos = self.stop_idx - idx0
            shift = len(self.new_chars) - (stoppos + 1 - startpos)
            return (seq[:startpos] + self.new_chars + seq[stoppos + 1:], shift)
        if self.operation == 'ext':
            return seq + self.new_chars + self.ext, 0

    @classmethod
    def from_str(cls, mutstring, sep='/'):
        m = parse_mutstring(mutstring)
        assert len(m) == 1, 'For multiple mutations, create a MutationSet instead'
        return m[0]


def _get_aligned_muts(AQS, ATS, gapchars='-.', zeroindex=1):
    """ starting with two sequences that have been aligned,
    returns a list of mutation string codes, such as:
    ['K2_K6del', 'D26_N28delinsR', E39_A40insGD', 'R44K', *56LextPVPW']
    """
    out = []
    ins_start_idx = None
    insertions = ''
    delstart = None
    lastchar = '*'
    ishift = 0
    numdel = 0
    delins = ''

    def clear_insertions(ins_start_idx, insertions, extension=False):
        if extension:
            out.append('*{}{}ext{}'.format(
                ins_start_idx + 1, insertions[0], insertions[1:]))
        else:
            out.append('{}{}_{}{}ins{}'.format(
                ins_start_char, ins_start_idx, before,
                ins_start_idx + 1, insertions))

    def clear_deletions(delstart, numdel, delins, idx):
        string = '{}{}del'.format(
            delstart,
            '_{}'.format(lastchar + str(idx - 1)) if numdel > 1 else '')
        if delins:
            string += 'ins' + delins
        out.append(string)

    for idx, (before, after) in enumerate(zip(AQS, ATS), zeroindex):
        idx += ishift
        if before in gapchars:  # we have an insertion
            insertions += after  # add AA to insertion string
            if ins_start_idx is None:  # this is a new insertion
                ins_start_idx = idx - 1
                ins_start_char = lastchar
            ishift -= 1
            lastchar = before
            continue
        if after in gapchars:  # we have a deletion
            if numdel == 0:
                delstart = before + str(idx)
            numdel += 1
            lastchar = before
            continue

        if insertions:  # not an insertions but insertions need to be processed
            clear_insertions(ins_start_idx, insertions)
            ins_start_idx = None
            insertions = ''

        if before != after:
            if delstart:  # we are mid deletion... this should be a delins
                delins += after
                numdel += 1
                lastchar = before
                continue
            else:
                # regular single substitution
                out.append("".join([str(n) for n in [before, idx, after]]))
        else:
            if delstart:  # not a deletion but deletions need to be processed
                delstart = clear_deletions(delstart, numdel, delins, idx)
                delstart, numdel, delins = (None, 0, '')

        lastchar = before

    if insertions:
        clear_insertions(ins_start_idx, insertions, extension=True)
    if delstart:
        clear_deletions(delstart, numdel, delins, idx + 1)

    return out


def find_mutations(seq1, seq2, **kwargs):
    """get complete mutation string for a pair of sequences"""
    algn = align_seqs(seq1, seq2, **kwargs)
    mutstring = "/".join(_get_aligned_muts(*algn))
    return MutationSet.from_str(mutstring)


def mutate_sequence(seq, mutstring):
    ms = MutationSet.from_str(mutstring)
    return ms(seq)


class MutationSet(object):

    def __init__(self, muts=None):
        if isinstance(muts, str):
            muts = parse_mutstring(muts)
        elif not isinstance(muts, (list, set, tuple)):
            raise ValueError('Mutations argument must be str, list, set, or tuple')
        if not all([isinstance(m, Mutation) for m in muts]):
            raise ValueError('All MutationSet items must be Mutation Instances')
        self.muts = set(muts)

    def __contains__(self, query):
        if isinstance(query, str):
            query = Mutation.from_str(query)
        return query in self.muts

    def __add__(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        newset = self.muts
        [newset.add(that) for that in other.muts]
        return MutationSet(newset)

    def __sub__(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        newset = self.muts
        [newset.remove(that) for that in other.muts]
        return MutationSet(newset)

    def __call__(self, seq, idx0=1):
        """ apply the full mutation set to a sequence """
        for mut in self.muts:
            mut._assert_position_consistency(seq, idx0)
        shift = idx0
        for mut in sorted(set(self.muts), key=lambda x: x.start_idx):
            seq, new_offset = mut(seq, shift)
            shift -= new_offset
        return seq

    def union(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.union(other.muts))

    def intersection(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.intersection(other.muts))

    def difference(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.difference(other.muts))

    def issubset(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.issubset(other.muts)

    def issuperset(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.issuperset(other.muts)

    def isdisjoint(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.isdisjoint(other.muts)

    def __iter__(self):
        return self.muts.__iter__()

    def __eq__(self, other):
        """ should be rather robust way to compare mutations to some string
        or other Mutations instance, even if there's a little offset between
        the positions"""
        otherm = False
        if isinstance(other, MutationSet):
            otherm = other.muts
        elif isinstance(other, str):
            otherm = MutationSet.from_str(other).muts
        elif isinstance(other, (set, list, tuple)):
            try:
                otherm = MutationSet(other).muts
            except Exception:
                raise ValueError(
                    'Could not compare MutationSet object with other: {}'.format(other))
        if not otherm:
            raise ValueError(
                'operation not valid between type MutationSet and {}'.format(type(other)))
        else:
            if self.muts == otherm:
                return True
            else:
                if len(self.muts) == len(otherm):
                    if len(self.detect_offset(other)) == 1:
                        return True
        return False

    def __hash__(self):
        return hash(tuple(sorted(self.muts, key=lambda x: x.start_idx)))

    def __len__(self):
        return len(self.muts)

    def __repr__(self):
        return '<MutationSet: {}>'.format(str(self))

    def __str__(self):
        delim = '/'
        return delim.join([str(m) for m in sorted(set(self.muts), key=lambda x: x.start_idx)])

    @classmethod
    def from_str(cls, mutstring, sep='/'):
        return cls(parse_mutstring(mutstring))

    @property
    def deletions(self):
        return MutationSet(i for i in self.muts if i.operation == 'del')

    @property
    def insertions(self):
        return MutationSet(i for i in self.muts if i.operation == 'ins')

    @property
    def delinsertions(self):
        return MutationSet(i for i in self.muts if i.operation == 'delins')

    @property
    def substitutions(self):
        return MutationSet(i for i in self.muts if i.operation == 'sub')

    def detect_offset(self, refseq, maxshift=20, idx0=1):
        """ looks for a probable equality with frame shift between
        a sequence and this mutation set

        returns offset if there is a match, otherwise None
        """

        offsets = [0]
        for i in range(1, maxshift + 1):
            offsets.extend([i, -i])
        seqlist = list(str(refseq))

        mutD = {m.start_idx: m.start_char for m in self.muts}
        for offset in offsets:
            try:
                if all([seqlist[k - idx0 + offset] == v for k, v in mutD.items()]):
                    return offset
            except Exception as e:
                print(e)
                continue
        return None


def rand_mut(seq):
    from random import randint, choices

    operation = choices(['sub', 'sub', 'sub', 'sub', 'del', 'del', 'ins',
                         'ins', 'ins', 'delins', 'delins', 'ext'])[0]
    AAs = 'ACDEFGHIKLMNPQRSTVWY'
    start_idx = randint(1, len(seq) - 1)
    start_char = seq[start_idx - 1]
    ext, stop_idx, stop_char, new_chars = ('', '', '', '')
    if operation in ('sub', 'ext'):
        new_chars = choices(AAs)[0]
        if operation == 'ext':
            start_char = '*'
            start_idx = len(seq) + 1
            ext = "".join(choices(AAs, k=randint(1, 6)))
    elif operation in ('ins', 'delins'):
        new_chars = "".join(choices(AAs, k=randint(1, 6)))
    if operation in ('del', 'ins', 'delins'):
        stop_idx = start_idx + randint(0, 6)
        while stop_idx > len(seq) - 2:
            stop_idx = start_idx + randint(0, 6)
        stop_char = seq[stop_idx - 1]
    return Mutation(start_char, start_idx, stop_char, stop_idx, operation, new_chars, ext)
