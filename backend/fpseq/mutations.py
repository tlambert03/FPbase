"""Mutations module

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
import warnings

import numpy as np

from .align import align_seqs, parental_numbering
from .skbio_protein import SkbSequence

# optional prefix could be added
# (?:(?P<prefix>[A-Za-z0-9]+)\.)?
mutpattern = re.compile(
    r"""
    (?P<start_char>[{0}*]{1})
    (?P<start_idx>\d+)
    (?:_(?P<stop_char>[{0}]{1})(?P<stop_idx>\d+))?
    (?P<operation>delins|del|ins|)?
    (?P<new_chars>[{0}*]+)?
    (?:ext(?P<ext>[{0}]+))?
    (?:,|$|/|\s)""".format(r"A-Z*", "{1}"),
    re.X,
)


NEUTRAL_MUTATIONS = ["K26R", "Q80R", "N146H", "H231L"]


def parse_mutstring(string):
    """Returns a list of Mutation objects found in a mutation string"""
    return [Mutation(*mut) for mut in mutpattern.findall(str(string))]


DEFAULT_ALPHABET = "".join(SkbSequence.definite_chars.union("X*"))


class Mutation:
    """Basic mutation object to represent a single mutation operation"""

    class SequenceMismatch(Exception):
        pass

    def __init__(
        self,
        start_char,
        start_idx,
        stop_char,
        stop_idx,
        operation,
        new_chars,
        ext,
        start_label=None,
        stop_label=None,
        idx0=1,
        alphabet=DEFAULT_ALPHABET,
    ):
        """
        start_char: single letter amino acid code for start position (e.g. A)
        start_idx:  position of start mutation (e.g. 206)
        stop_char:  single letter amino acid code for stop position (e.g. A)
        stop_idx:   range ending for mutations spanning multiple amino acids
        operation:  type of mutation.  Must be one of:
                    {'sub', 'del', 'ins', 'delins', 'ext'}
        new_chars:  single or multiple letter amino acid code(s) for new amino acids
        ext:        if operation == 'ext', extension characters to add (after new_chars)
        start_label:  optional string to label the start position (e.g. V1a)
                    this is used when the start character does not agree with the
                    start index... for instance, when displaying a mutation
                    with numbering relative to an ancestor other than parent.
        idx0:       index of first character (not well tested)

        """
        if alphabet and start_char not in alphabet:
            raise ValueError(f"Invalid Amino Acid code: {start_char}")
        self.start_char = start_char
        try:
            self.start_idx = int(start_idx)
        except ValueError as e:
            raise ValueError("Mutation must have integer start index") from e

        if alphabet and stop_char and stop_char not in alphabet:
            raise ValueError(f"Invalid Amino Acid code: {stop_char}")
        self.stop_char = stop_char
        try:
            self.stop_idx = int(stop_idx)
        except ValueError:
            self.stop_idx = None
        self.operation = "ext" if ext else (operation or "sub")
        if self.operation not in ("sub", "del", "ins", "delins", "ext"):
            raise ValueError(f"Unrecognized operation: {self.operation}")
        if self.operation == "sub" and (stop_char or self.stop_idx):
            raise ValueError("Substitution mutations cannot specify a range (or a stop_char/idx)")
        if stop_idx and (int(stop_idx) < int(start_idx)):
            raise ValueError(f"Stop position ({stop_idx}) must be greater than start position ({start_idx})")
        if self.operation.endswith("ins"):
            if not (stop_char and self.stop_idx):
                print(stop_char)
                print(self.stop_idx)
                raise ValueError("Insertion mutations must specify a range (with stop_char/idx)")
            if not len(new_chars):
                raise ValueError("Insertion mutations must specify new characters to insert")
            if self.operation == "ins" and (int(stop_idx) - int(start_idx) != 1):
                raise ValueError(
                    "Insertion range ({}-{}) is {} than 1 position".format(
                        start_idx,
                        stop_idx,
                        "greater" if (int(stop_idx) - int(start_idx) > 1) else "less",
                    )
                )
        if self.operation == "del" and new_chars:
            raise ValueError("Deletion mutations cannot specify new_chars (use delins instead)")
        self.new_chars = new_chars
        self.ext = ext
        self.idx0 = idx0
        self.start_label = start_label
        self.stop_label = stop_label
        if self.operation == "sub" and len(new_chars) != 1:
            raise ValueError(f"A substitution must have a single new character {self}")

    def __str__(self):
        out = f"{self.start_char}{self.start_label or self.start_idx}"
        if self.stop_char and self.stop_idx:
            out += f"_{self.stop_char}{self.stop_label or self.stop_idx}"
        if self.operation == "ext":
            out += self.new_chars + self.operation + self.ext
        elif self.operation == "sub":
            out += self.new_chars
        else:
            out += self.operation + self.new_chars
        return out

    def __repr__(self):
        return f"<Mutation: {self}>"

    def __eq__(self, other):
        if str(self) == str(other):
            return True
        return False

    def __hash__(self):
        return hash(str(self))

    def __call__(self, seq, idx0=1):
        # not calling this since the start-CHAR may have changed from a
        # previous mutation when stringing...
        # instead... this gets called at each step of MutationSet.apply
        # self._assert_position_consistency(seq, idx0)
        startpos = self.start_idx - idx0
        if startpos > len(seq):  # allowing one extra position for extensions
            raise IndexError(f"Starting position {self.start_idx} is outside of sequence with length {len(seq)}")
        if self.operation == "sub":
            end = startpos + 1
            return (seq[:startpos] + self.new_chars + seq[end:], 0)
        if self.operation == "ins":
            end = startpos + 1
            return (
                seq[: startpos + 1] + self.new_chars + seq[end:],
                len(self.new_chars),
            )
        if self.operation in ("del", "delins"):
            stoppos = startpos
            if self.stop_idx:
                stoppos = self.stop_idx - idx0
            nextpos = stoppos + 1
            shift = len(self.new_chars) - (nextpos - startpos)
            return (seq[:startpos] + self.new_chars + seq[nextpos:], shift)
        if self.operation == "ext":
            return seq + self.new_chars + self.ext, 0

    def _assert_position_consistency(self, seq, idx0=1):
        """test whether the mutation actually lines up with the
        provided sequence being mutated."""
        startpos = self.start_idx - idx0
        # leaving start_char out should prevent this check
        if self.operation == "ext":
            if self.start_idx != len(seq) - idx0 + 2:
                raise ValueError("Extension start char not consistent")
            return
        # if we're inserting at the very beginning, we don't need to check
        if startpos < 0:
            return
        if self.start_char and seq[startpos] != self.start_char:
            beg = startpos - 3
            beg2 = startpos + 1
            end = startpos + 4
            raise self.SequenceMismatch(
                "Mutation {} does not align with the parent seq: {}.".format(
                    self,
                    f"{seq[beg:startpos]}>{seq[startpos]}<{seq[beg2:end]}",
                )
            )
        if self.stop_idx and self.stop_char:
            stoppos = self.stop_idx - idx0
            if not seq[stoppos] == self.stop_char:
                beg = stoppos - 3
                beg2 = stoppos + 1
                end = stoppos + 4
                raise self.SequenceMismatch(
                    "Mutation {} does not match the sequence provided: {}".format(
                        self,
                        f"{seq[beg:stoppos]}>{seq[stoppos]}<{seq[beg2:end]}",
                    )
                )

    @classmethod
    def from_str(cls, mutstring, sep="/"):
        """Generate a Mutation object from a mutation string such as 'A206K'"""
        m = parse_mutstring(mutstring)
        if not m:
            raise ValueError(f"Mutation code invalid: {mutstring}")
        if len(m) > 1:
            raise ValueError("Multiple mutation codes found. For multiple mutations, create a MutationSet instead")
        return m[0]


def _get_aligned_muts(AQS, ATS, gapchars="-.", zeroindex=1):
    """starting with two sequences that have been aligned,
    returns a list of mutation string codes, such as:
    ['K2_K6del', 'D26_N28delinsR', E39_A40insGD', 'R44K', *56LextPVPW']
    """
    out = []
    ins_start_idx = None
    insertions = ""
    delstart = None
    lastchar = "*"
    ishift = 0
    numdel = 0
    delins = ""

    def clear_insertions(ins_start_idx, insertions, extension=False):
        if extension:
            out.append(f"*{ins_start_idx + 1}{insertions[0]}ext{insertions[1:]}")
        else:
            out.append(f"{ins_start_char}{ins_start_idx}_{before}{ins_start_idx + 1}ins{insertions}")

    def clear_deletions(delstart, numdel, delins, idx):
        string = "{}{}del".format(delstart, f"_{lastchar + str(idx - 1)}" if numdel > 1 else "")
        if delins:
            string += "ins" + delins
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
            insertions = ""

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
                delstart, numdel, delins = (None, 0, "")

        lastchar = before

    if insertions:
        clear_insertions(ins_start_idx, insertions, extension=True)
    if delstart:
        clear_deletions(delstart, numdel, delins, idx + 1)

    return out


def get_mutations(seq1, seq2, reference=None):
    """Detects mutations from seq1 to seq2, returns a MutationSet"""
    algn = align_seqs(seq1, seq2)
    return algn.as_mutations(reference=reference)


def mutate_sequence(seq, mutstring, idx0=1, correct_offset=False):
    """applies the provided mutstring to the provided sequence,
    returning a mutated string result"""
    ms = MutationSet.from_str(mutstring)
    return ms.apply(seq, idx0=idx0, correct_offset=correct_offset)


class MutationSet:
    """Class to hold a set of mutation objects, and apply them to a sequence.

    Mostyl a wrapper around a python set()
    """

    def __init__(self, muts=None, position_labels=None):
        """optional position_labels list will change the numbering of the
        muationset ... for instance, to match a reference sequence numbering"""
        if isinstance(muts, str):
            muts = parse_mutstring(muts)
        elif not isinstance(muts, list | set | tuple):
            raise ValueError("Mutations argument must be str, list, set, or tuple")
        if not all(isinstance(m, Mutation) for m in muts):
            raise ValueError("All MutationSet items must be Mutation Instances")
        self.muts = set(muts)
        if position_labels is not None:
            self.position_labels = position_labels
            for mut in self.muts:
                if len(position_labels) >= mut.start_idx:
                    # the - 1 assumes that the mutation positions are 1-indexed
                    mut.start_label = position_labels[mut.start_idx - 1]
                    if mut.stop_idx:
                        mut.stop_label = position_labels[mut.stop_idx - 1]
                else:
                    try:
                        mut.start_label = str(len(position_labels) - mut.start_idx + int(position_labels[-1]))
                    except Exception:
                        mut.start_label = str(mut.start_idx)

    def __contains__(self, query):
        if isinstance(query, str):
            query = Mutation.from_str(query)
        return query in self.muts

    def __add__(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        newset = self.muts.copy()
        [newset.add(that) for that in other.muts]
        return MutationSet(newset)

    def __sub__(self, other):
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        newset = self.muts.copy()
        [newset.discard(that) for that in other.muts]
        return MutationSet(newset)

    def apply(self, seq, idx0=1, correct_offset=False):
        """apply the full mutation set to a sequence"""

        shift = idx0
        for mut in self.muts:
            try:
                mut._assert_position_consistency(seq, shift)
            except Mutation.SequenceMismatch as e:
                offset = self.detect_offset(seq)
                if offset:
                    if correct_offset:
                        warnings.warn(
                            f"An offset of {offset} amino acids was detected"
                            " between the sequence and the mutation "
                            "set, and automatically corrected",
                            stacklevel=2,
                        )
                        shift -= offset
                    else:
                        raise Mutation.SequenceMismatch(
                            "{}. But a match was found {} position{} away: {}".format(
                                str(e),
                                offset,
                                "s" if abs(offset) > 1 else "",
                                self.shift(offset),
                            )
                        ) from e
                else:
                    raise e
        for mut in self:
            seq, new_offset = mut(seq, shift)
            shift -= new_offset
        if correct_offset:
            return seq, shift
        return seq

    def _has_adjacent(self):
        groups = self._consecutive_groups()
        for g in groups:
            ops = [m.operation for m in g]
            if "sub" in ops and ("del" in ops or "delins" in ops):
                return True
        return False

    def merge_delins(self, merge_subs=5):
        """Clean up mutation set to remove substitutions next to dels or delins"""
        groups = self._consecutive_groups()
        newgroups = []
        for g in groups:
            ops = [m.operation for m in g]
            if (len(g) >= merge_subs and all(m == "sub" for m in ops)) or (
                "sub" in ops and ("del" in ops or "delins" in ops)
            ):
                # should already be sorted
                start = f"{g[0].start_char}{g[0].start_idx}"
                if g[-1].stop_char:
                    stop = f"{g[-1].stop_char}{g[-1].stop_idx}"
                else:
                    stop = f"{g[-1].start_char}{g[-1].start_idx}"
                newchars = "".join([m.new_chars for m in g])
                newgroups.append(Mutation.from_str(f"{start}_{stop}delins{newchars}"))
            else:
                newgroups.extend(list(g))
        # can't decide whether to return new object or update this one
        self.muts = set(newgroups)

    def _consecutive_groups(self):
        """returns a list of np.arrays containin adjacent Mutation objects"""
        msl = list(self)
        g = np.split(msl, np.where(np.diff([m.start_idx for m in msl]) != 1)[0] + 1)
        return g

    def union(self, other):
        """Return the union of mutations as a new set.

        (i.e. all mutations that are in either set.)"""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.union(other.muts))

    def intersection(self, other):
        """Return the intersection of two mutation sets as a new set.

        (i.e. all mutations that are in both sets.)"""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.intersection(other.muts))

    def difference(self, other):
        """Return the difference of two or more MutationSets as a new set.

        (i.e. all mutations that are in this MutationSet but not the other.)"""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return MutationSet(self.muts.difference(other.muts))

    def issubset(self, other):
        """Report whether another MutationSet contains this MutationSet."""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.issubset(other.muts)

    def issuperset(self, other):
        """Report whether this MutationSet contains another MutationSet."""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.issuperset(other.muts)

    def isdisjoint(self, other):
        """Return True if two MutationSets have a null intersection"""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.isdisjoint(other.muts)

    def add(self, other):
        """Add a mutation to a MutationSet.

        This has no effect if the mutation is already present."""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.add(other.muts)

    def remove(self, other):
        """Remove a mutation from a MutationSet; it must be a member.

        If the mutation is not a member, raise a KeyError."""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.remove(other.muts)

    def discard(self, other):
        """Remove a mutation from a MutationSet if it is a member.

        If the mutation is not a member, do nothing."""
        if isinstance(other, str):
            other = MutationSet.from_str(other)
        return self.muts.discard(other.muts)

    def __iter__(self):
        sorted_muts = sorted(self.muts, key=lambda x: x.start_idx)
        while sorted_muts:
            yield sorted_muts.pop(0)

    def __eq__(self, other):
        """Determine whether two mutation sets are the same"""
        if not other:
            return False
        otherm = False
        if isinstance(other, MutationSet):
            otherm = other.muts
        elif isinstance(other, str):
            otherm = MutationSet.from_str(other).muts
        elif isinstance(other, set | list | tuple):
            try:
                otherm = MutationSet(other).muts
            except Exception as e:
                raise ValueError(f"Could not compare MutationSet object with other: {other}") from e
        if not otherm:
            raise ValueError(f"operation not valid between type MutationSet and {type(other)}")
        else:
            if self.muts == otherm:
                return True
            # else:
            #     if len(self.muts) == len(otherm):
            #         if len(self.detect_offset(other)) == 1:
            #             return True
        return False

    def __hash__(self):
        return hash(tuple(self))

    def __len__(self):
        return len(self.muts)

    def __repr__(self):
        return f"<MutationSet: {self}>"

    def __str__(self):
        delim = "/"
        return delim.join([str(m) for m in sorted(set(self.muts), key=lambda x: x.start_idx)])

    @classmethod
    def from_str(cls, mutstring, sep="/"):
        return cls(parse_mutstring(mutstring))

    @property
    def deletions(self):
        return MutationSet([i for i in self.muts if i.operation == "del"])

    @property
    def insertions(self):
        return MutationSet([i for i in self.muts if i.operation == "ins"])

    @property
    def delinsertions(self):
        return MutationSet([i for i in self.muts if i.operation == "delins"])

    @property
    def substitutions(self):
        return MutationSet([i for i in self.muts if i.operation == "sub"])

    @property
    def extensions(self):
        return MutationSet([i for i in self.muts if i.operation == "ext"])

    def shift(self, amount):
        """shift the position numbering of the mutation set by amount"""
        ms = MutationSet(str(self))
        for mut in ms.muts:
            mut.start_idx += amount
        return ms

    def detect_offset(self, refseq, maxshift=20, idx0=1):
        """looks for a probable equality with frame shift between
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
                if all(seqlist[pos - idx0 + offset] == letter for pos, letter in mutD.items()):
                    return offset
            except IndexError:
                continue
        return None

    def relative_to_root(self, parent, root):
        """display mutation string with parent amino acids, but with positioning
        relative to some other root sequence"""
        return str(MutationSet(str(self), parental_numbering(*align_seqs(root, parent))))


def rand_mut(seq):
    from random import choices, randint

    # make extensions less likely
    ch = ["sub"] * 10 + ["del"] * 5 + ["ins"] * 5 + ["delins"] * 3 + ["ext"]
    operation = choices(ch)[0]
    AAs = "ACDEFGHIKLMNPQRSTVWY"
    start_idx = randint(1, len(seq) - 1)
    start_char = seq[start_idx - 1]
    ext, stop_idx, stop_char, new_chars = ("", "", "", "")
    if operation in ("sub", "ext"):
        new_chars = choices(AAs)[0]
        if operation == "ext":
            start_char = "*"
            start_idx = len(seq) + 1
            ext = "".join(choices(AAs, k=randint(1, 6)))
    elif operation in ("ins", "delins"):
        new_chars = "".join(choices(AAs, k=randint(1, 6)))
    if operation in ("del", "ins", "delins"):
        if operation == "ins":
            stop_idx = start_idx + 1
        else:
            stop_idx = start_idx + randint(1, 6)
        while stop_idx > len(seq) - 2:
            stop_idx = start_idx + randint(0, 6)
        stop_char = seq[stop_idx - 1]
    return Mutation(
        start_char,
        start_idx,
        stop_char,
        stop_idx,
        operation,
        new_chars,
        ext,
        alphabet=None,
    )
