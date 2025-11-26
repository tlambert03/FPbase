"""Sequence alignment utilities.

This module provides functions for aligning protein sequences using
the Needleman-Wunsch global alignment algorithm with BLOSUM62 scoring matrix.
"""

from __future__ import annotations

from Bio import Align
from Bio.Align import substitution_matrices

from .util import chunked_lines


def align_seqs(query, target, gop=5, gep=1, band_size=0):
    """Perform global alignment of two sequences.

    Uses the Needleman-Wunsch algorithm with BLOSUM62 scoring matrix
    and affine gap penalties.

    Parameters
    ----------
    query : str
        First sequence to align
    target : str
        Second sequence to align
    gop : int, optional
        Gap open penalty (default: 5)
    gep : int, optional
        Gap extension penalty (default: 1)
    band_size : int, optional
        Not currently used (kept for backwards compatibility)

    Returns
    -------
    SequenceAlignment
        Alignment result wrapped in SequenceAlignment class
    """
    query = str(query)
    target = str(target)

    # Create aligner with BLOSUM62 and gap penalties
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # BioPython uses negative scores for penalties
    aligner.open_gap_score = -gop
    aligner.extend_gap_score = -gep

    # Perform alignment
    alignments = aligner.align(target, query)
    # Get the best alignment
    alignment = alignments[0]

    return SequenceAlignment(alignment, query, target)


def parental_numbering(aseq1, aseq2):
    """Generate position numbering for aligned sequences.

    Given two ALIGNED sequences, return a 'position list' for the second
    sequence based on the parental (first) sequence.

    Parameters
    ----------
    aseq1 : str
        Aligned parental sequence
    aseq2 : str
        Aligned child sequence

    Returns
    -------
    list[str]
        Position labels for each position in aseq2
    """
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


class SequenceAlignment:
    """Wrapper class for sequence alignment results.

    This class provides a consistent interface for alignment results,
    compatible with the previous parasail-based implementation.
    """

    def __init__(self, alignment, query, target):
        """Initialize from BioPython Alignment object.

        Parameters
        ----------
        alignment : Bio.Align.Alignment
            The BioPython alignment result
        query : str
            Original query sequence
        target : str
            Original target sequence
        """
        self._alignment = alignment
        # Store original sequences for compatibility
        self.query = query
        self.target = target
        self.score = alignment.score
        self._cigar = None
        self._mutations = None

    @property
    def cigar(self):
        """Generate CIGAR string from alignment coordinates.

        Returns
        -------
        str
            CIGAR string representation of the alignment
        """
        if self._cigar is not None:
            return self._cigar

        cigar_parts = []
        coords = self._alignment.coordinates

        # BioPython coordinates are [seq1_start, seq1_end] for each aligned region
        # coords[0] is target, coords[1] is query
        target_coords = coords[0]
        query_coords = coords[1]

        for i in range(len(target_coords) - 1):
            target_len = target_coords[i + 1] - target_coords[i]
            query_len = query_coords[i + 1] - query_coords[i]

            if target_len == query_len:
                # Match/mismatch region
                if target_len > 0:
                    # Use '=' for match region
                    # (Could be split into = for match and X for mismatch in future)
                    cigar_parts.append(f"{target_len}=")
            elif target_len > query_len:
                # Deletion from query (or insertion in target)
                cigar_parts.append(f"{target_len}D")
            else:
                # Insertion in query (or deletion from target)
                cigar_parts.append(f"{query_len}I")

        self._cigar = "".join(cigar_parts)
        return self._cigar

    @property
    def cigar_tuple(self):
        """Get CIGAR as list of (length, operation) tuples.

        Returns
        -------
        list[tuple[int, str]]
            List of (length, operation) tuples
        """
        if hasattr(self, "_cigar_tuple"):
            return self._cigar_tuple
        self._cigar_tuple = self._tuples_from_cigar()
        return self._cigar_tuple

    def _tuples_from_cigar(self):
        """Convert CIGAR string to list of tuples.

        Returns
        -------
        list[tuple[int, str]]
            List of (length, operation) tuples
        """
        tuples = []
        length_stack = []
        for character in self.cigar:
            if character.isdigit():
                length_stack.append(character)
            else:
                tuples.append((int("".join(length_stack)), character))
                length_stack = []
        return tuples

    def __repr__(self):
        """String representation."""
        return str(self)

    def __str__(self):
        """Pretty-print alignment with match indicators."""
        a = chunked_lines(self.aligned_target_sequence(), spacer="")
        b = chunked_lines(self.aligned_query_sequence(), spacer="")
        out = []
        for t, q in zip(a, b):
            out.append(q)
            out.append(
                "".join(["*" if x != y else (" " if x == " " else "|") for x, y in zip(t, q)])
            )
            out.append(t + "\n")
        return "\n".join(out)

    def __iter__(self):
        """Iterate over aligned sequences."""
        yield self.aligned_query_sequence()
        yield self.aligned_target_sequence()

    def as_mutations(self, reference=None):
        """Convert alignment to mutation set.

        Parameters
        ----------
        reference : str, optional
            Reference sequence for numbering

        Returns
        -------
        MutationSet
            Set of mutations between sequences
        """
        from .mutations import MutationSet, _get_aligned_muts

        seq1, seq2 = self
        mutstring = "/".join(_get_aligned_muts(seq1, seq2))
        if reference is not None:
            return MutationSet(mutstring, parental_numbering(*align_seqs(reference, seq1)))
        return MutationSet(mutstring)

    def print_alignment(self, max_length=80):
        """Print alignment to stdout.

        Parameters
        ----------
        max_length : int, optional
            Maximum line length (not currently used)
        """
        print(self.aligned_query_sequence() + "\n" + self.aligned_target_sequence())

    def aligned_query_sequence(self):
        """Get aligned query sequence with gaps.

        Returns
        -------
        str
            Query sequence with gaps inserted
        """
        # BioPython alignment format returns aligned sequences
        # alignment[0] is target (first arg to align()), alignment[1] is query
        return str(self._alignment[1])

    def aligned_target_sequence(self):
        """Get aligned target sequence with gaps.

        Returns
        -------
        str
            Target sequence with gaps inserted
        """
        # BioPython alignment format returns aligned sequences
        # alignment[0] is target (first arg to align()), alignment[1] is query
        return str(self._alignment[0])

    def _get_aligned_sequence(self, seq, gap_type, gap_char="-", eq_char="="):
        """Build aligned sequence from CIGAR string.

        Parameters
        ----------
        seq : str
            Original sequence
        gap_type : str
            'D' for query sequence, 'I' for target sequence
        gap_char : str, optional
            Character to use for gaps (default: '-')
        eq_char : str, optional
            Character representing matches in CIGAR (default: '=')

        Returns
        -------
        str
            Aligned sequence with gaps
        """
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
        """Create alignment from two sequences.

        Parameters
        ----------
        query : str
            First sequence
        target : str
            Second sequence
        **kwargs
            Additional arguments passed to align_seqs

        Returns
        -------
        SequenceAlignment
            Alignment result
        """
        return align_seqs(query, target, **kwargs)
