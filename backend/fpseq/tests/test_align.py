"""Tests for sequence alignment functionality.

These tests validate alignment behavior without depending on specific
alignment library implementation details.
"""

from __future__ import annotations

import pytest

from fpseq.align import align_seqs, parental_numbering

# Sample protein sequences for testing
# GFP sequence (partial)
GFP_SEQ = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCF"
    "SRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEY"
    "NYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKR"
    "DHMVLLEFVTAAGITHGMDELYK"
)

# mCherry sequence (partial)
MCHERRY_SEQ = (
    "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFM"
    "YGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQ"
    "KKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYT"
    "IVEQYERAEGRHSTGGMDELYK"
)

# Shorter sequences for testing edge cases
SHORT_SEQ1 = "MSKGEELFTGVVPILVELDGDVNGHKF"
SHORT_SEQ2 = "MVSKGEEDNMAIIKEFMRFKVHMEGSVNG"


class TestAlignSeqs:
    """Test the align_seqs function."""

    def test_identical_sequences(self):
        """Test alignment of identical sequences."""
        result = align_seqs(GFP_SEQ, GFP_SEQ)
        aligned_query, aligned_target = result

        # Aligned sequences should be identical
        assert aligned_query == aligned_target
        assert aligned_query == GFP_SEQ
        # Should have no gaps
        assert "-" not in aligned_query
        assert "-" not in aligned_target

    def test_different_sequences(self):
        """Test alignment of different sequences."""
        result = align_seqs(GFP_SEQ, MCHERRY_SEQ)
        aligned_query, aligned_target = result

        # Both aligned sequences should have the same length
        assert len(aligned_query) == len(aligned_target)
        # Should have some matches (these are both fluorescent proteins)
        matches = sum(1 for a, b in zip(aligned_query, aligned_target) if a == b and a != "-")
        assert matches > 50  # Should have significant similarity

    def test_short_sequences(self):
        """Test alignment of shorter sequences."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        aligned_query, aligned_target = result

        assert len(aligned_query) == len(aligned_target)
        # Both sequences start with M, so should have at least one match
        matches = sum(1 for a, b in zip(aligned_query, aligned_target) if a == b and a != "-")
        assert matches > 0

    def test_alignment_score(self):
        """Test that alignment has a score attribute."""
        result = align_seqs(GFP_SEQ, MCHERRY_SEQ)
        assert hasattr(result, "score")
        assert isinstance(result.score, int | float)
        # Score for similar sequences should be positive with BLOSUM62
        assert result.score > 0

    def test_alignment_cigar(self):
        """Test that alignment has CIGAR string."""
        result = align_seqs(GFP_SEQ, GFP_SEQ)
        assert hasattr(result, "cigar")
        assert isinstance(result.cigar, str)
        # For identical sequences, should be all matches
        assert "=" in result.cigar or "M" in result.cigar

    def test_aligned_query_sequence(self):
        """Test the aligned_query_sequence method."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        aligned_query = result.aligned_query_sequence()
        assert isinstance(aligned_query, str)
        # Should contain the original sequence characters (possibly with gaps)
        for char in SHORT_SEQ1:
            assert char in aligned_query

    def test_aligned_target_sequence(self):
        """Test the aligned_target_sequence method."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        aligned_target = result.aligned_target_sequence()
        assert isinstance(aligned_target, str)
        # Should contain the original sequence characters (possibly with gaps)
        for char in SHORT_SEQ2:
            assert char in aligned_target

    def test_string_representation(self):
        """Test string representation of alignment."""
        result = align_seqs(SHORT_SEQ1[:20], SHORT_SEQ2[:20])
        string_repr = str(result)
        assert isinstance(string_repr, str)
        # Should contain both sequences
        assert len(string_repr) > 0

    def test_iteration(self):
        """Test that alignment is iterable."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        aligned_query, aligned_target = result
        assert isinstance(aligned_query, str)
        assert isinstance(aligned_target, str)
        assert len(aligned_query) == len(aligned_target)

    def test_as_mutations(self):
        """Test conversion to mutations."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        mutations = result.as_mutations()
        # Should return a MutationSet
        assert hasattr(mutations, "muts")
        # Should have some mutations since sequences are different
        assert len(mutations) > 0

    def test_gap_penalties(self):
        """Test custom gap penalties."""
        # Default penalties
        result1 = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        # Higher gap open penalty
        result2 = align_seqs(SHORT_SEQ1, SHORT_SEQ2, gop=10)
        # Higher gap extension penalty
        result3 = align_seqs(SHORT_SEQ1, SHORT_SEQ2, gep=2)

        # All should produce valid alignments
        for result in [result1, result2, result3]:
            assert len(result.aligned_query_sequence()) == len(result.aligned_target_sequence())

    @pytest.mark.skip(reason="Banded alignment doesn't support traceback in current parasail implementation")
    def test_banded_alignment(self):
        """Test banded alignment."""
        # band_size > 0 should trigger banded alignment
        # Note: banded alignment may not support all features (like CIGAR)
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2, band_size=10)
        # Should still produce valid alignment even if CIGAR is not available
        assert hasattr(result, "query")
        assert hasattr(result, "target")

    def test_query_target_attributes(self):
        """Test that alignment has query and target attributes."""
        result = align_seqs(SHORT_SEQ1, SHORT_SEQ2)
        assert hasattr(result, "query")
        assert hasattr(result, "target")
        assert result.query == SHORT_SEQ1
        assert result.target == SHORT_SEQ2


class TestParentalNumbering:
    """Test the parental_numbering function."""

    def test_identical_sequences(self):
        """Test numbering with identical sequences."""
        aligned1 = "MSKGEELF"
        aligned2 = "MSKGEELF"
        numbering = parental_numbering(aligned1, aligned2)

        # Should be simple 1-based numbering
        assert numbering == ["1", "2", "3", "4", "5", "6", "7", "8"]

    def test_with_deletion(self):
        """Test numbering with deletion in child sequence."""
        aligned1 = "MSKGEELF"
        aligned2 = "MSK-EELF"  # deletion at position 4
        numbering = parental_numbering(aligned1, aligned2)

        # Should skip position 4
        assert numbering == ["1", "2", "3", "5", "6", "7", "8"]

    def test_with_insertion(self):
        """Test numbering with insertion in child sequence."""
        aligned1 = "MSK-EELF"
        aligned2 = "MSKGEELF"  # insertion after position 3
        numbering = parental_numbering(aligned1, aligned2)

        # Should have 3a for insertion
        assert numbering == ["1", "2", "3", "3a", "4", "5", "6", "7"]

    def test_multiple_insertions(self):
        """Test numbering with multiple consecutive insertions."""
        aligned1 = "MSK---EELF"
        aligned2 = "MSKABCEELF"  # three insertions after position 3
        numbering = parental_numbering(aligned1, aligned2)

        # Should have 3a, 3b, 3c for insertions
        assert numbering == ["1", "2", "3", "3a", "3b", "3c", "4", "5", "6", "7"]

    def test_mixed_operations(self):
        """Test numbering with both insertions and deletions."""
        aligned1 = "MSKGE-ELF"
        aligned2 = "MSK-EXELF"  # deletion at 4, insertion after 5
        numbering = parental_numbering(aligned1, aligned2)

        assert "5a" in numbering  # Should have insertion notation
        assert len(numbering) == len(aligned2.replace("-", ""))


class TestAlignmentIntegration:
    """Integration tests using align_seqs in realistic scenarios."""

    def test_mutation_detection(self):
        """Test using alignment to detect mutations."""
        # Create a sequence with known mutation
        original = "MSKGEELFTGVVPIL"
        mutated = "MSKGEELFTGVVAIL"  # P->A at position 13 (1-indexed)

        result = align_seqs(original, mutated)
        mutations = result.as_mutations()

        # Should detect the single substitution
        assert "P13A" in str(mutations)

    def test_deletion_detection(self):
        """Test detection of deletions."""
        original = "MSKGEELFTGVVPIL"
        deleted = "MSKGEELFTGVPIL"  # deletion of V at position 12

        result = align_seqs(original, deleted)
        mutations = result.as_mutations()

        # Should detect deletion
        assert "del" in str(mutations).lower()

    def test_insertion_detection(self):
        """Test detection of insertions."""
        original = "MSKGEELFTGVVPIL"
        inserted = "MSKGEELFTGVVAPIL"  # insertion of A

        result = align_seqs(original, inserted)
        mutations = result.as_mutations()

        # Should detect insertion
        assert "ins" in str(mutations).lower()

    def test_with_reference_sequence(self):
        """Test alignment with reference sequence for numbering."""
        seq1 = "MSKGEELFTGVVPIL"
        seq2 = "MSKGEELFTGVVAIL"
        reference = "MSKGEELFTGVVPIL"

        result = align_seqs(seq1, seq2)
        mutations = result.as_mutations(reference=reference)

        # Should return MutationSet with reference numbering
        assert hasattr(mutations, "muts")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
