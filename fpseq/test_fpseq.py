import pytest

from . import Mutation, mutate_sequence


class TestMutations:
    def test_deletion_position_consistency(self):
        mut = Mutation.from_str("D2_E3del")
        try:
            mut._assert_position_consistency("MDELYK")
        except Exception:
            pytest.fail("Unexpected Exception in position_consistency ..")
        with pytest.raises(Mutation.SequenceMismatch):
            mut._assert_position_consistency("MADELYK")
        with pytest.raises(Mutation.SequenceMismatch):
            mut._assert_position_consistency("MDALYK")

    def test_single_point_deletion(self):
        mut = Mutation.from_str("D2del")
        assert mut("MDELYK")[0] == "MELYK"

    def test_deletion_with_range(self):
        mut = Mutation.from_str("D2_E3del")
        assert mut("MDELYK")[0] == "MLYK"

    def test_delins_requires_range(self):
        with pytest.raises(ValueError):
            # Insertion mutations must specify a range (with stop_char/idx)
            Mutation.from_str("D2delinsVSKG")

    def test_delins_with_range(self):
        mut = Mutation.from_str("D2_E3delinsVSKG")
        assert mut("MDELYK")[0] == "MVSKGLYK"
        assert mut("MDELYK")[1] == 2
        mut = Mutation.from_str("D2_D2delinsVSKG")
        assert mut("MDELYK")[0] == "MVSKGELYK"
        assert mut("MDELYK")[1] == 3

    def test_substitution(self):
        mut = Mutation.from_str("S4T")
        assert mut("MDESLK")[0] == "MDETLK"

    def test_insertion(self):
        mut = Mutation.from_str("K3_L4insRSG")
        assert mut("MDKLRG")[0] == "MDKRSGLRG"
        with pytest.raises(ValueError):
            Mutation.from_str("K23insRSG")

    def test_extension(self):
        mut = Mutation.from_str("*5TextAKGT")
        assert mut("MDKL")[0] == "MDKLTAKGT"

    def test_mutation_set(self):
        mruby = (
            "MNSLIKENMRMKVVLEGSVNGHQFKCTGEGEGNPYMGTQTMRIKVIEGGPLPFAFDILATSFMYGSR"
            "TFIKYPKGIPDFFKQSFPEGFTWERVTRYEDGGVITVMQDTSLEDGCLVYHAQVRGVNFPSNGAVMQ"
            "KKTKGWEPNTEMMYPADGGLRGYTHMALKVDGGGHLSCSFVTTYRSKKTVGNIKMPGIHAVDHRLER"
            "LEESDNEMFVVQREHAVAKFAGLGGG"
        )
        mruby2 = (
            "MVSKGEELIKENMRMKVVMEGSVNGHQFKCTGEGEGNPYMGTQTMRIKVIEGGPLPFAFDILATSFM"
            "YGSRTFIKYPKGIPDFFKQSFPEGFTWERVTRYEDGGVVTVMQDTSLEDGCLVYHVQVRGVNFPSNG"
            "PVMQKKTKGWEPNTEMMYPADGGLRGYTHMALKVDGGGHLSCSFVTTYRSKKTVGNIKMPGIHAVDH"
            "RLERLEESDNEMFVVQREHAVAKFAGLGGGMDELYK"
        )
        result = mutate_sequence(mruby, "N2_S3delinsVSKGEE/L15M/I102V/A119V/A131P/*228MextDELYK")
        assert result == mruby2
