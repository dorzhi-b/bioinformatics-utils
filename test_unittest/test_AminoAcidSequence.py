import unittest
from bioinformatics_utils import AminoAcidSequence


class TestBiologicalSequences(unittest.TestCase):

    def test_amino_acid_sequence_len(self):
        aa = AminoAcidSequence("ACDE")
        self.assertEqual(len(aa), 4)

    def test_amino_acid_sequence_str(self):
        aa = AminoAcidSequence("ACDE")
        self.assertEqual(str(aa), "ACDE")

    def test_amino_acid_sequence_repr(self):
        aa = AminoAcidSequence("ACDE")
        self.assertEqual(repr(aa), "AminoAcidSequence('ACDE')")

    def test_amino_acid_sequence_is_valid_alphabet(self):
        aa = AminoAcidSequence("ACDE")
        self.assertTrue(aa.is_valid_alphabet())
        invalid_aa = AminoAcidSequence("ACDEX")
        self.assertFalse(invalid_aa.is_valid_alphabet())

    def test_amino_acid_sequence_count_hydrophobic_residues(self):
        aa = AminoAcidSequence("AVLI")
        self.assertEqual(aa.count_hydrophobic_residues(), 3)

    def test_amino_acid_sequence_count_hydrophilic_residues(self):
        aa = AminoAcidSequence("DEHKR")
        self.assertEqual(aa.count_hydrophobic_residues(), 0)

if __name__ == "__main__":
    unittest.main()
