import unittest
from bioinformatics_utils import RNASequence


class TestBiologicalSequences(unittest.TestCase):

    def test_rna_sequence_len(self):
        rna = RNASequence("AUCG")
        self.assertEqual(len(rna), 4)

    def test_rna_sequence_str(self):
        rna = RNASequence("AUCG")
        self.assertEqual(str(rna), "AUCG")

    def test_rna_sequence_repr(self):
        rna = RNASequence("AUCG")
        self.assertEqual(repr(rna), "RNASequence('AUCG')")

    def test_rna_sequence_is_valid_alphabet(self):
        rna = RNASequence("AUCG")
        self.assertTrue(rna.is_valid_alphabet())
        invalid_rna = RNASequence("AUCGX")
        self.assertFalse(invalid_rna.is_valid_alphabet())
        
if __name__ == "__main__":
    unittest.main()
