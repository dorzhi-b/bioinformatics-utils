import unittest
from bioinformatics_utils import DNASequence


class TestBiologicalSequences(unittest.TestCase):

    def test_dna_sequence_len(self):
        dna = DNASequence("ATCG")
        self.assertEqual(len(dna), 4)

    def test_dna_sequence_getitem(self):
        dna = DNASequence("ATCG")
        self.assertEqual(dna[0], "A")
        self.assertEqual(dna[1], "T")

    def test_dna_sequence_str(self):
        dna = DNASequence("ATCG")
        self.assertEqual(str(dna), "ATCG")

    def test_dna_sequence_repr(self):
        dna = DNASequence("ATCG")
        self.assertEqual(repr(dna), "DNASequence('ATCG')")

    def test_dna_sequence_complement(self):
        dna = DNASequence("ATCG")
        self.assertEqual(dna.complement(), "TAGC")

    def test_dna_sequence_gc_content(self):
        dna = DNASequence("ATCG")
        self.assertAlmostEqual(dna.gc_content(), 0.5) 

    def test_dna_sequence_transcribe(self):
        dna = DNASequence("ATCG")
        rna = dna.transcribe()
        self.assertEqual(str(rna), "AUCG")

    def test_dna_sequence_is_valid_alphabet(self):
        dna = DNASequence("ATCG")
        self.assertTrue(dna.is_valid_alphabet())
        invalid_dna = DNASequence("ATCZ")
        self.assertFalse(invalid_dna.is_valid_alphabet())

if __name__ == "__main__":
    unittest.main()
