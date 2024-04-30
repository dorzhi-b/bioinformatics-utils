import os
import unittest
from Bio import SeqIO
from Bio.SeqUtils import GC
from bioinformatics_utils import filter_fastq 

class TestFilterFastq(unittest.TestCase):

    def setUp(self):
        
        self.test_input_path = "test_input.fastq"
        with open(self.test_input_path, "w") as f:
            f.write(
                "@seq1\nACTG\n+\nIIII\n"
                "@seq2\nCGCG\n+\nJJJJ\n"
                "@seq3\nTTAACC\n+\nHHHGGG\n"
            )
        
    def tearDown(self):
        
        if os.path.exists(self.test_input_path):
            os.remove(self.test_input_path)

    def test_gc_bounds(self):
        
        result = filter_fastq(self.test_input_path, gc_bounds=(0, 40))
        self.assertEqual(len(result), 2)  
        self.assertIn("seq1", result)
        self.assertIn("seq3", result)

    def test_length_bounds(self):
        
        result = filter_fastq(self.test_input_path, length_bounds=(4, 4))
        self.assertEqual(len(result), 1)  
        self.assertIn("seq1", result)

    def test_quality_threshold(self):

        result = filter_fastq(self.test_input_path, quality_threshold=35)
        self.assertEqual(len(result), 1)
        self.assertIn("seq2", result)

    def test_combined_filtering(self):

        result = filter_fastq(self.test_input_path, gc_bounds=(30, 50), length_bounds=(3, 4), quality_threshold=30)
        self.assertEqual(len(result), 1)
        self.assertIn("seq1", result)

if __name__ == "__main__":
    unittest.main()
