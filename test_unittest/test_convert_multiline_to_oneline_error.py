import unittest
from bio_files_processor import convert_multiline_fasta_to_oneline

class TestSequenceFunctions(unittest.TestCase):

    def test_convert_multiline_fasta_to_oneline_error(self):
        
        input_path = "non_existing.fasta"
        with self.assertRaises(FileNotFoundError):
            convert_multiline_fasta_to_oneline(input_path)

if __name__ == '__main__':
    unittest.main()