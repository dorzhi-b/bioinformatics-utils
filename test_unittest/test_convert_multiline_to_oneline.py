import os
import unittest
from bio_files_processor import convert_multiline_fasta_to_oneline

class TestSequenceFunctions(unittest.TestCase):

    def test_convert_multiline_fasta_to_oneline(self):

        test_input_path = "test_input.fasta"
        test_output_path = "test_output.fasta"
        with open(test_input_path, "w") as f:
            f.write(">id1\nATCG\n>id2\nGCTA")

        convert_multiline_fasta_to_oneline(test_input_path, test_output_path)
        
        self.assertTrue(os.path.exists(test_output_path))
        with open(test_output_path, "r") as f:
            content = f.read()
            self.assertEqual(content, ">id1\nATCG>id2\nGCTA")

if __name__ == '__main__':
    unittest.main()