import os
import unittest
from bio_files_processor import change_fasta_start_pos

class TestSequenceFunctions(unittest.TestCase):

    def test_change_fasta_start_pos(self):

        test_input_path = "test_input.fasta"
        test_output_path = "test_output.fasta"
        with open(test_input_path, "w") as f:
            f.write(">id1\nATCG")

        change_fasta_start_pos(test_input_path, 2, test_output_path)

        self.assertTrue(os.path.exists(test_output_path))
        with open(test_output_path, "r") as f:
            content = f.read()
            self.assertEqual(content, ">id1\nCGAT")

if __name__ == '__main__':
    unittest.main()