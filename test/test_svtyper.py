import unittest
from svtyper import *

class TestCigarParsing(unittest.TestCase):
    def test_cigarstring_to_tuple(self):
        string1 = '5H3S2D1N5M3I2P2X1='
        self.assertEqual(cigarstring_to_tuple(string1),
                [(5, 5), (4, 3), (2, 2), (3, 1),
                    (0, 5), (1, 3), (6, 2), (8, 2),
                    (7, 1)])

    def test_get_query_pos_from_cigar(self):
        # forward
        cigar_string = '2S3M1D2M2I3M3S'
        cigar = cigarstring_to_tuple(cigar_string)
        query_pos = get_query_pos_from_cigar(cigar, True)
        self.assertEqual(query_pos.query_start, 3)
        self.assertEqual(query_pos.query_end, 13)
        self.assertEqual(query_pos.query_length, 15)

        # get_query_pos_from_cigar currently modifies the cigar list in place.
        # that's why the code below doesn't work as intended.
        query_pos = get_query_pos_from_cigar(cigar, False)
        self.assertEqual(query_pos.query_start, 2)
        self.assertEqual(query_pos.query_end, 12)
        self.assertEqual(query_pos.query_length, 15)

    def test_get_reference_end_from_cigar(self):
        cigar_string = '2S5M3D2M3S'
        cigar = cigarstring_to_tuple(cigar_string)
        self.assertEqual(get_reference_end_from_cigar(1, cigar), 11)

    def test_get_start_diagonal(self):
        cigar_string = '2S5M3D1I1M3S'
        split_piece = SplitRead.SplitPiece(1, 25, True, cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(get_start_diagonal(split_piece), 23)
        split_piece2 = SplitRead.SplitPiece(1, 25, False, cigarstring_to_tuple(cigar_string), 60)
        self.assertEqual(get_start_diagonal(split_piece2), 23)

    def test_get_end_diagonal(self):
        cigar_string = '2S5M3D2I1M3S'
        split_piece = SplitRead.SplitPiece(1, 25, True, cigarstring_to_tuple(cigar_string), 60)
        split_piece.set_reference_end(34)
        self.assertEqual(get_end_diagonal(split_piece), 34 - (2 + 8))
        split_piece2 = SplitRead.SplitPiece(1, 25, False, cigarstring_to_tuple(cigar_string), 60)
        split_piece2.set_reference_end(34)
        self.assertEqual(get_end_diagonal(split_piece2), 34 - (2 + 8))

if __name__ == '__main__':
    suite =  unittest.TestLoader().loadTestsFromTestCase(TestCigarParsing)
    unittest.TextTestRunner(verbosity=2).run(suite)
