from unittest import TestCase, main
from svtools.bedpe import Bedpe

class BedpeTests(TestCase):
    def test_parse_score(self):
        self.assertEqual(Bedpe.parse_score('20'), 20)
        self.assertEqual(Bedpe.parse_score('.'), '.')

    def test_malformed(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.malformedFlag, 1)
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'MISSING' ]
        b2 = Bedpe(entry2)
        self.assertEqual(b2.malformedFlag, 2)
        self.assertEqual(b2.info1, entry2[18]) 

    def test_info(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.info, 'SVTYPE=BND;AF=0.2')
        
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'MISSING' ]
        b2 = Bedpe(entry2)
        self.assertEqual(b2.info, 'SVTYPE=BND;AF=0.2')

        entry3 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'SECONDARY' ]
        b3 = Bedpe(entry3)
        self.assertEqual(b3.info, 'SVTYPE=BND;AF=0.2')

    def test_set_info(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND' ]
        b1 = Bedpe(entry1)
        b1.set_info('AF', '0.2')
        self.assertEqual(b1.info, 'SVTYPE=BND;AF=0.2')
        
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND', 'MISSING' ]
        b2 = Bedpe(entry2)
        b2.set_info('AF', '0.2')
        self.assertEqual(b2.info, 'SVTYPE=BND;AF=0.2')

        entry3 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND', 'SECONDARY' ]
        b3 = Bedpe(entry3)
        b3.set_info('AF', '0.2')
        self.assertEqual(b3.info1, 'SVTYPE=BND;AF=0.2')
        self.assertEqual(b3.info2, 'SECONDARY;AF=0.2')

        entry4 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND', '.' ]
        b4 = Bedpe(entry4)
        b4.set_info('PRESENT', None)
        self.assertEqual(b4.info, 'SVTYPE=BND;PRESENT')
        self.assertEqual(b4.info2, '.')


    def test_retrieve_svtype(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.retrieve_svtype(), 'BND')
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'AF=0.2', 'AF=0.2' ]
        with self.assertRaises(SystemExit):
            b = Bedpe(entry2)

    def test_retrieve_af(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.retrieve_af(), '0.2')
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND', 'SVTYPE=BND' ]
        b2 = Bedpe(entry2)
        self.assertIsNone(b2.retrieve_af())

    def test_str(self):
        # Note that we are testing float to float equivalence. Actually passing in an integer will result in it being converted to float with
        # with decimal place
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57.0', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(str(b1), '\t'.join(entry1))


if __name__ == "__main__":
    main()

