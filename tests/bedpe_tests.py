from unittest import TestCase, main
from svtools.bedpe import Bedpe

class BedpeTests(TestCase):
    def test_parse_score(self):
        self.assertEqual(Bedpe.parse_score('20'), 20)
        self.assertEqual(Bedpe.parse_score('.'), '.')

    def test_malformed(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.malformedFlag, 1)
        self.assertNotEqual(b1.info1, 'MISSING') 
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'SVTYPE=BND;AF=0.2', 'MISSING' ]
        b2 = Bedpe(entry2)
        self.assertEqual(b2.malformedFlag, 2)
        self.assertEqual(b2.info1, entry2[12]) 

    def test_retrieve_svtype(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.retrieve_svtype(), 'BND')
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'AF=0.2', 'AF=0.2' ]
        with self.assertRaises(SystemExit):
            b = Bedpe(entry2)

    def test_retrieve_af(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(b1.retrieve_af(), '0.2')
        entry2 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', 'SVTYPE=BND', 'SVTYPE=BND' ]
        b2 = Bedpe(entry2)
        self.assertIsNone(b2.retrieve_af())

if __name__ == "__main__":
    main()

