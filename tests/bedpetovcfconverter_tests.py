from unittest import TestCase
from svtools.bedpetovcfconverter import BedpeToVcfConverter
from svtools.vcf.file import Vcf
from svtools.bedpe import Bedpe

class TestBedpeToVcfConverter(TestCase):

    def setUp(self):
        vcf = Vcf()
        self.converter = BedpeToVcfConverter(vcf)

    def test_adjust_by_tag(self):
        entry1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'SVTYPE=BND;AF=0.2', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(entry1)
        self.assertEqual(self.converter.adjust_by_tag(b1, 'CIPOS', '+', 200), 200)
        self.assertEqual(self.converter.adjust_by_tag(b1, 'CIPOS', '-', 200), 201)
        b1.info1 = 'SVTYPE=BND;AF=0.2;CIPOS=-2,3'
        self.assertEqual(self.converter.adjust_by_tag(b1, 'CIPOS', '-', 200), 203)
        self.assertEqual(self.converter.adjust_by_tag(b1, 'CIPOS', '+', 200), 202)

    def test_determine_sep(self):
        self.assertEqual(self.converter.determine_sep('-'), '[')
        self.assertEqual(self.converter.determine_sep('+'), ']')

    def test_determine_flanks(self):
        self.assertEqual(self.converter.determine_flanks('-'), ('', 'N'))
        self.assertEqual(self.converter.determine_flanks('+'), ('N', ''))

    def test_bnd_alt_string(self):
        self.assertEqual(self.converter.bnd_alt_string('+', '-', '2', '22222'), 'N[2:22222[')
        self.assertEqual(self.converter.bnd_alt_string('-', '+', '2', '22222'), ']2:22222]N')
