from unittest import TestCase
from svtools.vcftobedpeconverter import VcfToBedpeConverter
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from svtools.bedpe import Bedpe

class TestVcfToBedpeConverter(TestCase):
    def setUp(self):
        self.converter = VcfToBedpeConverter()
        header_lines = [
                '##fileformat=VCFv4.2',
                '##fileDate=20090805',
                '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta',
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001']
        self.vcf = Vcf()
        self.vcf.add_header(header_lines)

    def test_bnd_alt_string(self):
        self.assertEqual(self.converter.parse_bnd_alt_string('A[1:6['), ('[', '1', 6))
        self.assertEqual(self.converter.parse_bnd_alt_string('A]1:6]'), (']', '1', 6))
        self.assertEqual(self.converter.parse_bnd_alt_string(']1:6]A'), (']', '1', 6))
        with self.assertRaises(AssertionError):
            self.converter.parse_bnd_alt_string(']1:6[A')
        with self.assertRaises(AssertionError):
            self.converter.parse_bnd_alt_string('1')

    def test_bnd_breakpoints(self):
        vcf_array1 = ['1', '20000', '235', 'T', 'A[1:6[', '0.00', '.', '.', 'GT', '0/0']
        v1 = Variant(vcf_array1, self.vcf)
        self.assertEqual(
                self.converter.bnd_breakpoints(v1),
                ('1', 20000, 20000, '1', 5, 5, '+', '-'))
        vcf_array2 = ['1', '20000', '235', 'T', ']1:6]N', '0.00', '.', '.', 'GT', '0/0']
        v2 = Variant(vcf_array2, self.vcf)
        self.assertEqual(
                self.converter.bnd_breakpoints(v2),
                ('1', 19999, 19999, '1', 6, 6, '-', '+'))

    def test_simple_breakpoints(self):
        vcf_array1 = ['1', '20000', '235', 'T', '<DEL>', '0.00', '.', 'END=20500', 'GT', '0/0']
        v1 = Variant(vcf_array1, self.vcf)
        self.assertEqual(
                self.converter.simple_breakpoints(v1),
                ('1', 20000, 20000, '1', 20500, 20500, '+', '-'))
        vcf_array2 = ['1', '20000', '235', 'T', '<DEL>', '0.00', '.', 'END=20500;STRANDS=-+:2', 'GT', '0/0']
        v2 = Variant(vcf_array2, self.vcf)
        self.assertEqual(
                self.converter.simple_breakpoints(v2),
                ('1', 20000, 20000, '1', 20500, 20500, '-', '+'))
        vcf_array3 = ['1', '20000', '235', 'T', '<DEL>', '0.00', '.', 'STRANDS=--:2', 'GT', '0/0']
        v3 = Variant(vcf_array3, self.vcf)
        with self.assertRaises(ValueError):
            self.converter.simple_breakpoints(v3)

    def test_adjust_coordinate(self):
        vcf_array1 = ['1', '20000', '235', 'T', '<DEL>', '0.00', '.', 'CIEND=-50,50', 'GT', '0/0']
        v1 = Variant(vcf_array1, self.vcf)
        self.assertEqual(
                self.converter.adjust_coordinate(v1, 'CIEND', 500, 1000),
                (450, 1050))
        self.assertEqual(
                self.converter.adjust_coordinate(v1, 'CIPOS', 500, 1000),
                (500, 1000))
        vcf_array2 = ['1', '20000', '235', 'T', '<DEL>', '0.00', '.', 'CIEND=50', 'GT', '0/0']
        v2 = Variant(vcf_array2, self.vcf)
        with self.assertRaises(ValueError):
            self.converter.adjust_coordinate(v2, 'CIEND', 500, 1000)
