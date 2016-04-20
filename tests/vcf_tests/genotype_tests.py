from unittest import TestCase, main
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from svtools.vcf.genotype import Genotype

class TestGenotype(TestCase):
    def setUp(self):
        header_lines = [
                '##fileformat=VCFv4.2',
                '##fileDate=20151202',
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
                '##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">',
                '##INFO=<ID=IMAFLAG,Number=.,Type=Flag,Description="Test Flag code">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">',
                '##FORMAT=<ID=INACTIVE,Number=1,Type=Integer,Description="A format not in use">',
                '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878' ]
        self.vcf = Vcf()
        self.vcf.add_header(header_lines)
        self.variant_line = '1	820915	5838_1	N	]GL000232.1:20940]N	0.00	.	SVTYPE=BND;STRANDS=-+:9;IMAFLAG	GT:SU	0/0:9'
        self.variant = Variant(self.variant_line.split('\t'), self.vcf)

    def test_equal(self):
        g1 = Genotype(self.variant, ['0/1'])
        g1.set_format('INACTIVE', 10)
        g2 = Genotype(self.variant, ['0/1'])
        g2.set_format('INACTIVE', 10)
        self.assertEqual(g1, g2)
    
    def test_set_format(self):
        g = Genotype(self.variant, ['0/1'])
        self.assertFalse('INACTIVE' in self.variant.format_dict)
        g.set_format('INACTIVE', 10)
        self.assertEqual(g.get_format('INACTIVE'), 10)
        self.assertTrue('INACTIVE' in self.variant.format_dict)

    def test_get_format(self):
        g = Genotype(self.variant, ['0/1'])
        g.set_format('INACTIVE', 10)
        self.assertEqual(g.get_format('INACTIVE'), 10)

    def test_get_gt_string(self):
        g = Genotype(self.variant, ['0/1'])
        g.set_format('INACTIVE', 10)
        self.assertEqual(g.get_gt_string(), '0/1:.:10')

if __name__ == "__main__":
    main()
