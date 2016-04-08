from unittest import TestCase, main
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant

class TestVariant(TestCase):
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
                '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878	NA0001' ]
        self.vcf = Vcf()
        self.vcf.add_header(header_lines)
        self.variant_line = '1	820915	5838_1	N	]GL000232.1:20940]N	0.00	.	SVTYPE=BND;STRANDS=-+:9;IMAFLAG	GT:SU	0/0:9	1/1:15'
        self.variant = Variant(self.variant_line.split('\t'), self.vcf)

    def test_set_info(self):
        self.variant.set_info('SVTYPE', 'INV')
        self.assertEqual(self.variant.info['SVTYPE'], 'INV')
        self.variant.set_info('IMAFLAG', False)
        self.assertEqual(self.variant.info['IMAFLAG'], False)
        with self.assertRaises(SystemExit) as cm:
            self.variant.set_info('SUPER', True)

    def test_get_info(self):
        self.assertEqual(self.variant.get_info('IMAFLAG'), True)
        self.assertEqual(self.variant.get_info('SVTYPE'), 'BND')
        with self.assertRaises(KeyError) as cm:
            self.variant.get_info('CALI')

    def test_get_info_string(self):
        self.assertEqual(self.variant.get_info_string(), 'SVTYPE=BND;STRANDS=-+:9;IMAFLAG')
        self.variant.set_info('IMAFLAG', False)
        self.assertEqual(self.variant.get_info_string(), 'SVTYPE=BND;STRANDS=-+:9')

    def test_get_format_string(self):
        self.assertEqual(self.variant.get_format_string(), 'GT:SU') 

    def test_get_gt_string(self):
        self.assertEqual(self.variant.get_gt_string(), '0/0:9	1/1:15')

    def test_genotype(self):
        self.assertEqual(self.variant.genotype('NA12878').get_gt_string(), '0/0:9')

    def test_var_string(self):
        self.assertEqual(self.variant.get_var_string(), self.variant_line)

if __name__ == "__main__":
    main()

