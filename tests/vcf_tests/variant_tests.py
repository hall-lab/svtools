from unittest import TestCase, main
from svtools.vcf.file import Vcf
from svtools.vcf.variant import Variant
from svtools.vcf.genotype import Genotype

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

    def test_parse_genotypes(self):
        genotype_field_strings = ['0/1:20', '0/0:15']
        parsed_dict = self.variant._parse_genotypes(genotype_field_strings)

        na12878_gt = Genotype(self.variant, genotype_field_strings[0].split(':'))
        na0001_gt = Genotype(self.variant, genotype_field_strings[1].split(':'))
        expected_genotype_dict = { 'NA12878': na12878_gt, 'NA0001': na0001_gt }

        self.assertEqual(parsed_dict, expected_genotype_dict)

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

    def test_genotypes(self):
        self.assertEqual([ x.get_gt_string() for x in self.variant.genotypes() ], ['0/0:9', '1/1:15'])

    def test_var_string(self):
        self.assertEqual(self.variant.get_var_string(), self.variant_line)
        self.variant.genotype('NA12878').set_format('GT', './.')
        self.assertEqual(self.variant.get_var_string(use_cached_gt_string=True), self.variant_line)
        self.assertNotEqual(self.variant.get_var_string(), self.variant_line)

if __name__ == "__main__":
    main()

