import time
from unittest import TestCase, main
from svtools.vcf.file import Vcf

class TestFormat(TestCase):
    def test_init(self):
        f = Vcf.Format('GT', 1, 'String', '"Genotype"')
        self.assertEqual(f.id, 'GT')
        self.assertEqual(f.number, '1')
        self.assertEqual(f.type, 'String')
        self.assertEqual(f.desc, 'Genotype')
        self.assertEqual(f.hstring, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    def test_eq(self):
        f = Vcf.Format('GT', 1, 'String', '"Genotype"')
        g = Vcf.Format('GT', 1, 'String', 'Genotype')
        self.assertEqual(f, g)

class TestInfo(TestCase):
    def test_init(self):
        i = Vcf.Info('NS', 1, 'Integer', '"Number of Samples With Data"')
        self.assertEqual(i.id, 'NS')
        self.assertEqual(i.number, '1')
        self.assertEqual(i.type, 'Integer')
        self.assertEqual(i.desc, 'Number of Samples With Data')
        self.assertEqual(i.hstring, '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
    def test_eq(self):
        i = Vcf.Info('NS', 1, 'Integer', '"Number of Samples With Data"')
        j = Vcf.Info('NS', 1, 'Integer', 'Number of Samples With Data')
        self.assertEqual(i, j)

class TestAlt(TestCase):
    def test_init(self):
        a = Vcf.Alt('DEL:ME:ALU', '"Deletion of ALU element"')
        self.assertEqual(a.id, 'DEL:ME:ALU')
        self.assertEqual(a.desc, 'Deletion of ALU element')
        self.assertEqual(a.hstring, '##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">')
    def test_eq(self):
        a = Vcf.Alt('DEL:ME:ALU', '"Deletion of ALU element"')
        b = Vcf.Alt('DEL:ME:ALU', 'Deletion of ALU element')
        self.assertEqual(a, b)

class TestVcf(TestCase):
    def test_init(self):
        f = Vcf.Format('GT', 1, 'String', 'Genotype')
        vcf = Vcf()
        self.assertEqual(vcf.file_format, 'VCFv4.2')
        self.assertEqual(vcf.format_list, [f])

    def test_add_header(self):
        header_lines = [
                '##fileformat=VCFv4.2',
                '##fileDate=20090805',
                '##source=myImputationProgramV3.1',
                '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta',
                '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>',
                '##phasing=partial',
                '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
                '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">',
                '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">',
                '##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">',
                '##FILTER=<ID=q10,Description="Quality below 10">',
                '##FILTER=<ID=s50,Description="Less than 50% of samples have data">',
                '##ALT=<ID=DEL,Description="DELETION">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">',
                '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003']

        v = Vcf()
        v.add_header(header_lines)
        expected_header_lines = header_lines[:2] + header_lines[3:4] + header_lines[6:12] + header_lines[14:]
        expected_header_lines[1] = '##fileDate=' + time.strftime('%Y%m%d')
        self.assertEqual(v.get_header(), '\n'.join(expected_header_lines))


if __name__ == "__main__":
    main()

