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

class TestAlt(TestCase):
    def test_init(self):
        a = Vcf.Alt('DEL:ME:ALU', '"Deletion of ALU element"')
        self.assertEqual(a.id, 'DEL:ME:ALU')
        self.assertEqual(a.desc, 'Deletion of ALU element')
        self.assertEqual(a.hstring, '##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">')

class TestVcf(TestCase):
    def tests_init(self):
        f = Vcf.Format('GT', 1, 'String', 'Genotype')
        vcf = Vcf()
        self.assertEqual(vcf.file_format, 'VCFv4.2')
        self.assertEqual(vcf.format_list, [f])


if __name__ == "__main__":
    main()

