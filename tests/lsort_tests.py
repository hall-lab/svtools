from unittest import TestCase, main
from tempfile import gettempdir
import svtools.lsort as lsort

class Test_lsort(TestCase):
    def test_parser(self):
        parser = lsort.command_parser()

        args = parser.parse_args('file1 file2 file3'.split())
        self.assertEqual(args.vcf_files, ['file1', 'file2', 'file3'])
        self.assertEqual(args.batchsize, 200)
        self.assertEqual(args.tempdir, gettempdir())

        args2 = parser.parse_args('-b 2 -t temp file1 file2'.split())
        self.assertEqual(args2.batchsize, 2)
        self.assertEqual(args2.tempdir, 'temp')
        self.assertEqual(args2.vcf_files, ['file1', 'file2'])

if __name__ == "__main__":
    main()
