from unittest import TestCase, main
import tempfile
import os
import sys
import difflib
import svtools.lsort as lsort

class Test_lsort(TestCase):
    def test_parser(self):
        parser = lsort.command_parser()

        args = parser.parse_args('file1 file2 file3'.split())
        self.assertEqual(args.vcf_files, ['file1', 'file2', 'file3'])
        self.assertEqual(args.batchsize, 200)
        self.assertEqual(args.tempdir, tempfile.gettempdir())

        args2 = parser.parse_args('-b 2 -t temp file1 file2'.split())
        self.assertEqual(args2.batchsize, 2)
        self.assertEqual(args2.tempdir, 'temp')
        self.assertEqual(args2.vcf_files, ['file1', 'file2'])
    
    def test_lsort_init_defaults(self):
        file_list = ['file1', 'file2']
        lsort_class = lsort.Lsort(file_list)
        self.assertEqual(lsort_class.vcf_file_names, file_list)
        self.assertEqual(lsort_class.batchsize, 200)
        self.assertEqual(lsort_class.tempdir, tempfile.gettempdir())

    def test_lsort_init_full(self):
        file_list = ['file1', 'file2']
        lsort_class = lsort.Lsort(file_list, tempdir='tempydir', batchsize=5 )
        self.assertEqual(lsort_class.vcf_file_names, file_list)
        self.assertEqual(lsort_class.batchsize, 5)
        self.assertEqual(lsort_class.tempdir, 'tempydir')

class LsortIntegrationTest(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'lsort')
        # glob vcfs
        vcfs = list()
        for sample in ('NA12878', 'NA12891', 'NA12892'):
            vcfs.append(os.path.join(self.test_data_dir, '{0}.vcf'.format(sample)))
        expected_result = os.path.join(self.test_data_dir, 'lsort_expected')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            sorter = lsort.Lsort(vcfs, tempdir=None, batchsize=2, output_handle=output_handle)
            sorter.execute()
            output_handle.flush()
            expected_lines = open(expected_result).readlines()
            produced_lines = open(temp_output_path).readlines()
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            result = ''.join(diff)
            if result != '':
                for line in result:
                    sys.stdout.write(line)
                self.assertFalse(result)
        os.remove(temp_output_path)

    def run_file_list_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'lsort')
        # glob vcfs
        vcfs = list()
        for sample in ('NA12878', 'NA12891', 'NA12892'):
            vcfs.append(os.path.join(self.test_data_dir, '{0}.vcf.gz'.format(sample)))
        expected_result = os.path.join(self.test_data_dir, 'lsort_expected')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            sorter = lsort.Lsort(vcfs, tempdir=None, batchsize=2, skip_ref=True, output_handle=output_handle)
            sorter.execute()
            output_handle.flush()
            expected_lines = open(expected_result).readlines()
            produced_lines = open(temp_output_path).readlines()
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            result = ''.join(diff)
            if result != '':
                for line in result:
                    sys.stdout.write(line)
                self.assertFalse(result)
        os.remove(temp_output_path)


if __name__ == "__main__":
    main()
