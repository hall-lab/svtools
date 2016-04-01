from unittest import TestCase, main
import os
import svtools.vcfpaste
import sys
import tempfile
import difflib

class IntegrationTest_vcfpaste(TestCase):
    # FIXME We really don't need to have this stuff run with every test. Run once...
    def setUp(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'vcfpaste')
        # glob vcfs
        vcfs = list()
        for sample in ('NA12878', 'NA12891', 'NA12892'):
            vcfs.append(os.path.join(self.test_data_dir, '{0}.vcf'.format(sample)))
        # write out list since we have the paths and have to get those right
        temp_descriptor, self.list_of_vcfs = tempfile.mkstemp()
        temp_handle = os.fdopen(temp_descriptor, 'w') 

        temp_descriptor2, self.list_of_vcfs_with_truncated = tempfile.mkstemp()
        temp_handle2 = os.fdopen(temp_descriptor2, 'w')

        temp_descriptor3, self.list_of_gz_vcfs = tempfile.mkstemp()
        temp_handle3 = os.fdopen(temp_descriptor3, 'w')

        temp_descriptor4, self.list_of_vcfs_with_8col = tempfile.mkstemp()
        temp_handle4 = os.fdopen(temp_descriptor4, 'w')
        
        for vcf_path in vcfs:
            temp_handle.write(vcf_path + '\n')
            temp_handle2.write(vcf_path + '\n')
            temp_handle4.write(vcf_path + '\n')
            temp_handle3.write(vcf_path + '.gz\n')
        temp_handle.close()
        temp_handle3.close()

        truncated_vcf = os.path.join(self.test_data_dir, 'truncated.vcf')
        temp_handle2.write(truncated_vcf + '\n')
        temp_handle2.close()

        too_skinny_vcf = os.path.join(self.test_data_dir, '8col_NA12878.vcf')
        temp_handle4.write(too_skinny_vcf + '\n')
        temp_handle4.close()

        self.master = os.path.join(self.test_data_dir, 'master.vcf')
        self.thin_master = os.path.join(self.test_data_dir, 'thin_master.vcf')
        self.too_thin_master = os.path.join(self.test_data_dir, '6col_master.vcf')

    def tearDown(self):
        os.remove(self.list_of_vcfs)
        os.remove(self.list_of_vcfs_with_truncated)
        os.remove(self.list_of_gz_vcfs)

    def run_integration_test_without_master(self):
        expected_result = os.path.join(self.test_data_dir, 'expected_no_master.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        try:
            paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs, master=None, sum_quals=True)
            paster.execute(output_handle)
        finally:
            output_handle.close()
        expected_lines = open(expected_result).readlines()
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = '\n'.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)

    def run_integration_test_without_master_gzipped(self):
        expected_result = os.path.join(self.test_data_dir, 'expected_no_master.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        try:
            paster = svtools.vcfpaste.Vcfpaste(self.list_of_gz_vcfs, master=None, sum_quals=True)
            paster.execute(output_handle)
        finally:
            output_handle.close()
        expected_lines = open(expected_result).readlines()
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = '\n'.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)

    def run_integration_test_with_master(self):
        master_file = os.path.join(self.test_data_dir, 'master.vcf')
        expected_result = os.path.join(self.test_data_dir, 'expected_master.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        try:
            paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs, master=master_file, sum_quals=True)
            paster.execute(output_handle)
        finally:
            output_handle.close()
        expected_lines = open(expected_result).readlines()
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = '\n'.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)

    def run_integration_test_with_thin_master(self):
        expected_result = os.path.join(self.test_data_dir, 'expected_master.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        try:
            paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs, master=self.thin_master, sum_quals=True)
            paster.execute(output_handle)
        finally:
            output_handle.close()
        expected_lines = open(expected_result).readlines()
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = '\n'.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)

    def run_integration_test_with_too_thin_master(self):
        expected_result = os.path.join(self.test_data_dir, 'expected_master.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs, master=self.too_thin_master, sum_quals=True)
        with self.assertRaises(SystemExit) as cm:
            paster.execute(output_handle)
            exception = cm.exception
            self.assertEqual(exception.error_code, 1)
            output_handle.close()
        os.remove(temp_output_path)

    def run_integration_test_with_truncated_vcf(self):
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs_with_truncated, master=None, sum_quals=True)
        with self.assertRaises(SystemExit) as cm:
            paster.execute(output_handle)
            exception = cm.exception
            self.assertEqual(exception.error_code, 1)
            output_handle.close()
        os.remove(temp_output_path)

    def run_integration_test_with_too_thin_vcf(self):
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        output_handle = os.fdopen(temp_descriptor, 'w')
        paster = svtools.vcfpaste.Vcfpaste(self.list_of_vcfs_with_8col, master=None, sum_quals=True)
        with self.assertRaises(SystemExit) as cm:
            paster.execute(output_handle)
            exception = cm.exception
            self.assertEqual(exception.error_code, 1)
            output_handle.close()
        os.remove(temp_output_path)

class Test_vcfpaste(TestCase):

    def test_init_w_defaults(self):
        paster = svtools.vcfpaste.Vcfpaste('a_file_o_vcf_filenames')
        self.assertEqual(paster.vcf_list, 'a_file_o_vcf_filenames')
        self.assertIsNone(paster.master)
        self.assertIsNone(paster.sum_quals)

    def test_init_w_specified(self):
        paster = svtools.vcfpaste.Vcfpaste('some_file', 'master_blaster', True)
        self.assertEqual(paster.vcf_list, 'some_file')
        self.assertEqual(paster.master, 'master_blaster')
        self.assertTrue(paster.sum_quals)

class Test_vcfpaste_ui(TestCase):

    def test_add_arguments(self):
        parser = svtools.vcfpaste.command_parser()
        args1 = parser.parse_args(['--vcf-list', 'some_list'])
        self.assertEqual(args1.vcf_list, 'some_list')
        self.assertFalse(args1.sum_quals)
        self.assertIsNone(args1.master)

        args2 = parser.parse_args(['-f', 'some_list', '-m', 'some_master', '-q'])
        self.assertEqual(args2.vcf_list, 'some_list')
        self.assertTrue(args2.sum_quals)
        self.assertEqual(args2.master, 'some_master')
        

if __name__ == "__main__":
    main()

