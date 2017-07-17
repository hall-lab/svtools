from unittest import TestCase, main
import tempfile
import os
import sys
import difflib
import svtools.lmerge as lmerge

class LmergeUnitTest(TestCase):
    def test_null_format_string(self):
        self.assertEqual(lmerge.null_format_string('GT:GQ:AD'), './.:.:.')
        self.assertEqual(lmerge.null_format_string('GQ:AD'), '.:.')

class LmergeIntegrationTest(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'lmerge')
        input_file = os.path.join(self.test_data_dir, 'input')
        expected_result = os.path.join(self.test_data_dir, 'lmerge_output.new.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            # FIXME this is pretty hacky
            temp_handle, sys.stdout = sys.stdout, output_handle

            lmerge.l_cluster_by_line(input_file, 0.0, 20, True)
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


    def run_integration_test_with_genotypes(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'lmerge')
        input_file = os.path.join(self.test_data_dir, 'input')
        expected_result = os.path.join(self.test_data_dir, 'lmerge_output.gt.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            # FIXME this is pretty hacky
            temp_handle, sys.stdout = sys.stdout, output_handle

            lmerge.l_cluster_by_line(input_file, 0.0, 20, True, True)
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
