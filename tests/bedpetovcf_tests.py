from unittest import TestCase, main
import os
import time
import sys
import tempfile
import difflib
import svtools.bedpetovcf

class IntegrationTest_bedpetovcf(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'bedpetovcf')
        input = os.path.join(test_data_dir, 'input.bed')
        expected_result = os.path.join(test_data_dir, 'expected.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with open(input) as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.bedpetovcf.bedpeToVcf(input_handle, output_handle)
        expected_lines = open(expected_result).readlines()
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = ''.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        #os.remove(temp_output_path)

if __name__ == "__main__":
    main()
