from unittest import TestCase, main
import os
import time
import sys
import tempfile
import difflib
import svtools.bedpetobed12

class IntegrationTest_bedpetobed12(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'bedpetobed12')
        input_file = os.path.join(test_data_dir, 'input.bed')
        expected_result = os.path.join(test_data_dir, 'expected.bed')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.bed')
        with open(input_file) as input_stream, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.bedpetobed12.processBEDPE(input_stream, 'BEDPE', 1000000, output_handle)
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
