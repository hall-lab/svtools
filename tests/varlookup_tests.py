from unittest import TestCase, main
import os
import time
import sys
import tempfile
import difflib
import svtools.varlookup

class IntegrationTest_varlookup(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'varlookup')
        input_a = os.path.join(test_data_dir, 'input_a.bed')
        input_b = os.path.join(test_data_dir, 'input_b.bed')
        expected_result = os.path.join(test_data_dir, 'expected.bed')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.bed')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.varlookup.varLookup(input_a, input_b, output_handle, 50, '#', 'TEST')
        expected_lines = open(expected_result).readlines()
        # set timestamp for diff
        expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
        produced_lines = open(temp_output_path).readlines()
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = ''.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)

    def run_issue_209_regression_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'varlookup')
        input_a = os.path.join(test_data_dir, 'input_a1.bed')
        input_b = os.path.join(test_data_dir, 'input_b1.bed')
        expected_result = os.path.join(test_data_dir, 'expected1.bed')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.bed')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.varlookup.varLookup(input_a, input_b, output_handle, 50, '#', 'TEST')
        expected_lines = open(expected_result).readlines()
        # set timestamp for diff
        expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
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
