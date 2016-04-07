from unittest import TestCase, main
import os
import sys
import tempfile
import difflib
import svtools.prune

class IntegrationTest_bedpetovcf(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'prune')
        input = os.path.join(test_data_dir, 'input.no_missing.bed')
        expected_result = os.path.join(test_data_dir, 'expected.no_missing.bed')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.bed')
        with open(input) as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            pruner = svtools.prune.Pruner(50, None)
            pruner.cluster_bedpe(input_handle, output_handle, False)
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
