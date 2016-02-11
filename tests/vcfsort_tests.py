from unittest import TestCase, main
import os
import sys
import tempfile
import difflib
import svtools.vcfsort

class FakeArgs(object):
    def __init__(self, input, output):
        self.input = input
        self.output = output


class IntegrationTest_vcfsort(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'vcfsort')
        input = os.path.join(test_data_dir, 'input.vcf')
        expected_result = os.path.join(test_data_dir, 'expected.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        os.close(temp_descriptor)
        svtools.vcfsort.run_from_args(FakeArgs(input, temp_output_path))

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


