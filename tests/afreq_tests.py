from unittest import TestCase, main
import os
import time
import svtools.afreq
import sys
import tempfile
import difflib

class IntegrationTest_afreq(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'afreq')
        input = os.path.join(test_data_dir, 'input.vcf')
        expected_result = os.path.join(test_data_dir, 'expected.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with open(input, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            updater = svtools.afreq.UpdateInfo(input_handle)
            updater.execute(output_handle)
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

class AfreqUiTest(TestCase):
    def test_parser(self):
        parser = svtools.afreq.command_parser()
        args1 = parser.parse_args([])
        self.assertIsNone(args1.input_vcf)
        args2 = parser.parse_args(['somefile'])
        self.assertEqual(args2.input_vcf, 'somefile')

class TestUpdateInfo(TestCase):
    def test_numeric_alleles(self):
        instance = svtools.afreq.UpdateInfo(None)
        self.assertEqual(instance.numeric_alleles('0/0'), [0, 0])
        self.assertEqual(instance.numeric_alleles('0|1'), [0, 1])
    

if __name__ == "__main__":
    main()

