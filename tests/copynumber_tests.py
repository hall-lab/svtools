from unittest import TestCase, main
import os
import time
import sys
import tempfile
import difflib
import svtools.copynumber

class IntegrationTest_copynumber(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'copynumber')
        input = os.path.join(test_data_dir, 'input.vcf')
        expected_result = os.path.join(test_data_dir, 'expected.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with open(input, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.copynumber.write_copynumber(input_handle, 'NA12878', output_handle, ['1.99', '0.13', '5.32', '2.76'])
            expected_lines = open(expected_result).readlines()
            # set timestamp for diff
            expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
            #remove reference line
            produced_lines = open(temp_output_path).readlines()
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            result = ''.join(diff)
            if result != '':
                for line in result:
                    sys.stdout.write(line)
                self.assertFalse(result)
        os.remove(temp_output_path)

#class CopynumberUiTest(TestCase):
#    def test_parser(self):
#        pass
#        parser = svtools.copynumber.command_parser()
#        args1 = parser.parse_args([])
#        self.assertIsNone(args1.input_vcf)
#        args2 = parser.parse_args(['somefile'])
#        self.assertEqual(args2.input_vcf, 'somefile')

if __name__ == "__main__":
    main()

