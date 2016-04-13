from unittest import TestCase, main
import os
import sys
import tempfile
import difflib
import re
import svtools.vcftobedpe
import svtools.bedpetovcf
import svtools.vcfsort

class FakeArgs(object):
    def __init__(self, input, output):
        self.input = input
        self.output = output

class RoundtripTest(TestCase):
    def test_roundtrip_file_conversion(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'roundtrip')
        input_filepath = os.path.join(test_data_dir, 'input.vcf')
        expected_result = input_filepath
        temp_bedpe_descriptor, temp_bedpe_path = tempfile.mkstemp(suffix='.bedpe')
        temp_vcf_descriptor, temp_vcf_path = tempfile.mkstemp(suffix='.vcf')
        temp_output_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with open(input_filepath) as input_handle, os.fdopen(temp_bedpe_descriptor, 'w') as output_handle:
            svtools.vcftobedpe.vcfToBedpe(input_handle, output_handle)

        #do conversion back to vcf
        with open(temp_bedpe_path) as input_handle, os.fdopen(temp_vcf_descriptor, 'w') as output_handle:
            svtools.bedpetovcf.bedpeToVcf(input_handle, output_handle)
        #sort output file
        svtools.vcfsort.run_from_args(FakeArgs(temp_vcf_path, temp_output_path))

        expected_lines = [re.sub('#fileDate=.*', '#fileDate=<DUMMY_DATE>', line) for line in open(expected_result)]
        produced_lines = [re.sub('#fileDate=.*', '#fileDate=<DUMMY_DATE>', line) for line in open(temp_output_path)]
        diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
        result = ''.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertFalse(result)
        os.remove(temp_output_path)
        os.remove(temp_vcf_path)
        os.remove(temp_bedpe_path)

if __name__ == "__main__":
    main()
