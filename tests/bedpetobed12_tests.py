from unittest import TestCase, main
import os
import sys
import tempfile
import difflib
import svtools.bedpetobed12

class TestBedpeToBlockedBedConverter(TestCase):

    def setUp(self):
        self.converter = svtools.bedpetobed12.BedpetoBlockedBedConverter('TESTCONV', 20)

    def test_track_name(self):
        self.assertEqual(self.converter.track_name(), 'track name=TESTCONV itemRgb=On\n')

    def test_get_color(self):
        self.assertEqual(self.converter.get_color('DEL', 10), '153,0,0')
        self.assertEqual(self.converter.get_color('DEL', 50), self.converter.distant_color)
        self.assertEqual(self.converter.get_color('ITX', 10), self.converter.unknown_close_color)

    def test_bed12_name(self):
        self.assertEqual(self.converter.bed12_name('ITX', '22', None), 'ITX;ID=22')
        self.assertEqual(self.converter.bed12_name('ITX', '22', '0.25'), 'ITX;ID=22;AF=0.25')
        self.assertEqual(self.converter.bed12_name('ITX', '22', '0.25', ('+', '-')), 'ITX;ID=22;AF=0.25;STR=+-')

    def test_distant_coordinates(self):
        self.assertEqual(self.converter.distant_coordinates('+', 600, 600), (600, 1100))
        self.assertEqual(self.converter.distant_coordinates('-', 600, 600), (100, 600))

    def test_distant_block_sizes(self):
        self.assertEqual(self.converter.distant_block_sizes('+', 600, 800), (200, 1))
        self.assertEqual(self.converter.distant_block_sizes('-', 600, 800), (1, 200))

    def test_distant_block_starts(self):
        self.assertEqual(self.converter.distant_block_starts('+', 600, 800), (0, 699))
        self.assertEqual(self.converter.distant_block_starts('-', 600, 800), (0, 500))

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
