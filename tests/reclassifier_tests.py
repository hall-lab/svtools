from unittest import TestCase, main
import os
import time
import sys
import tempfile
import difflib
import svtools.sv_classifier
import gzip

def decode(x):
    if isinstance(x, bytes):
        return x.decode()
    return x

class IntegrationTest_sv_classify(TestCase):

    def test_chromosome_prefix(self):
        self.assertEqual(svtools.sv_classifier.chromosome_prefix('chrBLAH'), 'BLAH')
        self.assertEqual(svtools.sv_classifier.chromosome_prefix('BLAH'), 'chrBLAH')

    def test_integration_nb(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'sv_classifier')
        input = os.path.join(test_data_dir, 'reclass.test.vcf.gz')

        expected_result = os.path.join(test_data_dir, 'output.nb.vcf.gz')
        annot=os.path.join(test_data_dir, 'repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz')
        sex_file=os.path.join(test_data_dir, 'ceph.sex.txt')
        train=os.path.join(test_data_dir, 'training.vars.vcf.gz')

        diags_handle, diags_file = tempfile.mkstemp(suffix='.txt')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        sex=open(sex_file, 'r')
        sex_chrom_names = set(('X', 'Y'))
        with gzip.open(input, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.sv_classifier.run_reclassifier(input_handle, output_handle, sex, sex_chrom_names, annot,  0.9, None, 1.0, 0.2, train, 'naive_bayes', diags_file)
            expected_lines = [decode(x) for x in gzip.open(expected_result, 'r').readlines()]
            expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
            produced_lines = [x.decode() for x in open(temp_output_path).readlines()]
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            os.remove(temp_output_path)
            os.remove(diags_file)
            result = ''.join(diff)
            self.assertEqual(result, '')

    def test_integration_ls(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'sv_classifier')
        input = os.path.join(test_data_dir, 'reclass.test.vcf.gz')

        expected_result = os.path.join(test_data_dir, 'output.ls.vcf.gz')

        annot=os.path.join(test_data_dir, 'repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz')
        sex_file=os.path.join(test_data_dir, 'ceph.sex.txt')
        train=os.path.join(test_data_dir, 'training.vars.vcf.gz')

        diags_handle, diags_file = tempfile.mkstemp(suffix='.txt')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        sex=open(sex_file, 'r')
        sex_chrom_names = set(('X', 'Y'))
        with gzip.open(input, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.sv_classifier.run_reclassifier(input_handle, output_handle, sex, sex_chrom_names, annot,  0.9, None, 1.0, 0.2, train, 'large_sample', diags_file)
            expected_lines = [decode(x) for x in gzip.open(expected_result, 'r').readlines()]
            expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
            produced_lines = [decode(x) for x in open(temp_output_path).readlines()]
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            os.remove(temp_output_path)
            os.remove(diags_file)
            result = ''.join(diff)
            self.assertEqual(result, '')



    def test_integration_hyb(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'sv_classifier')
        input = os.path.join(test_data_dir, 'reclass.test.vcf.gz')
        expected_result = os.path.join(test_data_dir, 'output.hyb.vcf.gz')

        annot=os.path.join(test_data_dir, 'repeatMasker.recent.lt200millidiv.LINE_SINE_SVA.b37.sorted.bed.gz')
        sex_file=os.path.join(test_data_dir, 'ceph.sex.txt')
        train=os.path.join(test_data_dir, 'training.vars.vcf.gz')

        diags_handle, diags_file = tempfile.mkstemp(suffix='.txt')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        sex=open(sex_file, 'r')
        sex_chrom_names = set(('X', 'Y'))
        with gzip.open(input, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            svtools.sv_classifier.run_reclassifier(input_handle, output_handle, sex, sex_chrom_names, annot,  0.9, None, 1.0, 0.2, train, 'hybrid', diags_file)
            expected_lines = [decode(x) for x in gzip.open(expected_result, 'r').readlines()]
            expected_lines[1] = '##fileDate=' + time.strftime('%Y%m%d') + '\n'
            produced_lines = [decode(x) for x in open(temp_output_path).readlines()]
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            os.remove(temp_output_path)
            os.remove(diags_file)
            result = ''.join(diff)
            self.assertEqual(result, '')


if __name__ == "__main__":
    main()

