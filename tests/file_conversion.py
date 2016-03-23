import abc
import difflib
import glob
import os
import re
import sys
import tempfile
import shutil

class FileConversionBase(object):
    __metaclass__ = abc.ABCMeta

    @property
    def substitutions(self):
        # of the form (<pattern>, <replacement>) where pattern is a regular
        # expression and replacement is a string.
        return []

    @abc.abstractmethod
    def forward_convert(self, input_file, output_file):
        pass

    @abc.abstractproperty
    def test_data_directory_name(self):
        pass

    def input_file_path(self, test_name):
        return self._glob_for_one(self._data_dir, test_name, "input.*")

    def expected_output_file_path(self, test_name):
        return self._glob_for_one(self._data_dir, test_name, "expected.*")

    def test_forward_conversions(self):
        for test_name in self._test_names:
            print test_name
            self.convert_and_diff_output(self.forward_convert,
                    self.input_file_path(test_name),
                    self.expected_output_file_path(test_name))

    def _glob_for_one(self, *args):
        glob_arg = os.path.join(*args)
        input_files = glob.glob(glob_arg)
        if len(input_files) < 1:
            raise RuntimeError("Found no input-files with search: %s" %
                    glob_arg)
        elif len(input_files) > 1:
            raise RuntimeError("Found more than one input-file: %s" %
                    str(input_files))
        else:
            return input_files[0]

    @property
    def _data_dir(self):
        return self._relative_path('test_data', self.test_data_directory_name)

    def _relative_path(self, *args):
        here = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(here, *args)

    @property
    def _test_names(self):
        return os.listdir(self._data_dir)

    def convert_and_diff_output(self, conversion_fn, input_file_path,
            expected_output_file_path):
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')

        with open(input_file_path, 'r') as input_handle, os.fdopen(temp_descriptor, 'w') as output_handle:
            conversion_fn(input_handle, output_handle)

        if os.environ.get("SVTOOLS_REGENERATE_TEST_DATA", False):
            shutil.copyfile(temp_output_path, expected_output_file_path)

        self._compare_files(temp_output_path, expected_output_file_path)
        os.remove(temp_output_path)

    def _compare_files(self, got_file_path, expected_file_path):
        got_lines = self._get_filtered_lines(got_file_path)
        expected_lines = self._get_filtered_lines(expected_file_path)

        diff = difflib.unified_diff(got_lines, expected_lines,
            fromfile=expected_file_path, tofile=got_file_path)

        result = ''.join(diff)
        if result != '':
            for line in result:
                sys.stdout.write(line)
            self.assertTrue(False,
                    msg= "%s and %s differ!" %
                    (got_file_path, expected_file_path))
        return result

    def _get_filtered_lines(self, file_path):
        lines = [self._apply_substitutions(line)
                for line in open(file_path, 'r')]
        return lines

    def _apply_substitutions(self, line):
        for pattern, replacement in self.substitutions:
            line = re.sub(pattern, replacement, line)
        return line
