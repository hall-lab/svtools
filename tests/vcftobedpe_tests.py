from .file_conversion import FileConversionBase
from unittest import TestCase
import svtools.vcftobedpe

class BedpeToVcfTest(TestCase, FileConversionBase):
    @property
    def substitutions(self):
        return [('#fileDate=.*', '#fileDate=<DUMMY_DATE>')]

    def forward_convert(self, input_file, output_file):
        return svtools.vcftobedpe.vcfToBedpe(input_file, output_file)

    @property
    def test_data_directory_name(self):
        return 'vcftobedpe'
