from unittest import TestCase, main
import sys
import os
import svtools.io.utils as sio

class InputStreamTest(TestCase):
    def test_init_hyphen(self):
        new_handle = sio.InputStream('-')
        self.assertIs(new_handle.handle, sys.stdin)

    def test_valid(self):
        self.assertTrue(sio.InputStream.valid('-'))
        with self.assertRaises(IOError):
            sio.InputStream.valid(None)

    def test_context_manager(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, '..', 'test_data', 'io', 'utils')
        test_input = os.path.join(test_data_dir, 'file.txt.gz')

        temporary_obj = None
        with sio.InputStream(test_input) as stream:
            temporary_obj = stream
            for line in stream:
                sys.stdout.write(line)
        self.assertTrue(temporary_obj.closed)




if __name__ == "__main__":
    main()

