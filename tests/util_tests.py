from unittest import TestCase, main
import sys
import os
import svtools.utils as su 

class InputStreamTest(TestCase):
    def test_init_hyphen(self):
        new_handle = su.InputStream('-')
        self.assertIs(new_handle.handle, sys.stdin)

    def test_init_stdin(self):
        new_handle = su.InputStream('stdin')
        self.assertIs(new_handle.handle, sys.stdin)

    def test_valid(self):
        self.assertTrue(su.InputStream.valid('-'))
        with self.assertRaises(IOError):
            su.InputStream.valid(None)

    def test_context_manager(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'utils')
        test_input = os.path.join(test_data_dir, 'file.txt.gz')

        temporary_obj = None
        with su.InputStream(test_input) as stream:
            temporary_obj = stream
            for line in stream:
                sys.stdout.write(line)
        self.assertTrue(temporary_obj.closed)

    def test_plain_iteration(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'utils')
        test_input = os.path.join(test_data_dir, 'file.txt')

        stream = su.InputStream(test_input)
        for line in stream:
            sys.stdout.write(line)
        stream.close()
        self.assertTrue(stream.handle.closed)


if __name__ == "__main__":
    main()

