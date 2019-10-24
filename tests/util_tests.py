from unittest import TestCase, main
import sys
import os
import svtools.utils as su 

def decode(x):
    if isinstance(x, bytes):
        return x.decode()
    return x

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
                sys.stdout.write(decode(line))
        self.assertTrue(temporary_obj.closed)

    def test_plain_iteration(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        test_data_dir = os.path.join(test_directory, 'test_data', 'utils')
        test_input = os.path.join(test_data_dir, 'file.txt')

        stream = su.InputStream(test_input)
        for line in stream:
            sys.stdout.write(decode(line))
        stream.close()
        self.assertTrue(stream.handle.closed)

class ParseBndAltStringTest(TestCase):
    def test_bnd_alt_string(self):
        self.assertEqual(su.parse_bnd_alt_string('A[1:6['), ('[', '1', '6'))
        self.assertEqual(su.parse_bnd_alt_string('A]1:6]'), (']', '1', '6'))
        self.assertEqual(su.parse_bnd_alt_string(']1:6]A'), (']', '1', '6'))
        self.assertEqual(su.parse_bnd_alt_string(']HLA-DQB1*06:09:01:6]A'), (']', 'HLA-DQB1*06:09:01', '6'))
        with self.assertRaises(AssertionError):
            su.parse_bnd_alt_string(']1:6[A')
        with self.assertRaises(AssertionError):
            su.parse_bnd_alt_string('1')

if __name__ == "__main__":
    main()

