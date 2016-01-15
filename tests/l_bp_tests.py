from unittest import TestCase, main
from l_bp import *

class Test_l_bp(TestCase):
    def test_find_all(self):
        test_string = 'ABBA'
        sub_string = 'A'

        result = [ x for x in find_all(test_string, sub_string) ]
        self.assertEqual(result, [0,3])

if __name__ == "__main__":
    main()

