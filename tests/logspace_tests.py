from unittest import TestCase
import svtools.logspace as ls
import math

class LogspaceTests(TestCase):
    def test_get_p(self):
        self.assertEqual(ls.get_p(0), 1.0)
        self.assertEqual(ls.get_p(1), math.e)

    def test_get_ls(self):
        self.assertEqual(ls.get_ls(0), float("-inf"))
        self.assertEqual(ls.get_ls(1), 0.0)

    def test_ls_multiply(self):
        x = float("-inf")
        y = 1
        self.assertEqual(ls.ls_multiply(x, y), float("-inf"))
        self.assertEqual(ls.ls_multiply(y, x), float("-inf"))
        self.assertEqual(ls.ls_multiply(y, y), 2.0)

    def test_ls_divide(self):
        x = 2
        y = 1
        self.assertEqual(ls.ls_divide(x, y), 1.0)
        self.assertEqual(ls.ls_divide(x, x), 0.0)

    def test_ls_add(self):
        x = float("-inf")
        y = 1
        z = 2
        self.assertEqual(ls.ls_add(x, y), 1.0)
        self.assertEqual(ls.ls_add(y, x), 1.0)
        # Using almost equal here as I think we are losing precision
        self.assertAlmostEqual(ls.ls_add(y, z), math.log(math.e + math.exp(2)))
        self.assertAlmostEqual(ls.ls_add(z, y), math.log(math.e + math.exp(2)))
