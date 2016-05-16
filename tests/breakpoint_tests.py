from unittest import TestCase
from svtools.breakpoint import Breakpoint

class BreakpointTests(TestCase):
    def setUp(self):
        self.entry = "1	9572383	4921_1	N	N]10:94079366]	97.89	.	SVTYPE=BND;STRANDS=++:5;CIPOS=-10,10;CIEND=-10,10;CIPOS95=0,0;CIEND95=0,0;MATEID=4921_2;EVENT=4921;SU=5;PE=0;SR=5;PRPOS=9.80198e-21,9.80198e-19,9.80198e-17,9.80198e-15,9.80198e-13,9.80198e-11,9.80198e-09,9.80198e-07,9.80198e-05,0.00980198,0.980198,0.00980198,9.80198e-05,9.80198e-07,9.80198e-09,9.80198e-11,9.80198e-13,9.80198e-15,9.80198e-17,9.80198e-19,9.80198e-21;PREND=9.80198e-21,9.80198e-19,9.80198e-17,9.80198e-15,9.80198e-13,9.80198e-11,9.80198e-09,9.80198e-07,9.80198e-05,0.00980198,0.980198,0.00980198,9.80198e-05,9.80198e-07,9.80198e-09,9.80198e-11,9.80198e-13,9.80198e-15,9.80198e-17,9.80198e-19,9.80198e-21;SNAME=H_OY-FR97_8708-FR97_8708"
        prpos = [1e-100, 9.80198e-21, 9.80198e-19, 9.80198e-17, 9.80198e-15, 9.80198e-13, 9.80198e-11, 9.80198e-09, 9.80198e-07, 9.80198e-05, 0.00980198, 0.980198, 0.00980198, 9.80198e-05, 9.80198e-07, 9.80198e-09, 9.80198e-11, 9.80198e-13, 9.80198e-15, 9.80198e-17, 9.80198e-19, 9.80198e-21, 1e-100]
        sum_prpos = sum(prpos)
        self.prpos = [float(x) / sum_prpos for x in prpos]
        prend = [1e-100, 9.80198e-21, 9.80198e-19, 9.80198e-17, 9.80198e-15, 9.80198e-13, 9.80198e-11, 9.80198e-09, 9.80198e-07, 9.80198e-05, 0.00980198, 0.980198, 0.00980198, 9.80198e-05, 9.80198e-07, 9.80198e-09, 9.80198e-11, 9.80198e-13, 9.80198e-15, 9.80198e-17, 9.80198e-19, 9.80198e-21, 1e-100]
        sum_prend = sum(prend)
        self.prend = [float(x) / sum_prend for x in prend]
    
    def test_init(self):
        bp = Breakpoint(self.entry, fixed_slop=1)
        self.assertEqual(bp.l, self.entry)
        self.assertEqual(bp.sv_type, 'BND')
        self.assertEqual(bp.chr_l, '1')
        self.assertEqual(bp.chr_r, '10')
        self.assertEqual(bp.strands, '++:5')
        self.assertEqual(bp.start_l, 9572383 - 11)
        self.assertEqual(bp.end_l, 9572383 + 11)
        self.assertEqual(bp.start_r, 94079366 - 11)
        self.assertEqual(bp.end_r, 94079366 + 11)
        self.assertEqual(bp.p_l, self.prpos)
        self.assertEqual(bp.p_r, self.prend)

    def test_str(self):
        bp = Breakpoint(self.entry, fixed_slop=1)
        expected = [
                '1', 
                str(9572383 - 11),
                str(9572383 + 11),
                '10',
                str(94079366 - 11),
                str(94079366 + 11),
                'BND',
                '++:5',
                str(self.prpos),
                str(self.prend)
                ]

        self.assertEqual(str(bp), '\t'.join(expected))
