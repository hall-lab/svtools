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
        self.assertEqual(bp.left.chrom, '1')
        self.assertEqual(bp.right.chrom, '10')
        self.assertEqual(bp.strands, '++:5')
        self.assertEqual(bp.left.start, 9572383 - 11)
        self.assertEqual(bp.left.end, 9572383 + 11)
        self.assertEqual(bp.right.start, 94079366 - 11)
        self.assertEqual(bp.right.end, 94079366 + 11)
        self.assertEqual(bp.left.p, self.prpos)
        self.assertEqual(bp.right.p, self.prend)

        # This was previously implemented in l_bp_tests, adding in here too
        test_line = '1	1000	2345_1	N	[2:1100[N	0.00	.	SVTYPE=BND;STRANDS=--:7;IMPRECISE;CIPOS=-2,2;CIEND=-2,2;CIPOS95=-1,1;CIEND95=-1,1;MATEID=2345_2;EVENT=2345;SU=7;PE=7;SR=0;PRPOS=0.025,0.25,0.45,0.25,0.025;PREND=0.025,0.25,0.45,0.25,0.025'
        no_slop = Breakpoint(test_line)
        self.assertEqual(no_slop.left.p, [0.025, 0.25, 0.45, 0.25, 0.025])
        self.assertEqual(no_slop.right.p, [0.025, 0.25, 0.45, 0.25, 0.025])

        fixed_slop = Breakpoint(test_line, fixed_slop = 1)
        self.assertEqual(fixed_slop.left.p, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])
        self.assertEqual(fixed_slop.right.p, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])

        percent_slop = Breakpoint(test_line, percent_slop = 0.2)
        print percent_slop
        self.assertEqual(percent_slop.left.p, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])
        self.assertEqual(percent_slop.right.p, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])

        percent_and_fixed_slop = Breakpoint(test_line, percent_slop = 0.2, fixed_slop = 2)
        self.assertEqual(percent_and_fixed_slop.left.p, [1e-100, 1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100, 1e-100])
        self.assertEqual(percent_and_fixed_slop.right.p, [1e-100, 1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100, 1e-100])

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

    def test_ovl(self):
        bp = Breakpoint(self.entry, fixed_slop=1)
        bp2 = Breakpoint(self.entry, fixed_slop=2)
        # Note that this is a regression test. This value was arrived at using the existing code.
        # It's correctness is unknown.
        self.assertEqual(bp.ovl(bp2), 1.0)

    def test_floats_from_tag(self):
        bp = Breakpoint(self.entry, fixed_slop=1)
        info = { 'TAG': '1.2,1.3'}
        self.assertEqual(bp.floats_from_tag(info, 'TAG'), [1.2, 1.3])
        with self.assertRaises(RuntimeError):
            bp.floats_from_tag(info, 'AG')


