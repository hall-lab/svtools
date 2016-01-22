from unittest import TestCase, main
from svtools.l_bp import *

class Test_l_bp(TestCase):
    def test_find_all(self):
        test_string = 'ABBA'
        sub_string = 'A'

        result = [ x for x in find_all(test_string, sub_string) ]
        self.assertEqual(result, [0,3])

    def test_to_map(self):
        string = 'NS=3;AF=0.5;DB'
        expected = { 'NS' : '3', 'AF' : '0.5', 'DB' : None }
        result = to_map(string)
        self.assertEqual(result, expected)

    def test_split_v(self):
        info_bnd_map = {
                'SVTYPE' : 'BND',
                'STRANDS' : '-+',
                'CIPOS' : '-100,10',
                'CIEND' : '-10,400',
                'END' : '1100'
                }

        var1 = '1	1000	2345	N	]2:1100]N	0	.	SVTYPE=BND;STRANDS=-+;CIPOS=-100,10;CIEND=-10,400'
        self.assertEqual(split_v(var1), ['BND', '1', '2', '-+', 900, 1010, 1090, 1500, info_bnd_map])

        var2 = '1	1000	2345	N	[2:1100[N	0	.	SVTYPE=BND;STRANDS=-+;CIPOS=-100,10;CIEND=-10,400'
        self.assertEqual(split_v(var2), ['BND', '1', '2', '-+', 900, 1010, 1090, 1500, info_bnd_map])

        var3 = '1	1000	2345	N	<DEL>	0	.	SVTYPE=DEL;STRANDS=+-;END=1100;CIPOS=-100,10;CIEND=-10,400'
        info_bnd_map['SVTYPE'] = 'DEL'
        info_bnd_map['STRANDS'] = '+-'
        self.assertEqual(split_v(var3), ['DEL', '1', '1', '+-', 900, 1010, 1090, 1500, info_bnd_map])

    def test_breakpoint_init(self):
        test_line = '1	1000	2345_1	N	[2:1100[N	0.00	.	SVTYPE=BND;STRANDS=--:7;IMPRECISE;CIPOS=-2,2;CIEND=-2,2;CIPOS95=-1,1;CIEND95=-1,1;MATEID=2345_2;EVENT=2345;SU=7;PE=7;SR=0;PRPOS=0.025,0.25,0.45,0.25,0.025;PREND=0.025,0.25,0.45,0.25,0.025'
        no_slop = breakpoint(test_line)
        # not testing parts that take split_v values directly as that is tested above
        self.assertEqual(no_slop.p_l, [0.025, 0.25, 0.45, 0.25, 0.025])
        self.assertEqual(no_slop.p_r, [0.025, 0.25, 0.45, 0.25, 0.025])

        fixed_slop = breakpoint(test_line, fixed_slop = 1)
        self.assertEqual(fixed_slop.p_l, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])
        self.assertEqual(fixed_slop.p_r, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])

        percent_slop = breakpoint(test_line, percent_slop = 0.2)
        print percent_slop
        self.assertEqual(percent_slop.p_l, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])
        self.assertEqual(percent_slop.p_r, [1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100])

        percent_and_fixed_slop = breakpoint(test_line, percent_slop = 0.2, fixed_slop = 2)
        self.assertEqual(percent_and_fixed_slop.p_l, [1e-100, 1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100, 1e-100])
        self.assertEqual(percent_and_fixed_slop.p_r, [1e-100, 1e-100, 0.025, 0.25, 0.45, 0.25, 0.025, 1e-100, 1e-100])

if __name__ == "__main__":
    main()

