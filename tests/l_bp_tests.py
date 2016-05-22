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

    def test_trim(self):
        self.assertEqual(trim([0, 0, 2, 3, 5, 0, 0, 0]), [2, 3])
        self.assertEqual(trim([0, 0, 2, 3, 5]), [2, 0])
        self.assertEqual(trim([2, 3, 5, 0, 0, 0]), [0, 3])
        self.assertEqual(trim([2, 3, 5]), [0, 0])

    def test_align_intervals(self):
        interval1 = [0, 4, [0.2, 0.2, 0.2, 0.2, 0.2]]
        interval2 = [1, 5, [0.2, 0.2, 0.2, 0.2, 0.2]]
        self.assertEqual(align_intervals([interval1, interval2]), [0, 5, [[0.2, 0.2, 0.2, 0.2, 0.2, 0], [0, 0.2, 0.2, 0.2, 0.2, 0.2]]])


if __name__ == "__main__":
    main()

