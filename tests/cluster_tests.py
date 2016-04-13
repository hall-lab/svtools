from unittest import TestCase
from svtools.bedpe import Bedpe
from svtools.cluster import Cluster

class ClusterTests(TestCase):
    def test_can_add(self):
        bedpe = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b = Bedpe(bedpe)

        c = Cluster()
        c.chrom_a = b.c1
        c.chrom_b = b.c2
        c.min_a = b.s1
        c.max_a = b.e1
        c.min_b = b.s2
        c.max_b = b.e2
        c.strand_a = b.o1
        c.strand_b = b.o2

        self.assertTrue(c.can_add(b, 1))
        c.size = 1

        c.sv_event = 'DEL'
        self.assertFalse(c.can_add(b, 1))

        c.sv_event = 'BND'
        self.assertTrue(c.can_add(b, 1))

        c.chrom_a = 'X'
        self.assertFalse(c.can_add(b, 1))

        c.chrom_a = b.c1
        c.chrom_b = 'X'
        self.assertFalse(c.can_add(b, 1))


        c.chrom_b = b.c2
        c.min_a = 305
        self.assertFalse(c.can_add(b, 1))

        c.min_a = b.s1
        c.max_a = 150
        self.assertFalse(c.can_add(b, 1))

        c.max_a = b.e1
        c.min_b = 405
        self.assertFalse(c.can_add(b, 1))

        c.min_b = b.s1
        c.max_b = 150
        self.assertFalse(c.can_add(b, 1))

    def test_add(self):
        bedpe1 = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b1 = Bedpe(bedpe1)

        bedpe2= [ '1', '195', '305', '2', '295', '405', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.3' ]
        b2 = Bedpe(bedpe2)

        c = Cluster()
        c.add(b1, None)
        self.assertEqual(c.size, 1)
        self.assertEqual(c.sv_event, 'BND')
        self.assertEqual(c.filter, '0.2')
        self.assertEqual(c.chrom_a, '1')
        self.assertEqual(c.min_a, 200)
        self.assertEqual(c.max_a, 300)
        self.assertEqual(c.chrom_b, '2')
        self.assertEqual(c.min_b, 300)
        self.assertEqual(c.max_b, 400)
        self.assertEqual(c.strand_a, '+')
        self.assertEqual(c.strand_b, '-')

        c.add(b2, None)
        self.assertEqual(c.size, 2)
        self.assertEqual(c.sv_event, 'BND')
        self.assertEqual(c.filter, '0.3')
        self.assertEqual(c.chrom_a, '1')
        self.assertEqual(c.min_a, 195)
        self.assertEqual(c.max_a, 305)
        self.assertEqual(c.chrom_b, '2')
        self.assertEqual(c.min_b, 295)
        self.assertEqual(c.max_b, 405)
        self.assertEqual(c.strand_a, '+')
        self.assertEqual(c.strand_b, '-')

    def test_get_cluster_string(self):
        bedpe = [ '1', '200', '300', '2', '300', '400', '777_1', '57', '+', '-', 'BND', 'PASS', '.', '.', '.', '.', '.', '.', 'MISSING', 'SVTYPE=BND;AF=0.2' ]
        b = Bedpe(bedpe)

        c = Cluster()

        with self.assertRaises(ValueError):
            c.get_cluster_string()
        
        c.add(b, None)
        self.assertEqual(c.get_cluster_string(), str(b))

