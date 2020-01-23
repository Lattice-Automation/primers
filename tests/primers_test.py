from unittest import TestCase

from primers import primers, Primer
from primers.primers import _rc, _parse_add_len


class TestPrimers(TestCase):
    """Test primer creation."""

    def test_primers(self):
        """Create primers without additional sequence."""

        seq = "CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT"

        ps = primers(seq)
        p1, p2 = ps

        self.assertIsInstance(p1, Primer)
        self.assertIsInstance(p2, Primer)
        self.assertIn(p1.seq, seq)
        self.assertIn(_rc(p2.seq), seq)
        self.assertTrue(p1.dg)
        self.assertTrue(p2.dg)
        self.assertTrue(p1.gc)
        self.assertTrue(p2.gc)
        self.assertTrue(p1.fwd)
        self.assertFalse(p2.fwd)
        self.assertTrue(p1.penalty)
        self.assertTrue(p2.penalty)
        self.assertEqual(p1.fwd, p1.fwd.upper())
        self.assertEqual(p2.fwd, p2.fwd.upper())

    def test_parse(self):
        """Throw an error on an invalid input sequence."""

        seq = "AUGGCCUUUCUCC"

        with self.assertRaises(ValueError):
            primers(seq)

    def test_parse_add_len(self):
        """Create a min/max added bp length range based on input."""

        s = "ATGGTAT"
        n = len(s)

        l1 = _parse_add_len(s, (-1, -1))
        l2 = _parse_add_len(s, (5, -1))
        l3 = _parse_add_len(s, (-1, 5))

        self.assertEqual((n, n), l1)
        self.assertEqual((5, n), l2)
        self.assertEqual((0, 5), l3)
