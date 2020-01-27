from unittest import TestCase

from primers import primers, Primer
from primers.primers import _rc, _parse_add_len, LEN_MIN, LEN_MAX


class TestPrimers(TestCase):
    """Test primer creation."""

    def test_primers(self):
        """Create primers without additional sequence."""

        seq = "CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT"

        ps = primers(seq)
        p1, p2 = ps

        self.assertIsInstance(p1, Primer)
        self.assertIsInstance(p2, Primer)
        self.assertEqual(p1.seq, p1.seq.upper())
        self.assertEqual(p2.seq, p2.seq.upper())
        self.assertIn(p1.seq, seq)
        self.assertIn(_rc(p2.seq), seq)
        self.assertTrue(seq.startswith(p1.seq))
        self.assertTrue(seq.endswith(_rc(p2.seq)))
        self.assertTrue(LEN_MIN <= len(p1.seq) <= LEN_MAX)
        self.assertTrue(LEN_MIN <= len(p2.seq) <= LEN_MAX)
        self.assertTrue(p1.tm)
        self.assertTrue(p2.tm)
        self.assertEqual(p1.tm, p1.tm_total)
        self.assertEqual(p2.tm, p2.tm_total)
        self.assertTrue(p1.gc)
        self.assertTrue(p2.gc)
        self.assertTrue(p1.fwd)
        self.assertFalse(p2.fwd)
        self.assertFalse(p1.offtargets)
        self.assertFalse(p2.offtargets)
        self.assertTrue(p1.penalty)
        self.assertTrue(p2.penalty)

    def test_primers_add_left_fixed(self):
        """Create primers with additional sequence on the left, fixed."""

        seq = "CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT"
        add_fwd = "GGTCTC"

        ps = primers(seq, add_fwd=add_fwd)
        p1, p2 = ps

        self.assertTrue(p1.seq.startswith(add_fwd))
        self.assertTrue((add_fwd + seq).startswith(p1.seq))
        self.assertEqual(p2.tm, p2.tm_total)
        self.assertNotEqual(p1.tm, p1.tm_total)

    def test_primers_add_right_fixed(self):
        """Create primers with additional sequence on the right, fixed."""

        seq = "CTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGTCAGGTACGCTAGTATGCGTTCTTCCTTCCAGAGGTATGTGGCTGCGTGGTCAAAAGTGCGGCATTCGTATTTGCTCCTCGTGTTTACTCTCACAAACTTGACCTGGAGATCAAGGAGATGCTTCTTGTGGAACTGGACAACGCATCAACGCAACGGATCTACGTTACAGCGT"
        add_rev = "GGTCTC"

        ps = primers(seq, add_rev=add_rev)
        p1, p2 = ps

        self.assertTrue(p2.seq.startswith(add_rev))
        self.assertTrue((add_rev + _rc(seq)).startswith(p2.seq))
        self.assertEqual(p1.tm, p1.tm_total)
        self.assertNotEqual(p2.tm, p2.tm_total)

    def test_parse(self):
        """Throw an error on an invalid input sequence."""

        seq = "AUGGCCUUUCUCGGGGGGCCUUGGUCUCUGAAAC"

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
