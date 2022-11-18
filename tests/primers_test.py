from unittest import TestCase

from primers import primers, score, Primer
from primers.primers import _rc, _parse_add_len, LEN_MIN, LEN_MAX, _binding_seq


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
        self.assertFalse(p1.off_target_count)
        self.assertFalse(p2.off_target_count)
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

    def test_primers_parent(self):
        """Create primers given a parent with diff-case sequence."""

        ps = primers(
            "AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA",
            offtarget_check="ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga",
        )

        self.assertTrue(ps)

    def test_primers_score(self):
        """Score an existing pair of primers."""

        fwd, rev = score(
            "GGTCTCAATGAGACAATA",
            "TTTCGTATGCTGACCTAG",
            offtarget_check="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA",
        )
        self.assertAlmostEqual(fwd.penalty, 19, delta=3)
        self.assertAlmostEqual(rev.penalty, 6, delta=3)

        fwd_seq = "GGTCTCAATGAGACAATA"
        rev_seq = "AAAAAATTTCGTATGCTGACCTAG"
        fwd, rev = score(
            fwd_seq,
            rev_seq,
            seq="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA",
        )
        self.assertEqual(fwd_seq, fwd.seq)
        self.assertEqual(rev_seq, rev.seq)
        self.assertGreater(fwd.tm_total, fwd.tm)
        self.assertGreater(rev.tm_total, rev.tm)

    def test_primers_score_only_fwd(self):
        """Score an existing pair of primers."""

        fwd, rev = score(
            "GGTCTCAATGAGACAATAGCACACAC",
            offtarget_check="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA",
        )

        self.assertAlmostEqual(fwd.penalty, 11, delta=3)
        self.assertIsNone(rev)

    def test_primers_score_exceptions(self):
        """Score an existing pair of primers."""

        with self.assertRaises(ValueError):
            score("ACGACTAC")  # too short fwd

        with self.assertRaises(ValueError):
            score("ACGACTACGACTACGATC", "GACTACG")  # too short rev

    def test_primers_binding_sites(self):
        """Subsequence the seq parameter via binding sites"""

        add_fwd, seq, add_rev = _binding_seq(
            "GGTCTCAATGAGACAATA",
            rev="AAAAAATTTCGTATGC",
            seq="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT",
        )
        self.assertEqual(6, add_fwd)
        self.assertEqual(3, add_rev)
        self.assertEqual("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT", seq)

        # GGTCTCAATGAGACAATA
        #       AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT
        #                             CTAGGTCAGCATACGAAATTTTTT
        #       AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA <- goal

        add_fwd, seq, add_rev = _binding_seq(
            "GGTCTCAATGAGACAATA",
            rev="TTTCGTATGCTGACCTAG",
            seq="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT",
        )
        self.assertEqual(6, add_fwd)
        self.assertEqual("AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA", seq)
        self.assertEqual(0, add_rev)

        with self.assertRaises(ValueError):
            # fwd primer not in seq
            _binding_seq("GGTCTCAATGAGACAATA", seq="GCATCGATCTCATCTACGACTAGCAT")

        with self.assertRaises(ValueError):
            # rev primer not in seq
            _binding_seq(
                "GGTCTCAATGAGACAATA",
                rev="TTTCGTATCCCCCCTAG",
                seq="AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAATTT",
            )

    def test_primers_binding_sites_zero_index(self):
        """validate that the binding seq search works across the zero-index of a plasmid

        for: https://github.com/Lattice-Automation/primers/issues/4
        """

        add_fwd, seq, add_rev = _binding_seq(
            "CCCGTAGAAAAGATCAAAGGATCTTC",
            rev="GTAAAAAGGCCGCGTTGCTG",
            seq="CGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCT",
        )

        self.assertEqual(
            "CCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTAC",
            seq,
        )
        self.assertEqual(0, add_fwd)
        self.assertEqual(0, add_rev)
