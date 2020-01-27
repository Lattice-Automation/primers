from unittest import TestCase

from primers.offtargets import offtargets


class TestOfftargets(TestCase):
    """Test offtarget detection."""

    def test_offtargets(self):
        """Find and cache offtarget binding sites."""

        # GTGGCTAGCC is one by removed from GTGGCTAGGC in seq
        parent = "CTGACTCTACTTGGAAATGTGGCTAGGCCTTTGCCCACGCACCTGATCGGTCCTGTGGCTAGCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCATTGCAAGAATC"
        seq = "CTGACTCTACTTGGAAATGTGGCTAGGCCTT"

        ot = offtargets(seq, parent)

        self.assertEqual(0, ot[0])
        self.assertEqual(len(seq), len(ot))
        self.assertTrue(any(o for o in ot))
