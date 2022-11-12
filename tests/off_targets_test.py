from unittest import TestCase

from primers.off_targets import off_targets


class TestOffTargets(TestCase):
    """Test offtarget detection."""

    def test_off_targets(self):
        """Find and cache offtarget binding sites."""

        # GTGGCTAGCC is one by removed from GTGGCTAGGC in seq
        parent = "CTGACTCTACTTGGAAATGTGGCTAGGCCTTTGCCCACGCACCTGATCGGTCCTGTGGCTAGCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCATTGCAAGAATC"
        seq = "CTGACTCTACTTGGAAATGTGGCTAGGCCTT"

        ot = off_targets(seq, parent)

        self.assertEqual(0, ot[0])
        self.assertEqual(len(seq), len(ot))
        self.assertTrue(any(o for o in ot))
