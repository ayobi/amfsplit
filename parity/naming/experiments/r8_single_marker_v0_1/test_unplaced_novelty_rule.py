#!/usr/bin/env python3
"""Unit tests for the experimental N6 unplaced-evidence gate."""
import unittest
from evaluate_unplaced_novelty_rule import experimental_decision


def row(*, status="novel lineage", rank="unplaced", matched=0, unmatched=3, no_call=0):
    return {
        "status": status,
        "novelty_claim": "candidate_only" if status.startswith("novel") else "not_asserted",
        "rank": rank,
        "matched_count": str(matched),
        "unmatched_count": str(unmatched),
        "no_call_count": str(no_call),
    }


class TestUnplacedNoveltyRule(unittest.TestCase):
    def test_all_three_unmatched_remains_candidate(self):
        status, claim, reason = experimental_decision(row(unmatched=3, no_call=0))
        self.assertEqual(status, "novel lineage")
        self.assertEqual(claim, "candidate_only")
        self.assertEqual(reason, "N6_AT_LEAST_TWO_UNMATCHED_MARKERS_NO_MATCHED_ARMS")

    def test_two_unmatched_one_no_call_remains_candidate(self):
        status, claim, _ = experimental_decision(row(unmatched=2, no_call=1))
        self.assertEqual((status, claim), ("novel lineage", "candidate_only"))

    def test_one_unmatched_two_no_call_is_insufficient(self):
        status, claim, reason = experimental_decision(row(unmatched=1, no_call=2))
        self.assertEqual((status, claim), ("insufficient evidence", "not_asserted"))
        self.assertEqual(reason, "N6_TOO_FEW_UNMATCHED_MARKERS_unmatched=1_no_call=2")

    def test_all_no_call_is_insufficient(self):
        status, claim, reason = experimental_decision(row(unmatched=0, no_call=3))
        self.assertEqual((status, claim), ("insufficient evidence", "not_asserted"))
        self.assertEqual(reason, "N6_TOO_FEW_UNMATCHED_MARKERS_unmatched=0_no_call=3")

    def test_matched_arm_conflict_does_not_assert_novelty(self):
        status, claim, reason = experimental_decision(row(matched=1, unmatched=2, no_call=0))
        self.assertEqual((status, claim), ("insufficient evidence", "not_asserted"))
        self.assertEqual(reason, "N6_CONFLICT_UNPLACED_BUT_MATCHED_ARM_PRESENT")

    def test_known_or_other_priority_status_is_preserved(self):
        original = row(status="divergent class", rank="unplaced", unmatched=3)
        status, claim, reason = experimental_decision(original)
        self.assertEqual(status, "divergent class")
        self.assertEqual(claim, "not_asserted")
        self.assertEqual(reason, "N6_NOT_APPLICABLE_BASELINE_PRESERVED")

    def test_non_unplaced_novel_lineage_does_not_get_rewritten(self):
        original = row(status="novel lineage", rank="genus", unmatched=1, no_call=2)
        status, claim, reason = experimental_decision(original)
        self.assertEqual((status, claim), ("novel lineage", "candidate_only"))
        self.assertEqual(reason, "N6_NOT_APPLICABLE_BASELINE_PRESERVED")


if __name__ == "__main__":
    unittest.main(verbosity=2)
