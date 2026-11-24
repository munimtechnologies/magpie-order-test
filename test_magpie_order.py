# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Unit tests for magpie_order.py; every expected number is traced by hand in the comments."""
import json
import math
import os
import tempfile
import unittest

import numpy as np

import magpie_order as mo


class Events(unittest.TestCase):
    def test_events(self):
        # START = END = K = 7: (S,S)->0, (S,0)->1, (0,1)->END
        self.assertEqual(mo.events((0, 1), 7), [(7, 7, 0), (7, 0, 1), (0, 1, 7)])

    def test_min_period(self):
        self.assertEqual(mo.min_period((0, 1, 0, 1)), 2)
        self.assertEqual(mo.min_period((0, 1, 2, 0, 1, 2, 0)), 3)
        self.assertEqual(mo.min_period((0, 1, 0)), 3)      # p must be <= len/2
        self.assertEqual(mo.min_period((0, 1, 0, 2)), 4)
        self.assertEqual(mo.min_period((5,)), 1)

    def test_pvalue(self):
        self.assertAlmostEqual(mo.pvalue([1, 2, 3], 2), 0.75)
        self.assertAlmostEqual(mo.pvalue([1, 2, 3], 5), 0.25)


class MLE(unittest.TestCase):
    def test_hand_example(self):
        # K=3, (0,1,0) and (2,1,2). Order 1: contexts S,0,1,2 each split 1/1 -> 4 * 2 ln .5.
        # Order 2: only (S,S) is split -> 2 ln .5; the other six contexts are deterministic.
        r = mo.mle_stats([(0, 1, 0), (2, 1, 2)], 3)
        self.assertEqual(r["events"], 8)
        self.assertAlmostEqual(r["order1"]["loglik"], 8 * math.log(0.5))
        self.assertAlmostEqual(r["order2"]["loglik"], 2 * math.log(0.5))
        self.assertEqual(r["order1"]["k_observed"], 4)
        self.assertEqual(r["order2"]["k_observed"], 1)
        self.assertEqual(r["order1"]["k_full"], 12)   # 4 contexts x 3
        self.assertEqual(r["order2"]["k_full"], 21)   # 7 contexts x 3
        self.assertAlmostEqual(r["order1"]["aic_observed"], 8 - 16 * math.log(0.5))
        self.assertAlmostEqual(r["order2"]["aic_observed"], 2 - 4 * math.log(0.5))
        self.assertTrue(r["fires"])

    def test_no_fire_when_equal(self):
        # (0,1) twice: every context deterministic, k_observed 0, AIC 0 = 0 -> does not fire
        r = mo.mle_stats([(0, 1), (0, 1)], 2)
        self.assertAlmostEqual(r["order1"]["loglik"], 0.0)
        self.assertAlmostEqual(r["order2"]["loglik"], 0.0)
        self.assertFalse(r["fires"])


class LOO(unittest.TestCase):
    def setUp(self):
        # K=2 (M=3 outcomes), one group: (0,1) twice and (1,0) once -> 9 events
        self.comp = mo.loo_components([(0, (0, 1)), (0, (0, 1)), (0, (1, 0))], 2, 1)

    def test_shape(self):
        self.assertEqual(self.comp.shape, (9, 4))

    def test_hand_values(self):
        # record 1, event (S,S)->0: after removing both copies, context S has only o=1 (count 1):
        # (0 + .5) / (1 + 1.5) = 0.2 for all four columns
        np.testing.assert_allclose(self.comp[0], [0.2, 0.2, 0.2, 0.2])
        # record 1, event (0,1)->END: order-1 context 1 keeps one 1->0 count: .5 / 2.5 = 0.2;
        # order-2 context (0,1) is empty after removal: .5 / 1.5
        np.testing.assert_allclose(self.comp[2], [0.2, 0.2, 1 / 3, 1 / 3])

    def test_gains_nonnegative(self):
        g = mo.gains(self.comp)
        self.assertGreaterEqual(g["G_group"], 0.0)
        self.assertGreaterEqual(g["G_pooled"], 0.0)
        self.assertEqual(g["events"], 9)


class Worlds(unittest.TestCase):
    def test_gen_world_keeps_length_and_period(self):
        # chunks (0,1) and (1,2,0) give the deterministic cycle 0->1->2->0
        by_group = [[(0, 1, 0, 1), (1, 2, 0)]]
        spec = mo.make_spec(by_group)
        self.assertEqual(spec, [[(4, 2, False), (3, 3, False)]])
        S, T = mo.first_order_tables([[s[:mo.min_period(s)] for s in g] for g in by_group], 3)
        np.testing.assert_allclose(S[0], [0.5, 0.5, 0.0])
        np.testing.assert_allclose(T[0], np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float))
        rng = np.random.default_rng(1)
        for _ in range(20):
            recs, _f = mo.gen_world(spec, S, T, rng, True)
            a, b = recs[0][1], recs[1][1]
            self.assertEqual(len(a), 4)
            self.assertEqual(len(b), 3)
            self.assertEqual(a[2], a[0])
            self.assertEqual(a[3], a[1])
            for s in (a, b):
                self.assertFalse(any(s[i] == s[i - 1] for i in range(1, len(s))))

    def test_planted_rules(self):
        rules = mo.planted_rules(3, 1.0, np.random.default_rng(2))
        self.assertEqual(len(rules), 9)   # (K+1)*K contexts minus K with prev2 == prev1
        for (a, b), row in rules.items():
            self.assertAlmostEqual(float(row.sum()), 1.0)
            self.assertEqual(row[b], 0.0)


class EndToEnd(unittest.TestCase):
    def test_build(self):
        head = ("ID,group,mother,closest_contact,closest_contact_FL0_A1,closest_adult_contact,"
                "total_number_of_contacts,proportion_spent_in_contact,call,fl01,group01,exgroup01")
        rows = ["F1,G1,m,c,0,a,3,0.5,A B,1,0,0", "F1,G1,m,c,0,a,3,0.5,B A,0,0,1",
                "F2,G2,m,c,0,a,3,0.5,A B,0,1,0", "F2,G2,m,c,0,a,3,0.5,B A,0,0,0"]
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.csv")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(head + "\n" + "\n".join(rows) + "\n")
            units, P, sets = mo.build(mo.parse_rows(p))
        self.assertEqual(units, ["A", "B"])
        self.assertEqual(P, ["A B", "B A"])
        self.assertEqual(sets, {"G1": ["A B"], "G2": ["A B"]})

    def test_analyse_smoke(self):
        sets = {"G1": ["A B", "A B A B", "B C A"], "G2": ["B C", "C A B", "A B C A"]}
        P = sorted({c for v in sets.values() for c in v} | {"C B"})
        res = mo.analyse(["A", "B", "C"], P, sets, 5, 3, seed=7, log=lambda m: None)
        self.assertEqual(res["n_types_P"], 7)
        self.assertEqual(res["n_types_U"], 6)
        self.assertEqual(res["unmatched_types_not_in_any_group"], ["C B"])
        self.assertEqual(res["periodic_types_U"], 1)
        for k in ("p_T1_beyond_pooling", "p_T2_beyond_pooling_and_loops"):
            self.assertTrue(0.0 < res["tests"][k] <= 1.0)
        self.assertIn(res["ledger"]["fires_H0b"], range(0, 6))
        self.assertIsInstance(res["verdict"], str)
        json.dumps(res)


if __name__ == "__main__":
    unittest.main()
