# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Unit tests for magpie_order_a1.py; every expected number is traced by hand in the comments."""
import json
import math
import unittest

import numpy as np

import magpie_order as mo
import magpie_order_a1 as a1


class LogLik(unittest.TestCase):
    def test_hand_example(self):
        # K=2, M=3, one sequence (0,1). Events: (S,S)->0, (S,0)->1, (0,1)->END.
        # Every order-1 and order-2 context holds exactly one count, so every event has
        # probability (1+k)/(1+3k) under both orders. LL1 = LL2_matched = 3 x, LL2_triplets = 2 x.
        for k in (1e-4, 0.1, 1.0):
            x = math.log((1 + k) / (1 + 3 * k))
            l1, l2m, l2t = a1.smoothed_loglik([(0, 1)], 2, k)
            self.assertAlmostEqual(l1, 3 * x)
            self.assertAlmostEqual(l2m, 3 * x)
            self.assertAlmostEqual(l2t, 2 * x)

    def test_t5_fire_rules(self):
        # matched gain is exactly 0 -> does not fire; triplets gain = -x > 0 -> fires
        r = a1.t5([(0, 1)], 2)
        self.assertEqual(r["mean_gain_matched"], 0.0)
        self.assertFalse(r["fires_matched"])
        self.assertTrue(r["fires_triplets"])
        self.assertEqual(len(r["by_k"]), 5)

    def test_order2_never_worse_near_ml(self):
        # nested maximum likelihood: in-sample order 2 >= order 1 when k -> 0
        rng = np.random.default_rng(3)
        for _ in range(20):
            seqs = [tuple(int(u) for u in rng.integers(0, 4, size=int(rng.integers(2, 7))))
                    for _ in range(15)]
            l1, l2m, _l2t = a1.smoothed_loglik(seqs, 4, 1e-12)
            self.assertGreaterEqual(l2m, l1 - 1e-6)


class BH(unittest.TestCase):
    def test_bh(self):
        # sorted .01 .03 .04 .2 vs .0125 .025 .0375 .05 -> only the first passes
        self.assertEqual(list(a1.bh([0.01, 0.04, 0.03, 0.2], 0.05)), [True, False, False, False])
        # step-up: .045 <= .05 at rank 3 rejects all three
        self.assertEqual(list(a1.bh([0.03, 0.04, 0.045], 0.05)), [True, True, True])
        self.assertEqual(list(a1.bh([0.01, 0.02, 0.9], 0.05)), [True, True, False])
        self.assertEqual(a1.bh([], 0.05).size, 0)


class Structure(unittest.TestCase):
    def test_structure_hand(self):
        # (0,1), K=2: tokens [0,1], START index 2, END index 3; events t=1: (S,0,1), t=2: (0,1,E)
        tok, ia, ib, io = a1.structure([(0, 1)], 2)
        self.assertEqual(list(tok), [0, 1])
        self.assertEqual(list(ia), [2, 0])
        self.assertEqual(list(ib), [0, 1])
        self.assertEqual(list(io), [1, 3])

    def test_counts_match_naive(self):
        rng = np.random.default_rng(5)
        K = 4
        M = K + 1
        seqs = [tuple(int(u) for u in rng.integers(0, K, size=int(rng.integers(2, 8))))
                for _ in range(25)]
        tok, ia, ib, io = a1.structure(seqs, K)
        got = a1.tri_counts(tok[None, :], ia, ib, io, K)[0]
        want = np.zeros((M * M, M))
        for s in seqs:
            for t, (a, b, o) in enumerate(mo.events(s, K)):
                if t >= 1:
                    want[a * M + b, o] += 1
        np.testing.assert_array_equal(got, want)
        # two rows: the second row is the same tokens again -> identical counts
        got2 = a1.tri_counts(np.vstack([tok, tok]), ia, ib, io, K)
        np.testing.assert_array_equal(got2[1], want)


class ZeroOrder(unittest.TestCase):
    def test_constant_tokens_never_significant(self):
        # all tokens 0: every permutation equals the data, so every p = 1 and nothing is significant
        r = a1.zero_order_test([(0, 0), (0, 0, 0)], 1, 300, np.random.default_rng(1), labels=["A"])
        self.assertEqual(r["S2"], 0)
        self.assertFalse(r["fires"])
        self.assertAlmostEqual(r["min_p"], 1.0)
        # contexts (S,0) and (0,0), each with M = 2 outcome cells
        self.assertEqual(r["contexts"], 2)
        self.assertEqual(r["cells"], 4)
        self.assertEqual(r["significant_cells"], [])

    def test_strong_structure_is_found(self):
        # 60 copies of the cycle (0,1,2): after 0-1 always comes 2, a permuted world rarely
        # repeats that; the smallest p is about 2/(B+1) and many cells share it, so BH rejects
        seqs = [(0, 1, 2)] * 30 + [(1, 2, 0)] * 30
        r = a1.zero_order_test(seqs, 3, 2000, np.random.default_rng(2), labels=["A", "B", "C"])
        self.assertTrue(r["fires"])
        self.assertGreaterEqual(r["S2_above"], 1)
        self.assertLess(r["min_p"], 0.01)

    def test_grade(self):
        self.assertEqual(a1.grade(0.30), "UNSAFE")
        self.assertEqual(a1.grade(0.10), "SAFE")
        self.assertEqual(a1.grade(0.2), "MARGINAL")


class Smoke(unittest.TestCase):
    def test_analyse_smoke(self):
        sets = {"G1": ["A B", "A B A B", "B C A"], "G2": ["B C", "C A B", "A B C A"]}
        P = sorted({c for v in sets.values() for c in v} | {"C B"})
        res = a1.analyse(["A", "B", "C"], P, sets, 4, 2, 200, seed=7, log=lambda m: None)
        self.assertEqual(res["n_types_U"], 6)
        self.assertEqual(res["n_types_P"], 7)
        self.assertEqual(set(res["nulls"]), {"H0p_pooled_first_order", "H0a_first_order",
                                             "H0b_first_order_loops"})
        for v in res["nulls"].values():
            self.assertEqual(v["worlds"], 4)
            self.assertEqual(v["perm_worlds"], 2)
            self.assertTrue(0.0 < v["p_T6_S2_real"] <= 1.0)
        self.assertEqual(res["ledger"]["magpie-zero-order-perm"]["worlds"], 2)
        json.dumps(res)


if __name__ == "__main__":
    unittest.main()
