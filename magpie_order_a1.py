# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""magpie-order-test, AMENDMENT-1: calibrate the preprint's own two procedures.

T5 (procedure "magpie-insample-loglik"): first- versus second-order in-sample log-likelihood with
additive smoothing k in KS, averaged over k; fires if order 2 is higher.
T6 (procedure "magpie-zero-order-perm"): every second-order transition probability compared with a
zero-order permutation null (call tokens shuffled across all positions, sequence lengths and
START/END kept), two-tailed, Benjamini-Hochberg at 0.05; fires if any cell is significant.

Both are run on the real type-level data and in worlds that have NO second-order rule:
H0p (one pooled first-order chain), H0a (each group its own first-order chain) and H0b (H0a plus
looped chunks). Pre-registered in AMENDMENT-1.md. numpy only.

Usage: python magpie_order_a1.py --data-dir DIR --out-dir DIR [--worlds 999] [--perm-worlds 200]
       [--perms 10000]
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

import magpie_order as mo

KS = (1e-4, 1e-3, 1e-2, 0.1, 1.0)
SEED_A1 = 20261125
FDR = 0.05
CHUNK = 500


# ----------------------------------------------------------------------------- T5

def smoothed_loglik(seqs, K, k):
    """In-sample log-likelihoods with additive smoothing k over M = K + 1 outcomes (units + END).
    Returns (LL1, LL2_matched, LL2_triplets). LL1: order 1 over every pair (START -> first unit,
    ..., last unit -> END). LL2_matched: order 2 over the same events, (START, START) being the
    first context. LL2_triplets: order 2 over the events at t >= 1 only (the preprint's triplets
    Start-A-B ... A-B-Stop)."""
    M = K + 1
    c1 = np.zeros((M, M))
    c2 = np.zeros((M * M, M))
    ev = []
    for s in seqs:
        for t, (a, b, o) in enumerate(mo.events(s, K)):
            c1[b, o] += 1
            c2[a * M + b, o] += 1
            ev.append((t, a, b, o))
    r1 = c1.sum(axis=1)
    r2 = c2.sum(axis=1)
    ll1 = ll2m = ll2t = 0.0
    for t, a, b, o in ev:
        ll1 += math.log((c1[b, o] + k) / (r1[b] + k * M))
        x = math.log((c2[a * M + b, o] + k) / (r2[a * M + b] + k * M))
        ll2m += x
        if t >= 1:
            ll2t += x
    return ll1, ll2m, ll2t


def t5(seqs, K):
    rows = []
    for k in KS:
        l1, l2m, l2t = smoothed_loglik(seqs, K, k)
        rows.append({"k": k, "LL1": l1, "LL2_matched": l2m, "LL2_triplets": l2t})
    dm = float(np.mean([r["LL2_matched"] - r["LL1"] for r in rows]))
    dt = float(np.mean([r["LL2_triplets"] - r["LL1"] for r in rows]))
    return {"by_k": rows, "mean_gain_matched": dm, "mean_gain_triplets": dt,
            "fires_matched": bool(dm > 1e-9), "fires_triplets": bool(dt > 1e-9)}


# ----------------------------------------------------------------------------- T6

def bh(p, q):
    """Benjamini-Hochberg step-up: boolean mask of rejected hypotheses at level q."""
    p = np.asarray(p, dtype=float)
    m = p.size
    out = np.zeros(m, dtype=bool)
    if m == 0:
        return out
    order = np.argsort(p, kind="mergesort")
    ok = p[order] <= q * np.arange(1, m + 1) / m
    if ok.any():
        kmax = int(np.nonzero(ok)[0].max())
        out[order[:kmax + 1]] = True
    return out


def structure(seqs, K):
    """Token vector and, for every second-order event (t >= 1), the indices of (prev2, prev1,
    outcome) into the extended vector [tokens..., START, END] (both coded K)."""
    n = sum(len(s) for s in seqs)
    s_idx, e_idx = n, n + 1
    tokens, ia, ib, io = [], [], [], []
    pos = 0
    for s in seqs:
        L = len(s)
        for t in range(1, L + 1):
            ia.append(pos + t - 2 if t >= 2 else s_idx)
            ib.append(pos + t - 1)
            io.append(pos + t if t < L else e_idx)
        tokens.extend(int(u) for u in s)
        pos += L
    return (np.array(tokens, dtype=np.int64), np.array(ia, dtype=np.int64),
            np.array(ib, dtype=np.int64), np.array(io, dtype=np.int64))


def tri_counts(tok_mat, ia, ib, io, K):
    """tok_mat (R, n) token rows -> second-order counts (R, M*M, M)."""
    R = tok_mat.shape[0]
    M = K + 1
    ext = np.concatenate([tok_mat, np.full((R, 2), K, dtype=np.int64)], axis=1)
    code = (ext[:, ia] * M + ext[:, ib]) * M + ext[:, io]
    code = code + (np.arange(R, dtype=np.int64) * (M ** 3))[:, None]
    cnt = np.bincount(code.ravel(), minlength=R * M ** 3)
    return cnt.reshape(R, M * M, M).astype(float)


def zero_order_test(seqs, K, B, rng, labels=None, chunk=CHUNK):
    """The preprint's second-order significance test on type-level sequences. Returns counts of
    BH-significant cells (all outcome cells of every observed context) and, if labels are given,
    the list of significant cells."""
    M = K + 1
    tokens, ia, ib, io = structure(seqs, K)
    obs = tri_counts(tokens[None, :], ia, ib, io, K)[0]
    rs = obs.sum(axis=1)
    ctx = np.nonzero(rs > 0)[0]
    otp = obs[ctx] / rs[ctx][:, None]
    hi = np.zeros_like(otp)
    lo = np.zeros_like(otp)
    msum = np.zeros_like(otp)
    nv = np.zeros(len(ctx))
    done = 0
    with np.errstate(invalid="ignore", divide="ignore"):
        while done < B:
            r = min(chunk, B - done)
            mat = rng.permuted(np.tile(tokens, (r, 1)), axis=1)
            c = tri_counts(mat, ia, ib, io, K)[:, ctx, :]
            rr = c.sum(axis=2)
            valid = rr > 0
            v3 = valid[:, :, None]
            tp = c / np.where(valid, rr, 1.0)[:, :, None]
            hi += np.sum((tp >= otp[None] - 1e-12) & v3, axis=0)
            lo += np.sum((tp <= otp[None] + 1e-12) & v3, axis=0)
            msum += np.sum(np.where(v3, tp, 0.0), axis=0)
            nv += valid.sum(axis=0)
            done += r
    n = nv[:, None]
    p = np.minimum(1.0, 2.0 * np.minimum((1.0 + hi) / (1.0 + n), (1.0 + lo) / (1.0 + n)))
    mean = np.where(n > 0, msum / np.maximum(n, 1.0), np.nan)
    sig = bh(p.ravel(), FDR).reshape(p.shape)
    with np.errstate(invalid="ignore"):
        above = sig & (otp > mean)
        below = sig & (otp < mean)
    res = {"contexts": int(len(ctx)), "cells": int(p.size), "perms": int(B),
           "min_p": float(p.min()) if p.size else 1.0,
           "S2": int(sig.sum()), "S2_above": int(above.sum()), "S2_below": int(below.sum()),
           "fires": bool(sig.sum() >= 1)}
    if labels is not None:
        names = list(labels) + ["START/END"]
        cells = []
        for ci, cx in enumerate(ctx):
            a, b = divmod(int(cx), M)
            for o in range(M):
                if sig[ci, o]:
                    cells.append({"prev2": names[a], "prev1": names[b], "next": names[o],
                                  "observed": float(otp[ci, o]), "null_mean": float(mean[ci, o]),
                                  "p": float(p[ci, o])})
        res["significant_cells"] = cells
    return res


# ----------------------------------------------------------------------------- analysis

def _rate(x, n):
    return x / n if n else float("nan")


def grade(rate):
    if rate >= 0.30:
        return "UNSAFE"
    if rate <= 0.10:
        return "SAFE"
    return "MARGINAL"


def analyse(units, P, sets, worlds, perm_worlds, perms, seed=SEED_A1, log=print):
    rng = np.random.default_rng(seed)
    K = len(units)
    idx = {u: i for i, u in enumerate(units)}

    def enc(c):
        return tuple(idx[u] for u in c.split())

    groups = [g for g in sets if sets[g]]
    by_group = [[enc(c) for c in sets[g]] for g in groups]
    U = sorted({s for seqs in by_group for s in seqs})
    Penc = [enc(c) for c in P]
    res = {"amendment": 1, "units": units, "n_types_U": len(U), "n_types_P": len(Penc),
           "tokens_U": int(sum(len(s) for s in U)), "ks": list(KS), "fdr": FDR}
    res["T5_U"] = t5(U, K)
    res["T5_P"] = t5(Penc, K)
    t = time.time()
    res["T6_U"] = zero_order_test(U, K, perms, rng, labels=units)
    log("real: T5 matched gain %.4f fires %s; T6 S2 %d (%.0f s)" % (
        res["T5_U"]["mean_gain_matched"], res["T5_U"]["fires_matched"], res["T6_U"]["S2"],
        time.time() - t))

    nulls_def = [
        ("H0p_pooled_first_order", mo.make_spec([U]), mo.first_order_tables([U], K), False),
        ("H0a_first_order", mo.make_spec(by_group), mo.first_order_tables(by_group, K), False),
        ("H0b_first_order_loops", mo.make_spec(by_group),
         mo.first_order_tables([[s[:mo.min_period(s)] for s in seqs] for seqs in by_group], K),
         True),
    ]
    nulls = {}
    for name, spec, (S, T), loop in nulls_def:
        fm = ft = f6 = f6a = fails = 0
        gains_m, s2 = [], []
        t = time.time()
        for w in range(worlds):
            recs, f = mo.gen_world(spec, S, T, rng, loop)
            fails += f
            Uw = sorted({s for _, s in recs})
            r5 = t5(Uw, K)
            fm += int(r5["fires_matched"])
            ft += int(r5["fires_triplets"])
            gains_m.append(r5["mean_gain_matched"])
            if w < perm_worlds:
                r6 = zero_order_test(Uw, K, perms, rng)
                f6 += int(r6["fires"])
                f6a += int(r6["S2_above"] >= 1)
                s2.append(r6["S2"])
            if (w + 1) % 100 == 0:
                log("%s: %d worlds, %.0f s" % (name, w + 1, time.time() - t))
        nw = min(perm_worlds, worlds)
        s2a = np.asarray(s2, dtype=float)
        nulls[name] = {
            "worlds": worlds, "generation_failures": fails,
            "T5_fires_matched": fm, "T5_rate_matched": _rate(fm, worlds),
            "T5_fires_triplets": ft, "T5_rate_triplets": _rate(ft, worlds),
            "T5_gain_matched": {"mean": float(np.mean(gains_m)),
                                "q05": float(np.quantile(gains_m, 0.05)),
                                "q50": float(np.quantile(gains_m, 0.5)),
                                "q95": float(np.quantile(gains_m, 0.95))},
            "perm_worlds": nw, "T6_fires": f6, "T6_rate": _rate(f6, nw),
            "T6_fires_above": f6a,
            "T6_S2": {"mean": float(s2a.mean()), "q50": float(np.quantile(s2a, 0.5)),
                      "q95": float(np.quantile(s2a, 0.95)), "max": float(s2a.max())},
            "p_T6_S2_real": float((1 + np.sum(s2a >= res["T6_U"]["S2"])) / (1 + s2a.size)),
        }
    res["nulls"] = nulls
    h0p = nulls["H0p_pooled_first_order"]
    res["grades"] = {
        "magpie-insample-loglik": grade(h0p["T5_rate_matched"]),
        "magpie-insample-loglik-triplets": grade(h0p["T5_rate_triplets"]),
        "magpie-zero-order-perm": grade(h0p["T6_rate"]),
    }
    res["reproduced_on_real_U"] = {"T5_matched": res["T5_U"]["fires_matched"],
                                   "T5_triplets": res["T5_U"]["fires_triplets"],
                                   "T6": res["T6_U"]["fires"]}
    res["ledger"] = {
        "magpie-insample-loglik": {"worlds": h0p["worlds"], "hits": h0p["T5_fires_matched"]},
        "magpie-zero-order-perm": {"worlds": h0p["perm_worlds"], "hits": h0p["T6_fires"]},
    }
    return res


def write_md(res, path):
    n = res["nulls"]
    lines = ["# magpie-order-test, amendment 1 results", "",
             "Made with AI agents working in MV2's name (MV2, Munim, Inc.). CC BY 4.0.", "",
             "| item | value |", "|---|---|",
             "| types U / tokens | %d / %d |" % (res["n_types_U"], res["tokens_U"]),
             "| T5 real: mean gain matched / triplets (nats) | %.4f / %.4f |" % (
                 res["T5_U"]["mean_gain_matched"], res["T5_U"]["mean_gain_triplets"]),
             "| T6 real: significant cells (above / below) | %d (%d / %d) of %d |" % (
                 res["T6_U"]["S2"], res["T6_U"]["S2_above"], res["T6_U"]["S2_below"],
                 res["T6_U"]["cells"])]
    for name, v in n.items():
        lines.append("| %s: T5 rate matched / triplets; T6 rate | %.3f / %.3f; %.3f |" % (
            name, v["T5_rate_matched"], v["T5_rate_triplets"], v["T6_rate"]))
    for k, v in res["grades"].items():
        lines.append("| grade %s (pooled first-order worlds) | %s |" % (k, v))
    lines.append("")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--worlds", type=int, default=999)
    ap.add_argument("--perm-worlds", type=int, default=200)
    ap.add_argument("--perms", type=int, default=10000)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    units, P, sets = mo.build(mo.parse_rows(os.path.join(a.data_dir, mo.DATA_FILE)))
    t = time.time()
    res = analyse(units, P, sets, a.worlds, a.perm_worlds, a.perms,
                  log=lambda m: print(m, flush=True))
    res["seconds"] = time.time() - t
    with open(os.path.join(a.out_dir, "order_a1_results.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    write_md(res, os.path.join(a.out_dir, "RESULTS-A1.md"))
    print("grades", json.dumps(res["grades"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
