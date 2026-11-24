# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""magpie-order-test: is the second-order structure reported for Western Australian magpie call
sequences more than (a) pooling group-specific repertoires and (b) looping a short chunk?
Pre-registered in PLAN.md. numpy only.

Usage: python magpie_order.py --data-dir DIR --out-dir DIR [--worlds 999] [--power-worlds 200]
"""
import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
ALPHA = 0.5
SEED = 20261124
RHOS = (0.3, 0.6)
DATA_FILE = "COMBINATION_repertoire.csv"


# ----------------------------------------------------------------------------- data

def _one(x):
    try:
        return float(str(x).strip()) == 1.0
    except ValueError:
        return False


def parse_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def build(rows):
    """units (sorted), P (sorted distinct type labels), sets {group: sorted labels the group's
    fledglings produced (fl01) or that were recorded from the natal group (group01)}."""
    P, sets = set(), {}
    for r in rows:
        call = " ".join(str(r.get("call", "")).split())
        if not call:
            continue
        P.add(call)
        g = str(r.get("group", "")).strip()
        sets.setdefault(g, set())
        if _one(r.get("fl01")) or _one(r.get("group01")):
            sets[g].add(call)
    units = sorted({u for c in P for u in c.split()})
    return units, sorted(P), {g: sorted(s) for g, s in sorted(sets.items())}


# ----------------------------------------------------------------------------- models

def events(seq, K):
    """(prev2, prev1, outcome) for each position; START = K as context, END = K as outcome."""
    out = []
    for t in range(len(seq) + 1):
        b = seq[t - 1] if t >= 1 else K
        a = seq[t - 2] if t >= 2 else K
        o = seq[t] if t < len(seq) else K
        out.append((a, b, o))
    return out


def mle_stats(seqs, K):
    """Maximum-likelihood first- and second-order chains on a list of sequences, with AIC/BIC
    under two parameter counts: 'observed' (non-zero outcomes minus one per seen context) and
    'full' (K free outcomes per seen context). T0 fires iff order-2 AIC (observed) is lower."""
    M = K + 1
    c1 = np.zeros((M, M))
    c2 = np.zeros((M * M, M))
    n = 0
    for s in seqs:
        for a, b, o in events(s, K):
            c1[b, o] += 1
            c2[a * M + b, o] += 1
            n += 1
    res = {"events": n, "sequences": len(seqs)}
    for name, c in (("order1", c1), ("order2", c2)):
        rows = c.sum(axis=1)
        ll, k_obs, seen = 0.0, 0, 0
        for i in range(c.shape[0]):
            if rows[i] > 0:
                v = c[i][c[i] > 0]
                ll += float(np.sum(v * np.log(v / rows[i])))
                k_obs += int(v.size) - 1
                seen += 1
        k_full = seen * K
        ln_n = math.log(n) if n > 0 else 0.0
        res[name] = {"loglik": ll, "contexts_seen": seen, "k_observed": k_obs, "k_full": k_full,
                     "aic_observed": 2 * k_obs - 2 * ll, "bic_observed": k_obs * ln_n - 2 * ll,
                     "aic_full": 2 * k_full - 2 * ll, "bic_full": k_full * ln_n - 2 * ll}
    o1, o2 = res["order1"], res["order2"]
    res["prefers_order2"] = {key: bool(o2[key] < o1[key]) for key in
                             ("aic_observed", "bic_observed", "aic_full", "bic_full")}
    res["fires"] = res["prefers_order2"]["aic_observed"]
    return res


def _counts(records, K, G):
    M = K + 1
    c1 = np.zeros((G, M, M))
    c2 = np.zeros((G, M * M, M))
    for g, s in records:
        for a, b, o in events(s, K):
            c1[g, b, o] += 1
            c2[g, a * M + b, o] += 1
    return c1, c2


def loo_components(records, K, G):
    """Leave-one-type-out predictive probabilities (add-ALPHA) for every event of every record:
    columns p1_group, p1_pool, p2_group, p2_pool. The held-out type is removed from its group
    (all copies there) and from the pool (all copies in every group)."""
    M = K + 1
    c1g, c2g = _counts(records, K, G)
    c1p, c2p = c1g.sum(axis=0), c2g.sum(axis=0)
    pool_w, grp_w = {}, {}
    for g, s in records:
        pool_w[s] = pool_w.get(s, 0) + 1
        grp_w[(g, s)] = grp_w.get((g, s), 0) + 1
    out = []
    for g, s in records:
        ev = events(s, K)
        d1 = np.zeros((M, M))
        d2 = np.zeros((M * M, M))
        for a, b, o in ev:
            d1[b, o] += 1
            d2[a * M + b, o] += 1
        g1 = c1g[g] - grp_w[(g, s)] * d1
        g2 = c2g[g] - grp_w[(g, s)] * d2
        p1 = c1p - pool_w[s] * d1
        p2 = c2p - pool_w[s] * d2
        for a, b, o in ev:
            ctx = a * M + b
            out.append(((g1[b, o] + ALPHA) / (g1[b].sum() + ALPHA * M),
                        (p1[b, o] + ALPHA) / (p1[b].sum() + ALPHA * M),
                        (g2[ctx, o] + ALPHA) / (g2[ctx].sum() + ALPHA * M),
                        (p2[ctx, o] + ALPHA) / (p2[ctx].sum() + ALPHA * M)))
    return np.array(out, dtype=float).reshape(-1, 4)


def gains(comp):
    """G_group: best cross-validated log-likelihood (nats) of the group-aware second-order model
    minus the best group-aware first-order model, mixing weights chosen on GRID (lambda = 0 is
    the first-order model, so G >= 0). G_pooled: the same without group terms."""
    p1g, p1p, p2g, p2p = comp[:, 0], comp[:, 1], comp[:, 2], comp[:, 3]
    L1 = {}
    for mu in GRID:
        L1[mu] = float(np.sum(np.log(mu * p1g + (1.0 - mu) * p1p)))
    mu1 = max(GRID, key=lambda m: L1[m])
    best = None
    for lam in GRID:
        for mu in GRID:
            v = float(np.sum(np.log(lam * p2g + (1.0 - lam) * (mu * p1g + (1.0 - mu) * p1p))))
            if best is None or v > best[0]:
                best = (v, lam, mu)
    lp1 = float(np.sum(np.log(p1p)))
    lp2, lamp = None, None
    for lam in GRID:
        v = float(np.sum(np.log(lam * p2p + (1.0 - lam) * p1p)))
        if lp2 is None or v > lp2:
            lp2, lamp = v, lam
    return {"events": int(comp.shape[0]),
            "L1_group": L1[mu1], "mu_order1": mu1,
            "L2_group": best[0], "lambda_order2": best[1], "mu_order2": best[2],
            "G_group": best[0] - L1[mu1],
            "L1_pooled": lp1, "L2_pooled": lp2, "lambda_pooled": lamp, "G_pooled": lp2 - lp1}


# ----------------------------------------------------------------------------- null worlds

def min_period(s):
    """Smallest p <= len/2 with s[i] == s[i-p] for all i >= p (a chunk looped at least twice);
    len(s) if none."""
    L = len(s)
    for p in range(1, L // 2 + 1):
        if all(s[i] == s[i - p] for i in range(p, L)):
            return p
    return L


def _row(own, pooled, exclude=None):
    for v in (own, pooled):
        if v.sum() > 0:
            return v / v.sum()
    u = np.ones(len(own))
    if exclude is not None and len(own) > 1:
        u[exclude] = 0.0
    return u / u.sum()


def first_order_tables(seqs_by_group, K):
    """Per group: start distribution over units and unit-to-unit transition rows (maximum
    likelihood, END dropped so lengths can be fixed), with pooled then uniform fallback."""
    G = len(seqs_by_group)
    st = np.zeros((G, K))
    tr = np.zeros((G, K, K))
    for g, seqs in enumerate(seqs_by_group):
        for s in seqs:
            st[g, s[0]] += 1
            for x, y in zip(s[:-1], s[1:]):
                tr[g, x, y] += 1
    pst, ptr = st.sum(axis=0), tr.sum(axis=0)
    S = np.zeros((G, K))
    T = np.zeros((G, K, K))
    for g in range(G):
        S[g] = _row(st[g], pst)
        for x in range(K):
            T[g, x] = _row(tr[g, x], ptr[x], exclude=x)
    return S, T


def sample_path(n, S, T, rng, rules=None):
    K = len(S)
    x = [int(rng.choice(K, p=S))]
    while len(x) < n:
        b = x[-1]
        a = x[-2] if len(x) >= 2 else K
        r = rules.get((a, b)) if rules is not None else None
        if r is None:
            r = T[b]
        x.append(int(rng.choice(K, p=r)))
    return x


def make_spec(seqs_by_group):
    return [[(len(s), min_period(s), any(s[i] == s[i - 1] for i in range(1, len(s))))
             for s in seqs] for seqs in seqs_by_group]


def gen_world(spec, S, T, rng, loop, rules=None, tries=200):
    """One null (or planted) world: every real type is replaced, group by group, by a sequence of
    the same length; with loop=True a chunk of the type's own period is drawn and looped, so the
    period is kept. Types stay distinct within a group and avoid self-transitions when the real
    type has none (up to `tries` draws; the count of failures is returned)."""
    recs, fails = [], 0
    for g, items in enumerate(spec):
        seen = set()
        for (L, p, has_self) in items:
            n = p if loop else L
            seq = None
            for _ in range(tries):
                chunk = sample_path(n, S[g], T[g], rng, rules)
                seq = tuple(chunk[i % n] for i in range(L))
                if not has_self and any(seq[i] == seq[i - 1] for i in range(1, L)):
                    continue
                if seq in seen:
                    continue
                break
            else:
                fails += 1
            seen.add(seq)
            recs.append((g, seq))
    return recs, fails


def planted_rules(K, rho, rng):
    """A random share rho of the second-order contexts (prev2 incl. START, prev1 a unit,
    prev2 != prev1) get their own next-unit distribution, Dirichlet(0.5) over units != prev1."""
    ctx = [(a, b) for a in range(K + 1) for b in range(K) if a != b]
    m = int(round(rho * len(ctx)))
    pick = rng.choice(len(ctx), size=m, replace=False)
    rules = {}
    for i in pick:
        a, b = ctx[int(i)]
        w = rng.dirichlet(np.full(K - 1, 0.5))
        row = np.zeros(K)
        row[[u for u in range(K) if u != b]] = w
        rules[(a, b)] = row / row.sum()
    return rules


def return_rate(seqs, K):
    """Share of unit positions t >= 2 whose unit equals the unit two back (A B -> A)."""
    hit = tot = 0
    for s in seqs:
        for t in range(2, len(s)):
            tot += 1
            hit += int(s[t] == s[t - 2])
    return hit / tot if tot else float("nan")


def pvalue(null, obs):
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null >= obs - 1e-9)) / (1 + null.size))


def world_stats(recs, K, G):
    U = sorted({s for _, s in recs})
    t0 = mle_stats(U, K)
    gn = gains(loo_components(recs, K, G))
    return {"fires": t0["fires"], "G_group": gn["G_group"], "G_pooled": gn["G_pooled"],
            "return_rate": return_rate(U, K)}


def _summ(vals):
    v = np.asarray(vals, dtype=float)
    return {"n": int(v.size), "mean": float(v.mean()), "q05": float(np.quantile(v, 0.05)),
            "q50": float(np.quantile(v, 0.5)), "q95": float(np.quantile(v, 0.95)),
            "max": float(v.max())}


def verdict(reproduced, p_t1, p_t2, power):
    if p_t2 <= 0.05:
        v = "SECOND_ORDER_SURVIVES"
    elif p_t1 <= 0.05:
        v = "LOOPING_EXPLAINS" if power >= 0.8 else "BEYOND_POOLING_LOOPING_UNDECIDED"
    else:
        v = "POOLING_EXPLAINS" if power >= 0.8 else "UNDECIDED_LOW_POWER"
    return v if reproduced else v + "_T0_NOT_REPRODUCED"


# ----------------------------------------------------------------------------- analysis

def analyse(units, P, sets, worlds, power_worlds, seed=SEED, log=print):
    rng = np.random.default_rng(seed)
    K = len(units)
    idx = {u: i for i, u in enumerate(units)}

    def enc(c):
        return tuple(idx[u] for u in c.split())

    groups = [g for g in sets if sets[g]]
    G = len(groups)
    by_group = [[enc(c) for c in sets[g]] for g in groups]
    D = [(gi, s) for gi, seqs in enumerate(by_group) for s in seqs]
    U = sorted({s for _, s in D})
    Penc = [enc(c) for c in P]
    spec = make_spec(by_group)
    periodic = [s for s in U if min_period(s) < len(s)]
    loop_pos = sum(len(s) - min_period(s) for s in periodic)

    res = {"units": units, "groups": groups, "group_sizes": [len(x) for x in by_group],
           "n_types_P": len(Penc), "n_types_U": len(U), "n_records_D": len(D),
           "unmatched_types_not_in_any_group": sorted(set(P) - {" ".join(units[i] for i in s)
                                                               for s in U}),
           "periodic_types_U": len(periodic),
           "loop_positions_share_U": loop_pos / max(1, sum(len(s) for s in U)),
           "self_transitions_U": sum(any(s[i] == s[i - 1] for i in range(1, len(s))) for s in U),
           "return_rate_U": return_rate(U, K)}
    res["T0_P"] = mle_stats(Penc, K)
    res["T0_U"] = mle_stats(U, K)
    res["cv_P"] = gains(loo_components([(0, s) for s in Penc], K, 1))
    res["cv_D"] = gains(loo_components(D, K, G))
    obs = res["cv_D"]["G_group"]
    log("real: T0_U fires %s, G_group %.4f, G_pooled %.4f" % (res["T0_U"]["fires"], obs,
                                                             res["cv_D"]["G_pooled"]))

    tables = {"H0a_first_order": first_order_tables(by_group, K),
              "H0b_first_order_loops": first_order_tables(
                  [[s[:min_period(s)] for s in seqs] for seqs in by_group], K)}
    nulls = {}
    for name, loop in (("H0a_first_order", False), ("H0b_first_order_loops", True)):
        S, T = tables[name]
        rows, fails, t_start = [], 0, time.time()
        for w in range(worlds):
            recs, f = gen_world(spec, S, T, rng, loop)
            fails += f
            rows.append(world_stats(recs, K, G))
            if (w + 1) % 100 == 0:
                log("%s: %d worlds, %.0f s" % (name, w + 1, time.time() - t_start))
        gvals = [r["G_group"] for r in rows]
        nulls[name] = {"worlds": worlds, "generation_failures": fails,
                       "T0_fires": int(sum(r["fires"] for r in rows)),
                       "T0_fire_rate": sum(r["fires"] for r in rows) / worlds,
                       "G_group": _summ(gvals),
                       "G_pooled": _summ([r["G_pooled"] for r in rows]),
                       "return_rate": _summ([r["return_rate"] for r in rows]),
                       "p_G_group": pvalue(gvals, obs),
                       "_g": gvals}
    res["nulls"] = {k: {kk: vv for kk, vv in v.items() if kk != "_g"} for k, v in nulls.items()}

    q95 = float(np.quantile(nulls["H0b_first_order_loops"]["_g"], 0.95))
    S, T = tables["H0b_first_order_loops"]
    power = {}
    for rho in RHOS:
        hits, fails, fires = 0, 0, 0
        for _ in range(power_worlds):
            rules = planted_rules(K, rho, rng)
            recs, f = gen_world(spec, S, T, rng, True, rules)
            fails += f
            st = world_stats(recs, K, G)
            hits += int(st["G_group"] > q95)
            fires += int(st["fires"])
        power[str(rho)] = {"worlds": power_worlds, "hits": hits, "rate": hits / power_worlds,
                           "T0_fires": fires, "generation_failures": fails}
        log("power rho %.1f: %d/%d" % (rho, hits, power_worlds))
    res["power"] = {"threshold_G_q95_H0b": q95, "by_rho": power}

    p1 = res["nulls"]["H0a_first_order"]["p_G_group"]
    p2 = res["nulls"]["H0b_first_order_loops"]["p_G_group"]
    pw = power[str(RHOS[0])]["rate"]
    res["tests"] = {"T0_reproduced": bool(res["T0_U"]["fires"]), "p_T1_beyond_pooling": p1,
                    "p_T2_beyond_pooling_and_loops": p2, "power_rho_0.3": pw}
    res["verdict"] = verdict(res["tests"]["T0_reproduced"], p1, p2, pw)
    res["ledger"] = {"procedure": "magpie-markov-aic", "worlds_H0b": worlds,
                     "fires_H0b": res["nulls"]["H0b_first_order_loops"]["T0_fires"],
                     "worlds_H0a": worlds, "fires_H0a": res["nulls"]["H0a_first_order"]["T0_fires"]}
    return res


def write_md(res, path):
    t0, cv, n = res["T0_U"], res["cv_D"], res["nulls"]
    lines = ["# magpie-order-test results", "",
             "Made with AI agents working in MV2's name (MV2, Munim, Inc.). CC BY 4.0.", "",
             "Verdict: **%s**" % res["verdict"], "",
             "| item | value |", "|---|---|",
             "| types (all labels / in some group) | %d / %d |" % (res["n_types_P"], res["n_types_U"]),
             "| group records | %d in %d groups |" % (res["n_records_D"], len(res["groups"])),
             "| looped types | %d |" % res["periodic_types_U"],
             "| T0 order-2 preferred by AIC (observed k) | %s |" % t0["fires"],
             "| G_group (nats) | %.4f |" % cv["G_group"],
             "| p T1 (vs group first-order) | %.4f |" % n["H0a_first_order"]["p_G_group"],
             "| p T2 (vs group first-order + loops) | %.4f |" % n["H0b_first_order_loops"]["p_G_group"],
             "| T0 fire rate, H0a / H0b | %.3f / %.3f |" % (n["H0a_first_order"]["T0_fire_rate"],
                                                         n["H0b_first_order_loops"]["T0_fire_rate"]),
             "| power rho 0.3 / 0.6 | %.3f / %.3f |" % (res["power"]["by_rho"]["0.3"]["rate"],
                                                     res["power"]["by_rho"]["0.6"]["rate"]),
             ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--worlds", type=int, default=999)
    ap.add_argument("--power-worlds", type=int, default=200)
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    units, P, sets = build(parse_rows(os.path.join(a.data_dir, DATA_FILE)))
    t = time.time()
    res = analyse(units, P, sets, a.worlds, a.power_worlds, log=lambda m: print(m, flush=True))
    res["seconds"] = time.time() - t
    with open(os.path.join(a.out_dir, "order_results.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    write_md(res, os.path.join(a.out_dir, "RESULTS-ORDER.md"))
    print("verdict", res["verdict"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
