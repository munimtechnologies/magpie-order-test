# magpie-order-test, AMENDMENT-1 (2026-11-24, before any result)

MV2 of Munim, Inc. Made with AI agents working in his name; no human expert has checked this
amendment. Text CC BY 4.0, code Apache-2.0.

## Why an amendment

PLAN.md said that if the preprint's methods showed a procedure different from T0, an amendment
would add it before any result was read. They do. This amendment adds the preprint's **own two
procedures** and measures how often each one reports second-order structure in worlds that have
none. Everything in PLAN.md stays as registered and is still run by the job magpie-order-v1; this
amendment is a separate job (magpie-order-a1-v1) and changes nothing in v1.

## What has been seen (full disclosure)

- **No result of ours exists.** The v1 job is still queued; no output of it, and no statistic on the
  data, has been seen by anyone working on this study.
- A plain-text copy of the preprint (Mason, Walsh, King and Ridley, bioRxiv 2026,
  doi:10.64898/2026.06.15.732496, CC BY-NC-ND 4.0, not redistributed here) was made on 2026-11-24.
  We read its Methods in full. From its Results we saw, through a keyword search, sentences
  reporting: 11 first-order transitions above chance and 17 below in adults; many second-order
  transitions significant against the permutation null; and in-sample log-likelihood improving
  from zero to first to second order in both age groups.
- The preprint has no data-availability statement; its token-level data (1,892 sequence productions
  with caller identity) are not public. The type-level table from the group's Proc R Soc B paper
  (figshare 28251398 v2, CC BY 4.0) is still the only data used.

## The preprint's two procedures (from its Methods)

1. **Model comparison by in-sample log-likelihood.** First- and second-order transition matrices
   with additive (Dirichlet) smoothing k in {0.0001, 0.001, 0.01, 0.1, 1}, log-likelihood summed
   over all observed transitions, results averaged over k; a higher (less negative) second-order
   log-likelihood is read as second-order structure. No penalty for the extra parameters and no
   held-out data. For nested maximum-likelihood chains the second-order in-sample log-likelihood
   can never be lower than the first-order one, so as k goes to 0 this comparison cannot fail.
2. **Significance of second-order transitions against a zero-order null.** 10,000 permutations
   that rebuild sequences from the base frequencies of the calls, keeping Start/Stop and the
   sequence lengths; two-tailed permutation p for each transition probability; Benjamini-Hochberg
   at 0.05. A zero-order null has no first-order structure either, so a second-order cell can be
   "significant" purely because of first-order structure.

## Added tests

The real data are the distinct types U of PLAN.md (type level; the preprint's frequency weights
are not public, so every type counts once).

- **T5 (procedure magpie-insample-loglik).** For each k above: LL1 = order-1 log-likelihood over
  every pair (Start -> first call ... last call -> Stop); LL2 over the same events with
  (Start, Start) as the first context ("matched", primary), and LL2 over the triplets only
  (Start-A-B ... A-B-Stop, as the preprint's text describes; reported). Smoothing over the 8
  outcomes (7 calls and Stop). T5 fires if the mean over k of LL2 - LL1 is above 1e-9.
- **T6 (procedure magpie-zero-order-perm).** Second-order transition probabilities for every
  outcome cell of every observed context (prev2 may be Start; context (Start, Start) excluded).
  Null: 10,000 permutations shuffling all call tokens across all positions (lengths, Start and Stop
  kept; base frequencies kept exactly). Two-tailed p = min(1, 2 min(p_hi, p_lo)) with
  p_hi = (1 + #null >= obs) / (1 + n), counting only permutations in which the context occurs.
  Benjamini-Hochberg at 0.05 over all cells. T6 fires if at least one cell is significant; the
  number of significant cells S2 is also recorded. With 10,000 permutations the smallest p is
  about 0.0002, so BH can reject when two or more cells reach it (checked in the unit tests).

Null worlds (seed 20261125): 999 each for T5; the first 200 of them also get T6 (10,000
permutations each).

- **H0p, pooled first-order (new, primary for calibration):** each type of U is replaced by a
  sequence of the same length from ONE first-order chain fitted on U, distinct, without a
  directly repeated call when the real type has none (PLAN.md's generator with a single group).
  This world has no second-order rule, no groups and no loops.
- **H0a and H0b:** exactly as in PLAN.md.

## Readings fixed now

- Calibration grade of each procedure from its fire rate in H0p: **UNSAFE** if >= 0.30, **SAFE**
  if <= 0.10, otherwise MARGINAL. Rates in H0a and H0b are reported beside it.
- Reproduction: whether T5 (matched) and T6 fire on the real U.
- p_T6_S2: share of H0b worlds (with one added) whose S2 is at least the real S2. Descriptive
  only; the registered test of second order beyond pooling and loops remains T2 of PLAN.md.
- Ledger rows (MV2's calibrated-animal-claims): procedure "magpie-insample-loglik": worlds =
  nulls.H0p_pooled_first_order.worlds, hits = .T5_fires_matched; procedure
  "magpie-zero-order-perm": worlds = .perm_worlds, hits = .T6_fires (file order_a1_results.json).

## Predictions (priors, for a Brier score)

| id | prediction | prior |
|---|---|---|
| Q6 | T5 (matched) fires on the real U | 0.80 |
| Q7 | T5 (matched) fire rate in H0p >= 0.30 | 0.60 |
| Q8 | T6 fires on the real U | 0.85 |
| Q9 | T6 fire rate in H0p >= 0.30 | 0.80 |
| Q10 | p_T6_S2 <= 0.05 (real S2 beyond H0b worlds) | 0.35 |

Q7 is uncertain because type-level data are sparse: at k = 1 the second-order rows, which hold few
counts each, are smoothed much more than the first-order rows, which can pull the average gain
below zero.

## Limits

Type level only, as in PLAN.md; the preprint's own token-level data could give different fire
rates, and its frequency weighting cannot be reproduced. The permutation scheme follows the
Methods text; the authors' MATLAB code is not public, so details (for example how contexts absent
from a permutation are treated) are our reading. The generator's distinctness and
no-repeat conditions add slight dependence beyond first order in H0p.

## Code and run

magpie_order_a1.py (imports the registered magpie_order.py), test_magpie_order_a1.py,
get_v1_code.py (fetches magpie_order.py and fetch_data.py and checks their git blob SHA-1 against
the values registered on 2026-11-24), gate_a1.py (nothing is computed until this amendment and its
code on GitHub equal the job's copies). One CPU job, 1 GB.
