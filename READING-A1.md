# magpie-order-test, amendment 1: graded reading

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human
expert has checked this reading yet. Text CC BY 4.0.

Job: magpie-order-a1-v1 (Blade, 2026-11-25 07:09-07:54 UTC, 2 CPUs, peak 0.04 GB). The gate
(`gate_a1.py`) passed, so AMENDMENT-1.md and its code on GitHub were identical to the job's copies
before anything was computed; `get_v1_code.py` checked the registered `magpie_order.py` and
`fetch_data.py` by git blob hash; the 10 unit tests passed. Every number below is copied from
[results/a1/order_a1_results.json](results/a1/order_a1_results.json) (summary table:
[results/a1/RESULTS-A1.md](results/a1/RESULTS-A1.md)).

## What was asked

The magpie preprint (Mason, Walsh, King and Ridley, bioRxiv 2026, doi:10.64898/2026.06.15.732496)
supports second-order call sequences with two procedures. AMENDMENT-1.md rebuilt both on the
public type-level data (117 distinct sequence types, 471 call tokens, 7 call units) and measured
how often each reports second-order structure in worlds that have none:

- **T5, in-sample log-likelihood** (procedure `magpie-insample-loglik`): order-2 against order-1
  log-likelihood on the same data, additive smoothing k from 0.0001 to 1, averaged over k.
- **T6, zero-order permutation** (procedure `magpie-zero-order-perm`): every second-order
  transition probability against 10,000 permutations of all call tokens, Benjamini-Hochberg at 0.05.

Calibration grade from the fire rate in H0p (one pooled first-order chain, no second-order rule, no
groups, no loops): UNSAFE if at least 0.30, SAFE if at most 0.10.

## Results

| | real data | H0p pooled first-order | H0a group first-order | H0b groups + loops |
|---|---|---|---|---|
| T5 fires (matched) | yes | 999 of 999 (1.000) | 999 of 999 (1.000) | 999 of 999 (1.000) |
| T5 fires (triplets only) | yes | 999 of 999 (1.000) | 999 of 999 (1.000) | 999 of 999 (1.000) |
| T5 mean gain, matched (nats) | 78.2856 | median 89.6912 (5%: 78.6704) | median 100.0106 (5%: 87.9852) | median 104.7240 (5%: 91.2963) |
| T6 fires | yes (S2 = 3: 1 above, 2 below) | 98 of 200 (0.490) | 187 of 200 (0.935) | 188 of 200 (0.940) |
| T6 significant cells S2 | 3 | median 0, 95%: 5, max 8 | median 5, 95%: 9, max 11 | median 5, 95%: 9, max 10 |
| share of worlds with S2 >= real 3 | | 0.388 | 0.886 | 0.886 |

**Grades (registered rule): `magpie-insample-loglik` UNSAFE (1.000), `magpie-insample-loglik-triplets`
UNSAFE (1.000), `magpie-zero-order-perm` UNSAFE (0.490).**

Both of the preprint's procedures reproduce on the real data: the in-sample comparison prefers
order 2 (gain 78.29 nats matched, 244.43 over triplets only) and the permutation test finds 3
significant second-order cells. Both also fire in worlds that contain no second-order rule at all:
the in-sample comparison in every one of 2,997 first-order worlds, the permutation test in about half
of the pooled first-order worlds and in more than 9 of 10 worlds built from each social group's own
first-order chain.

Two descriptive points, not registered tests:

1. The real in-sample gain (78.2856 nats) is **below the 5th percentile of every null**
   (78.6704 in H0p, 87.9852 in H0a, 91.2963 in H0b). On the measure the preprint reads as evidence
   for second order, the real sequence types show less apparent second-order structure than
   sequences drawn from first-order chains fitted to them. This matches v1, where the real pooled
   held-out gain was below the H0a 5th percentile and the within-group gain was 0.0000 nats.
2. Two of the three real significant cells have START as the call two back and are cells where the
   observed probability is 0 and the null mean about 0.25-0.27 (START, NLDS -> NLDS and START, SHDS
   -> SHDS). They say that a type never repeats its first call directly. A token shuffle breaks that
   first-order (no-direct-repeat) property, so the zero-order null flags it as a second-order cell.
   The third cell (NLDS, SHDS -> LH, observed 0.304 against a null mean of 0.048) is the only one
   above chance; the real count of 3 cells is typical of all three null families (p 0.388 to 0.886).

## Predictions (priors registered in AMENDMENT-1.md)

| id | prediction | prior | outcome | Brier |
|---|---|---|---|---|
| Q6 | T5 (matched) fires on the real U | 0.80 | true | 0.0400 |
| Q7 | T5 (matched) fire rate in H0p >= 0.30 | 0.60 | true (1.000) | 0.1600 |
| Q8 | T6 fires on the real U | 0.85 | true | 0.0225 |
| Q9 | T6 fire rate in H0p >= 0.30 | 0.80 | true (0.490) | 0.0400 |
| Q10 | p_T6_S2 <= 0.05 (real S2 beyond H0b worlds) | 0.35 | false (0.886) | 0.1225 |

Mean Brier 0.0770 over 5 predictions. The biggest miss was Q7: we feared heavy smoothing at k = 1
might pull the average gain below zero in sparse worlds; it never did (the smallest gains are still
near 79 nats).

## Reading

The amendment closes the gap that v1 left open. v1 used a stand-in (AIC with observed-transition
counting) and found it fired in 1,998 of 1,998 first-order worlds; the preprint's **own** two
procedures, rebuilt from its Methods, are no better calibrated: 1.000 and 0.490 in pooled
first-order worlds, 1.000 and 0.94 once group repertoires and loops are kept. Neither can separate
second-order sequences from first-order ones in data of this kind. Together with v1's calibrated
result (no within-group gain from the call two back, power 0.935), the type-level public data give
no support for the second-order claim, and the procedures that produced it would have reported it
whether or not it was true.

The verdict of PLAN.md (POOLING_EXPLAINS) is unchanged; this amendment adds the calibration of the
preprint's own methods. The contrast with the marmoset benchmark the preprint cites (Bosshard et al.
2024) is direct: that study's held-out test fired in 0 of 1,996 first-order worlds and its effect
survived (see marmoset-order-test).

## Limits

Type level only: the preprint's token-level sequences (1,892 productions with caller identity) and
frequency weights are not public, so its exact numbers cannot be reproduced and fire rates on token
data could differ (in-sample nested likelihoods cannot favour order 1 as k goes to 0 on any data,
so T5's rate is near 1 by construction). The permutation scheme follows the Methods text; the
authors' MATLAB code is not public. Across all 999 H0p worlds the generator used its fallback
sampler 17 times. T6 was run on the first 200 worlds of each null. The work was carried out with AI
agents working in MV2's name; no human expert has checked it yet.
