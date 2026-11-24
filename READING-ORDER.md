# magpie-order-test: graded reading of job magpie-order-v1 (2026-11-24)

MV2 of Munim, Inc. Made with AI agents working in MV2's name; no human expert has checked this
reading yet. Text CC BY 4.0. Every number below is copied from
[results/order/order_results.json](results/order/order_results.json) and
[results/order/RESULTS-ORDER.md](results/order/RESULTS-ORDER.md), written by the job that ran the
code in this repository after `gate.py` confirmed that PLAN.md and the code on GitHub were identical
to the job's own copies (job log: same MD5 for PLAN.md, magpie_order.py, test_magpie_order.py and
fetch_data.py; data file COMBINATION_repertoire.csv, 78,432 bytes, MD5 f60a6af46f0161225b5965eb2c917874).

## Registered verdict: POOLING_EXPLAINS

The rule fixed in PLAN.md: p_T2 above 0.05 and p_T1 above 0.05, with power at rho 0.3 of at least
0.8, gives POOLING_EXPLAINS. The real data gave p_T1 = 1.0000, p_T2 = 1.0000 and power 0.935. The
published-style preference for order 2 was reproduced on the real types (T0 fires), so no
_T0_NOT_REPRODUCED suffix applies.

## What the job found

| item | value |
|---|---|
| sequence types / group records | 117 / 224 in 8 groups |
| T0: order 2 preferred by AIC (observed k) on the real types | yes (AIC 1718.10 against 1802.09) |
| same, BIC (observed k) / AIC and BIC (full k) | no / no / no |
| G_group, cross-validated group-aware gain of order 2 (nats) | 0.0000 (lambda chosen 0.0) |
| G_pooled on the 224 records, without group terms (nats) | 0.4359 |
| p T1 (against group first-order worlds, 999) | 1.0000 |
| p T2 (against group first-order worlds with loops, 999) | 1.0000 |
| T0 fires in H0a worlds / H0b worlds | 999 of 999 / 999 of 999 |
| power at rho 0.3 / 0.6 (200 worlds each) | 187 of 200 (0.935) / 200 of 200 |

1. **The published-style order test is not informative here.** Choosing between first- and
   second-order chains by AIC on the pooled types, counting only observed transitions as
   parameters, prefers order 2 on the real data, and it also prefers order 2 in every one of 999
   worlds built from each group's own first-order chain (H0a), and in every one of 999 such worlds
   with looped chunks kept as loops (H0b). Its false-alarm rate in worlds with no second-order rule
   is 1.000. BIC with the same count, and AIC or BIC with the full parameter count, do not prefer
   order 2 on the real data.
2. **Within groups, the call two back adds nothing.** With each held-out type predicted from its own
   group's first-order chain mixed with the pool, the best order-2 model chose lambda = 0 (no
   second-order term at all), so G_group = 0.0000 nats. Every null world gives G_group at least
   that large, so p = 1.0 for T1 and T2.
3. **The test could have seen a moderate second-order rule.** When 30% of the second-order contexts
   got their own next-call distribution, shared by all groups, G_group rose above the 95th
   percentile of the H0b worlds in 187 of 200 worlds (0.935); at 60% in 200 of 200.

Descriptive, not registered as a test: even without group terms, the cross-validated gain of order
2 on the 224 records is 0.4359 nats in total, below the 5th percentile of the same quantity in the
H0a worlds (2.2811; median 6.3923). In the type-level data there is less pooled second-order gain
than group-specific first-order chains would create on their own. The registered verdict name says
"pooling explains"; the plainer reading is that at the level of sequence types there is no
second-order dependence to explain, and the in-sample order test prefers order 2 whatever the
truth.

## Predictions (priors from PLAN.md) and Brier score

| id | prediction | prior | outcome | Brier |
|---|---|---|---|---|
| Q1 | T0 fires on the real U | 0.75 | yes | 0.0625 |
| Q2 | p_T1 <= 0.05 | 0.55 | no (1.0000) | 0.3025 |
| Q3 | p_T2 <= 0.05 | 0.35 | no (1.0000) | 0.1225 |
| Q4 | T0 fires in at least 30% of H0b worlds | 0.60 | yes (100%) | 0.1600 |
| Q5 | power at rho 0.3 >= 0.8 | 0.60 | yes (0.935) | 0.1600 |

Mean Brier score 0.1615. We expected a 55% chance that group-aware second-order structure would
survive pooling; it did not.

## Limits

- **Type level only.** The public table records which sequence types each group's fledglings
  produced (presence or absence), not how often each sequence was said. The preprint's claim is
  about token-level sequences, which are not public; it stays untested here. If the authors publish
  token data, the same tests will be registered on them first.
- **T0 is our stand-in.** PLAN.md was fixed before we could read the preprint's methods. The
  preprint compares in-sample likelihoods across smoothing values and permutes against a
  zero-order null; [AMENDMENT-1.md](AMENDMENT-1.md), registered before this job returned, re-builds
  both and is queued as job magpie-order-a1-v1. Its result will be added here.
- Group repertoires come from fledgling follows and may miss rare types. The loop definition is
  strict (at least two full copies of a chunk). The null generator redraws a sequence up to 200
  times to keep types distinct within a group; 12 and 8 power worlds (rho 0.3, 0.6) needed the
  fallback, none of the null worlds did.
- The power worlds plant one kind of second-order rule (shared by all groups). A rule that differs
  from group to group was not planted.

## Prior work

Mason, Walsh, King and Ridley, bioRxiv 2026, doi:10.64898/2026.06.15.732496 (the claim); Mason,
King and Ridley 2026, Proc R Soc B 293:20251620 (the data, figshare doi:10.6084/m9.figshare.28251398,
CC BY 4.0). Searches on 2026-11-24 (Google Scholar-style web search for the preprint title, for
"magpie call sequences second-order Markov" with critique, replication or pooling, and bioRxiv)
found no published re-analysis or calibration of this claim.
