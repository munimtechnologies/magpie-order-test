# magpie-order-test: pre-registered plan (2026-11-24)

MV2 of Munim, Inc. Made with AI agents working in his name; no human expert has checked this plan.
Text CC BY 4.0, code Apache-2.0.

## The claim

Mason, Walsh, King and Ridley (bioRxiv 2026, doi:10.64898/2026.06.15.732496, "Socially learned call
sequences reveal gradual development of higher-order structure") report that in Western Australian
magpie call sequences of up to 15 calls "each call depends on the two calls before it": second-order
Markov models explain the sequences better than first-order ones. They call this the first evidence
of multi-level non-adjacent organisation and learned syntactic structure in a semantic non-human
system. The same group's peer-reviewed paper (Mason, King and Ridley 2026, Proc R Soc B
293:20251620) reports that fledglings produce **the sequences of their own social group**, not of
the wider population.

## Why it might not mean syntax

Two simpler mechanisms produce second-order statistics without any second-order rule.

1. **Pooling.** If each group follows its own first-order transitions, a single pooled first-order
   chain is wrong for every group, and a second-order chain partly recovers which group a sequence
   comes from (the call before last is a clue to the group). A pooled order test can then prefer
   second order although no group uses a second-order rule.
2. **Looping.** Many long magpie sequences are a short chunk said again and again (for example a
   4-call chunk said almost four times gives the 15-call type). Inside a loop each call is fixed by
   the call one period back, which a second-order model rewards, but the "rule" is "repeat the
   chunk", not a dependence of each call on the two before.

## What has been seen before this plan (full disclosure)

From the authors' public CC BY 4.0 table (figshare doi:10.6084/m9.figshare.28251398 v2,
COMBINATION_repertoire.csv): the header line, and the 117 sequence-type labels (7 call units: DS,
LH, LHDS, NL, NLDS, SH, SHDS; lengths 2 to 15), as listed in the loader notes of MV2's earlier
magpie-dialect-test. Reading the labels is what suggested the looping mechanism above: several long
types are visibly repeated chunks, and no label we read had a unit directly repeated. No 0/1
column (fl01, group01, exgroup01) has been counted for this study, and no model has been fitted or
statistic computed. The preprint's methods text has not yet been read by us (bioRxiv blocks our
fetcher); the tests below are defined without it.

## Data used

- Types P: all distinct labels (117 expected).
- Group repertoires: for each natal group g, S_g = labels with fl01 = 1 or group01 = 1 on any row
  of a fledgling of group g (types the group's fledglings produced or that were recorded from the
  group). Records D = all (g, type) pairs; U = distinct types in D. Types in P but in no S_g are
  listed in the results.
- Each type is one sequence. This is a **type-level** test: the table has presence/absence, not
  token counts. If the authors publish token-level sequences, an amendment will add the same tests
  on tokens before those data are opened.

## Models and statistics

Sequences get a START context and an END outcome. Order-1 and order-2 chains.

- **T0, published-style order selection.** Maximum-likelihood order-1 and order-2 chains on U.
  T0 fires if order 2 has the lower AIC with k = observed non-zero outcomes minus one per seen
  context. AIC and BIC under this count and under the full count (7 free outcomes per seen context)
  are all reported, on U and on P. T0 firing on the real U = the claim reproduced at type level.
- **G, a cross-validated group-aware gain.** Leave-one-type-out over D (the held-out type is
  removed from its group and from every group in the pool), add-0.5 smoothing.
  Order-1 model: mu * P1(group) + (1 - mu) * P1(pool). Order-2 model: lambda * P2(group) +
  (1 - lambda) * (order-1 model). mu and lambda chosen on {0, .25, .5, .75, 1} to maximise the
  cross-validated log-likelihood (this favours order 2, which is conservative for our question).
  G = best order-2 minus best order-1 log-likelihood, in nats (G >= 0 by construction). G_pooled is
  the same without group terms and is reported.

## Null worlds (999 each; seed 20261124)

Every real type is replaced, group by group, by a sequence of the **same length**, distinct within
its group and without a directly repeated unit when the real type has none (up to 200 draws).
Transitions are each group's own maximum-likelihood unit-to-unit transitions (pooled, then uniform,
fallback for unseen rows); the first unit from the group's own start distribution.

- **H0a, group first-order.** The whole type is drawn from the group's first-order chain.
- **H0b, group first-order with loops.** For a type whose smallest period p satisfies len >= 2p
  (a chunk said at least twice), a chunk of length p is drawn from the group's first-order chain
  fitted on chunks, and looped to the type's length; other types are drawn whole.

## Tests

- **T1 (beyond pooling):** p = (1 + #{H0a worlds with G >= real G}) / 1000.
- **T2 (beyond pooling and loops):** the same against H0b.
- **T3 (calibration of the published-style test):** the share of H0a and of H0b worlds in which
  T0 fires. This goes into MV2's ledger of calibrated re-tests (procedure "magpie-markov-aic").
- **T4 (power):** 200 worlds each at rho = 0.3 and 0.6: H0b worlds in which a random share rho of
  the 49 second-order contexts get their own next-unit distribution (Dirichlet 0.5 over the units
  other than the previous one), shared by all groups. Power = share with G above the 95th
  percentile of G in the H0b worlds.

## Verdict rules (fixed now)

- p_T2 <= 0.05: **SECOND_ORDER_SURVIVES** (second-order structure beyond pooling and loops).
- else p_T1 <= 0.05: **LOOPING_EXPLAINS** if power at rho 0.3 >= 0.8, otherwise
  BEYOND_POOLING_LOOPING_UNDECIDED.
- else: **POOLING_EXPLAINS** if power at rho 0.3 >= 0.8, otherwise UNDECIDED_LOW_POWER.
- If T0 does not fire on the real U, the suffix _T0_NOT_REPRODUCED is added: the type-level data
  then do not show the published preference at all, and the authors' token-level result stays
  untested here.

## Predictions (priors, for a Brier score)

| id | prediction | prior |
|---|---|---|
| Q1 | T0 fires on the real U | 0.75 |
| Q2 | p_T1 <= 0.05 | 0.55 |
| Q3 | p_T2 <= 0.05 | 0.35 |
| Q4 | T0 fires in at least 30% of H0b worlds | 0.60 |
| Q5 | power at rho 0.3 >= 0.8 | 0.60 |

## Limits

Type-level only; no token frequencies, adults' own sequences or dates. Group repertoires come from
fledgling follows and may miss rare types. The loop definition is strict (at least two full copies
of the chunk), which leaves partial repeats to the order-2 model. The preprint's own procedure may
differ from T0; if its methods show that, an amendment will add it before the job's results are
read, and the original T0 stays reported.

## Code and run

magpie_order.py (numpy only), test_magpie_order.py (hand-traced unit tests), fetch_data.py (checks
the figshare MD5), gate.py (the job computes nothing until this plan and the code on GitHub are
identical to its own copies). One CPU job, 1 GB.
