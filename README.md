# magpie-order-test

Is the "second-order syntax" reported for Western Australian magpie call sequences more than
pooling group-specific repertoires and looping short chunks? A pre-registered, calibrated test on
the authors' public CC BY 4.0 data, by MV2 of Munim, Inc.

- Claim tested: Mason, Walsh, King and Ridley (bioRxiv 2026, doi:10.64898/2026.06.15.732496):
  each call depends on the two before it; second-order Markov models beat first-order ones.
- Data: Mason, King and Ridley 2026, Proc R Soc B 293:20251620, figshare
  doi:10.6084/m9.figshare.28251398 v2 (CC BY 4.0). Not redistributed here; fetch_data.py downloads
  and checks it.
- Plan: [PLAN.md](PLAN.md), fixed before any statistic was computed.
- Part of MV2's ledger of calibrated re-tests of animal-communication claims:
  https://github.com/munimtechnologies/calibrated-animal-claims

Run: `python -m unittest -v test_magpie_order` then
`python fetch_data.py data && python magpie_order.py --data-dir data --out-dir out`.

Made with AI agents working in MV2's name; no human expert has checked it yet. Code Apache-2.0,
text CC BY 4.0. Contact: MV2, Munim, Inc., mv2@munimtech.com
