# Next steps (29 Oct, evening)

Done and written up in findings 14: small worlds to 10M, transplants, introduction size,
introductions without clearing (traced), replays of seed21021, W26 vs S21, and the three big
replicates to 10M (none made the transition). Three conclusions set the direction: a small
world can make the processive transition spontaneously (1 of 24 clean replays from seed21021
at 9.9M); "permissive versus resistant worlds" is retired (a cleared-patch artefact); and more
evolutionary activity alone does not make the transition routine (big replicates).

## Done 29 Oct: disturbance (findings 15)

Clustered and dispersed disturbance (one modest dose, chosen in advance) did not change how often
processive variants appeared or established from seed21021 at 9.9M (base 2/40, clustered 1/40,
dispersed 0/40); hypothesis closed. Baseline now 3 of 64 clean replays (~5%), all winners short
S21-family replicators that persisted to 400k.

## Proposed test of the findings-16 prediction (not started)

Neighbourhood x genotype, matched: from seed21021 at 9.9M, place an adjacent pair of W26 or of
E1, replacing residents (no clearing), at 40 'family' sites (>= 70% of the 10 neighbouring cells
hold a `>C$==` carrier) and 40 'non-family' sites (< 30%), chosen at random from the ~3,400 and
~1,900 such sites; each site gets W26 and E1 with the same random stream. 160 runs of 200k,
traced; ~5 h on aida plus scoring. Readouts: establishment (>= 10%, introduced family) and,
mechanistically, the share of residents' copies of the newcomer that are new members (first
5k steps). Prediction: W26 family > W26 non-family (from the data so far ~30% vs ~10%), with
E1 showing no such difference or the reverse (the interaction is the test). If W26 does no
better at family sites, the copier-compatibility account is wrong.

## Deferred

- More big replicates (frequency at scale; ~4 days per 3).
- S21 dissection (halves; length vs copy blocks vs susceptibility to residents). Padding with
  "inert" symbols would not isolate length: every symbol affects binding and search labels.
- Replays of the original big run from its 2.5-3.0M checkpoints (in the S3 archive).
- Partnerships accumulating further organisation; measurement (observed-reaction percentages,
  a drift null with space).

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- Memory: watch `free -g` while more than ~6 runs share aida; new long runs should set
  `SPECIESFORGET 1` unless species-id continuity matters.
- The `smol` repo is on GitHub as davidmold/smol (private, 15 Nov), with a README for outside
  readers (Simon) and the GPL-3.0 licence text.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
