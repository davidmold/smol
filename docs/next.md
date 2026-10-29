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

## Open: what next (not yet decided)

The same question, a genotype-blind ecological change that raises establishment, can be asked
with a different lever (one at a time, each with a decision rule set in advance), e.g. lower
overall density (more empty space everywhere) or a gradient in decay. Or step back to the
deferred list below.

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
- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
