# Next steps (29 Oct)

Done and written up in findings 14: small worlds to 10M, transplants, introduction size,
introductions without clearing (traced), replays of seed21021, W26 vs S21, and the three big
replicates to 10M (none made the transition). Three conclusions set the direction: a small
world can make the processive transition spontaneously (1 of 24 clean replays from seed21021
at 9.9M); "permissive versus resistant worlds" is retired (a cleared-patch artefact); and more
evolutionary activity alone does not make the transition routine (big replicates).

## Now: does spatially clustered disturbance help innovations establish? (running from 29 Oct)

Hypothesis: occasional local openings let new reproductive strategies expand alongside
compatible neighbours before resident interactions eliminate them (motivated by the traces:
residents often alter newcomers, and many products find no room). It could equally destroy
useful associations or favour ordinary fast copiers.

Design (`experiments/disturb/disturb.py`, scored by `disturb_agg.py`): from seed21021 at 9.9M, 40 reseeded runs per treatment, chemistry
unchanged, nothing introduced, disturbance blind to genotype and to any detected novelty:
- baseline: existing dynamics;
- clustered: every 1,000 steps remove everything in one random 7x7 patch (49 cells);
- dispersed: every 1,000 steps remove whatever is in 49 random cells across the world (same
  expected mortality and freed space, no spatial clustering).
Disturbance runs for the first 150k steps, then stops; all runs go to 200k. One preselected
schedule, no parameter sweep. Measures: appearance of processive variants (`--sample 1000`
every 10k); establishment sustained over >= 3 consecutive samples; persistence after
disturbance stops; extinction, occupancy and placed copies, so that turnover is not mistaken
for evolutionary potential. Established populations are continued to 400k to tell a bloom from
a lasting change, and histories are kept.

Decision rule: clustered beats dispersed -> spatial opportunity is a supported mechanism; test
it in another starting community and, later, with another kind of innovation. Both similar ->
available space or turnover is the lead. Neither substantial -> close this hypothesis (no
sweep). This is a screen for a large effect (baseline ~1/24), not an estimate of a small one.

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
