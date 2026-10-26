# Next steps (26 Oct, afternoon)

Done and written up in findings 14: all 38 small runs to 10M (21 reached it, 17 extinct; no
processive species established, three brief sightings below 1%), a random-sample scan that
closes the rare-variant blind spot, and the memory fix that let the killed runs finish.

## Running (on aida)

- Big classic replicates `runs/big/classic_rep{19019,20020,21021}` to 10M, all now on the
  hashed-species binary with `SPECIESFORGET` (at 4.4M, 5.2M, 4.0M on 26 Oct 14:30; ~23
  steps/s, so ~3 days). No processive species yet. Scan each every 500k with both
  `transition_scan.py --min 20` and `--sample 1000`.

## Then

1. Follow pr0xx seed21021 past its 9.98M detection (running to 13M, checkpoints every 50k):
   scan finely (`--sample 1000`, every 10-20k) and see whether the processive variants spread,
   persist at low level or vanish.
2. Introduction size (`experiments/processive/sizes.py`, running from 26 Oct; 270 runs of
   200k steps into communities seed36036 and seed21021 at 2.5M; sequences E1, S66, S21, S44):
   1, 2 adjacent, 5 and 20 molecules, each with a matched cleared-space control. A lone processive molecule cannot react with a copy of itself until one appears;
   the adjacent pair separates encounter limitation from establishment failure. Use the
   small-world sightings' own sequences alongside the big-world ones: if they fail even in
   groups, genotype or community fit matters, not just starting abundance.
3. Transition frequency: small worlds none established in 38 (fresh and survivor cohorts kept
   separate; 21 reached 10M), big worlds 1 of 1 so far with 3 replicates pending. Compare at
   matched activity in attempts. Three negative replicates would mean "not reproduced in three
   runs", not "rare".
4. Big replicates: does a parasitic crisis precede any transition? seed20020 has a different,
   low-novelty history.
5. Why seed21021 resists processive invasion.
6. Do partnerships accumulate further organisation?
7. Measurement: observed-reaction percentages; a drift null with space.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- Memory: watch `free -g` while more than ~6 runs share aida; new long runs should set
  `SPECIESFORGET 1` unless species-id continuity matters.
- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
