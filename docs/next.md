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

1. Transition frequency: small worlds 0 of 38 (with brief sightings), big worlds 1 of 1 so far,
   3 replicates pending. Compare at matched activity in attempts (the original big run did not
   log placed copies). If none of the replicates transitions, the original was a rare event.
2. Foothold vs persistence: the sightings say processive variants do arise in small worlds and
   vanish. Introduce a single processive molecule (not 20) into evolved small communities to
   measure the establishment probability from one copy.
3. Big replicates: does a parasitic crisis precede any transition? seed20020 has a different,
   low-novelty history.
4. Why seed21021 resists processive invasion.
5. Do partnerships accumulate further organisation?
6. Measurement: observed-reaction percentages; a drift null with space.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- Memory: watch `free -g` while more than ~6 runs share aida; new long runs should set
  `SPECIESFORGET 1` unless species-id continuity matters.
- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
