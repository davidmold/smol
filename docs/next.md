# Next steps (25 Oct, evening)

Done since the last update and written up in findings 14: the processive-species definition
and screen, the transplants, the first long small-world results (no processive species in 38
runs; extinction common), and a memory fix after five runs were killed for lack of memory.

## Running (25 Oct evening, on aida)

- Big classic replicates `runs/big/classic_rep{19019,20020,21021}` to 10M. 19019 and 21021
  continue from checkpoints (2.2M, 1.7M) with `SPECIESFORGET` on; 20020 (2.9M) is on the old
  binary and stays under 2 GB, so it is left alone. At ~23 steps/s each, 10M is ~4 days away.
- Small worlds to 10M: long seed64064 and seed68068 (old binary, finishing), seed66066
  (resumed at 6.1M), and the four OOM-killed survivors pr0xx 21021, 24024, 29029, 44044
  (resumed by `/tmp/queue3.sh` as memory frees). Each is scanned when it finishes;
  `~/published/scan_long/DONE2` marks the end.

## Then

1. Transition frequency, small vs big, matched by steps **and** by cumulative copying. The small
   runs to 10M have had a third to two thirds of the copying (attempts) the big world had at
   first detection; compare at matched activity, in attempts, since the original big run did
   not log placed copies.
2. Big replicates: scan every 500k past 3M; does a parasitic crisis precede any transition?
   seed20020 already has a very different, low-novelty history.
3. Why seed21021 resists processive invasion.
4. Do partnerships accumulate further organisation?
5. Measurement: observed-reaction percentages; a drift null with space.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- Memory: watch `free -g` while more than ~6 runs share aida; new long runs should set
  `SPECIESFORGET 1` unless species-id continuity matters.
- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
