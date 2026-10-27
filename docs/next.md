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

1. Native transitions in seed21021: from 9.9M 3 of 56 runs established a short (26-32 symbol)
   S21-family replicator; from 10.0M, with the original 95-symbol lineage present, 0 of 88. Is
   the short single-block form what makes the difference? Introduce the 26-symbol winner
   (`Q>BO^$O%H>C$===?O>A$$T$TQ>`) and S21 into the 9.9M and 10.0M states, singly and as pairs
   (as in natural.py), and trace them; and replay from 9.95M/9.97M (needs checkpoints: rerun
   9.9M -> 10.0M once with CKPTEVERY 10000) to see how often the family arises at all.
2. Why S21 fails in its own late community (0/64) but not in seed36036 (9/16 as pairs): its
   founders there are mostly used as templates by residents and often altered. Compare the
   residents that bind newcomers at 2.5M and at 9.9M (who they are, what their code does to
   the newcomer), from the traces.
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
