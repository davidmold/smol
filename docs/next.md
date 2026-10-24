# Next steps (24 Oct, afternoon)

Nothing is running on aida. Done since the last update and written up in findings:
descendant tracking and removal tests on evolved partnerships (13), the big runs to ~10M
including the classic world's second transition to processive replicators (10b), and the
novelty recheck with a placed-birth null, which withdraws the "selection suppresses novelty"
claim (1).

## Running (24 Oct evening)

- Three big classic replicates, `runs/big/classic_rep{19019,20020,21021}`, to 10M with
  OFFSPRINGLOG, GZIPLOGS and a checkpoint kept every 100k (2-4 days).
- `experiments/processive/long_small.sh`: the 18 small survivors continued 2.5M -> 10M
  (`runs/pr0xx`, survivor-conditioned) and 20 fresh small runs 0 -> 10M (`runs/long`), 9 at a time
  (about 3 days). Score both with `transition_scan.py` (detection, establishment, by steps and
  by cumulative placed copies; extinction as a competing outcome).

## Then

1. Transition frequency: small vs big worlds over matched steps and matched reproductive
   activity; does the parasitic crisis precede emergence and destroy small worlds?
2. Why seed21021 resists processive invasion.
3. Do partnerships accumulate further organisation?
4. Measurement: observed-reaction percentages; a drift null with space.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
