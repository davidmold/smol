# Confirming the tail-check arms race: plan, written before the runs (22 Nov 2026)

The result being tested is `docs/findings.md` section 35: in gene-and-machine worlds whose machine
also checks the gene's tail tag (gmt), parasites from the future are increasingly harmful through
2M steps, an arms race that the same rules without the check (gmq) do not show.

This file is committed before any of the runs below exist, so the tests and thresholds cannot be
tuned to the data.

## Why it needs confirming

1. **The pooled statistic overstates the evidence.** Section 35's rho = -0.38 (t = -12.4) treats
   all 900 matrix cells as independent, but cells of one run share a world and each row shares a
   baseline. Recomputed with runs as the units (`confirm_stats.py`, output `confirm_existing.txt`):
   - gmt: per-run rho(lag, d) negative in 9/10 runs, median -0.60 (sign test p = 0.02, Wilcoxon
     p = 0.05). One run (s305) goes strongly the other way.
   - gmq: negative in 3/10, median +0.30. gmt vs gmq: Mann-Whitney p = 0.06.
   - gm (section 33, no QSEARCH): negative in 6/10, median -0.20, with two runs (s301, s307) as
     escalating as any gmt run. gmt vs gm: p = 0.5.
   So within gmt the sign is clear, but the difference from the controls is not yet established.
2. **A confound in the design.** MT is 6 symbols longer than the seed machine and slower (510
   steps per copy against 455). Any effect could come from a longer, slower machine rather than
   from excluding parasites.
3. **10 runs, one design, 2M steps.** Whether the race continues is untested.
4. **"Arms race" needs both sides to move.** Future parasites being more harmful also fits
   parasites improving against static hosts. The two-way decomposition in `confirm_stats.py`
   asks separately whether later parasites hurt any fixed host more (beta_P < 0) and whether
   later hosts do better against any fixed parasites (beta_H > 0).

## The control for the confound

MD = MT with its tail test `?OMSJ}` replaced by a second copy of the head test `?WGEW}`. Checked
before the runs (`gmd_check.py`, output `gmd_check.txt`), against the gene-parasites of the
section 33 runs (50 windows, QSEARCH on):

| machine | symbols | loop closes | steps per copy | share of bound parasites copied (median over runs) |
|---|---|---|---|---|
| M0 (seed) | 42 | yes | 455 | 0.92 |
| MT (tail check) | 48 | yes | 510 | 0.08 |
| MD (dummy check) | 48 | yes | 510 | 0.89 |

MD matches MT in length and speed and M0 in what it lets through.

## Runs

- cf_gmt, cf_gmq, cf_gmd: 30 new seeds each (311-340), section 33 layout, QSEARCH 1, 2M steps;
  time-shift matrix as section 35 (times 0.2-2M, 2 seeds per cell, gene-parasites swapped).
- gmt6, gmq6: the original 10 gmt and 10 gmq runs continued from 2M to 6M; time-shift matrix on
  2-6M every 0.5M.
- Mechanism probe (`gmt_mech.py`) on the three new sets.

Script: `confirm_pipeline.sh` (setup `confirm_setup.py`). Estimated about two days on aida.

## Tests, decided now

Unit: the run. Statistic: rho_run = Spearman rank correlation of d with the lag over that run's
off-diagonal cells. Runs that go extinct before 2M are reported and excluded from rho_run.

**Primary (all three must hold to call the result confirmed):**
- P1. cf_gmt escalates: rho_run < 0 in most runs, Wilcoxon signed-rank two-sided p < 0.01.
- P2. cf_gmt escalates more than cf_gmq: Mann-Whitney on rho_run, p < 0.025.
- P3. cf_gmt escalates more than cf_gmd: Mann-Whitney on rho_run, p < 0.025.
  (0.025 = 0.05 split over the two control comparisons.)

**Secondary:**
- S1. Both sides move: in cf_gmt, beta_P < 0 and beta_H > 0 each with Wilcoxon p < 0.05, and the
  share of runs with both is higher than in cf_gmq and cf_gmd.
- S2. The race continues: in gmt6 (2-6M), rho_run < 0 in most runs (with 10 runs, sign test
  p < 0.05 needs 9/10); gmq6 shown alongside.
- S3. Mechanism replicates: in cf_gmt, late machines copy early parasites less than early machines
  do, and late parasites carry the tail tag (share >= 0.9).

**Readings fixed in advance:**
- P1-P3 hold: the arms race is confirmed and is caused by the recognition step, not by the
  machine's length or speed.
- P1 holds, P3 fails: escalation comes with the longer machine, not with recognition.
- P1 fails: section 35 was a small-sample result; it is withdrawn in the doc.
- S2 fails with P1-P3 holding: a real but transient race (it settles after 2M, like the turnover
  in the published runs).

Not tested here: whether a recognition step like this can evolve rather than being designed in.
