# Experiments

Every run lives on aida under `~/smol/runs/<condition>/<run>/` with its `run.conf`, `log.txt`
(or `nohup.out`), `species.tsv`, `popdy001.dat`, `copies.tsv`, `grid.snap` and checkpoints.
Short descriptions of each condition are also in `~/smol/runs/conditions.json`, which the
dashboard shows. Finished conditions have their text outputs gzipped (see tools.md).

**Archived to S3.** To free aida's disk, finished run directories were moved to
`s3://lantic-smol-archive/runs/<condition>/<run>.tar` (AWS profile `lantic`, us-east-1,
Glacier Instant Retrieval: readable at once, ~$4/TB-month; 85 objects, 217 GB on 27 Oct):
- 26 Oct: the three original big runs (`big/classic_400x300`, `gradient_400x300`,
  `strand_400x300`), all of `long/`, `pr0xx/` except seed21021, and `transplant/` (one tar per
  recipient).
- 27 Oct: `sizes`, `pr0xx/seed21021` (continued to 13M), and the finished conditions of 16-21 Oct:
  b2, bindT3, cd, cell, cellhi, cp, cross, dep, disp0.01, disp0.03, disp0.1, disp0.3, disp1.0,
  dmg, evo21, evo57, gm, gmc, half, new, old, partners, pr, pr0, pr0log, ps, pub, replay,
  strand0, strandTS, vert.
- 27 Oct: `natural/<recipient>` archived in full; on aida each run directory keeps only
  `trace.tsv`, `scan.json`, `intro.txt` and `log.txt` (`natural/_src` is kept: `replays.py`
  reads it).
Each tar's SHA-256 is in `s3://lantic-smol-archive/manifest/` and in
`experiments/archive_manifest.tsv`; it was checked against the S3 object size before the aida
copy was deleted. The scan and analysis outputs derived from them stay in this repo and in
`~/published` on aida. Restore one run:
`aws --profile lantic s3 cp s3://lantic-smol-archive/runs/long/seed59059.tar - | tar -xf -`

**Deleted, not archived (27 Oct, at David's request):** the original C++ Stringmol runs in
`aida:~/stringmol/runs` (14-17 Oct; ~172 GB, mostly a full snapshot every 100 steps), superseded
by smol. Their names are in `experiments/stringmol_runs_deleted.txt`; `~/stringmol/runs/batch.log`
and the C++ source were kept.

(`tools/archive_s3.sh` does the uploads, from chasca, so no keys live on aida).

Seeds: 40-seed batches use RANDSEED 18018, 19019, ..., 57057 (run dirs `seed18018` ...);
17-seed batches use the first 17 of those. "Reached 1M" counts runs whose log shows step
1,000,000.

## Rule checks on the classic chemistry (Oct 16-17)

| Condition | Rules | Reached 1M |
|---|---|---|
| `pr0` | **Baseline**: random partner, pair decay, failed encounters cost nothing | 24 / 40 |
| `old` | Old code's rules: scan-order partner, failed encounter uses partner's turn | 11 / 17 |
| `ps` | Old rules, 40 seeds | 19 / 40 |
| `pr` | Random partner, pair decay, old turn rule | 23 / 40 |
| `cd` | Random partner, pair decay, old turn rule (17 seeds) | 7 / 17 |
| `cp` | Scan-order partner, independent decay | 1 / 17 |
| `new` | Random partner, independent decay | 1 / 17 |
| `half` | Independent decay at half the decay rate | 3 / 17 |

Seed: the original 65-symbol replicase
`WWGEWLHHHRLUEUWJJJRJXUUUDYGRHJLRWWRE$BLUBO^B>C$=?>$$BLUBO%}OYHOB`, 121 copies in an 11x11
block, 125x100 grid, DECAY 0.0005, MUTATE 0.0002.

## Chemistry changes on the classic chemistry (Oct 18)

| Condition | Change | Reached 1M |
|---|---|---|
| `disp0.01` ... `disp1.0` | Strand displacement at rate 0.01 / 0.03 / 0.1 / 0.3 / 1.0 | 16, 11, 10, 1, 0 of 40 |
| (restore tests) | `DISPLACE_RESTORE 1` at rates 0.1 / 0.3 / 1.0, 20 seeds, to 300k | 8, 8, 3 of 20 at 300k |
| `bindT3` | `BIND_T 3` absolute-score binding (run on chasca, results copied) | 17 / 40 |
| `evo57` | Seeded with evolved code-first replicator from `pr0/seed57057` (26 symbols) | 34 / 40 |
| `evo21` | Seeded with evolved letter-prefix replicator from `pr0/seed21021` (28 symbols) | 24 / 40 |
| `cross` | Designed cross-replicating pair, `designs/cross.txt` | 40 / 40 |

## Copy-strand chemistry (Oct 18-19)

Seed `RERERERE^BQQ$=?>%}` unless stated (`designs/strand.conf`).

| Condition | Change | Reached 1M |
|---|---|---|
| `strand0` | `STRANDCOPY 1`: **baseline for this chemistry** (run on chasca) | 27 / 40 |
| `strandTS` | `STRANDCOPY 1` + `TWOSIDED 1` | 40 / 40 |
| `gm` | + `TRANSLATE 1`, seeded with gene G + machine M (`designs/gm.*`) | 40 / 40 (G+M persists) |
| `gmc` | as `gm` plus 25 plain copiers 36 cells away | 20 / 40 alive, G+M in 0 / 40 |
| `dmg/{cop,gm,gmc}_<rate>` | `DAMAGE` at 1e-6, 3e-6, 1e-5, 3e-5 | see findings 9 |
| `cell/`, `cellhi/` | cell-cycle machine (`designs/cell*.`), with and without copiers, damage up to 1e-3, `DAMAGE_SHIELD_BOUND` | cell loses to copiers at every rate |

## Big worlds (Oct 20-23; stopped at ~10M steps)

400x300 grid (120,000 cells), RANDSEED 18018, 121-molecule seed block in the centre,
NSTEPS 100M (stopped by hand at ~10M: classic 9.8M, strand 10.0M, gradient 12.9M), LIVEEVERY 200, CKPTEVERY 100k (5 kept). Live: `/live?run=big/<name>`.

| Run | Chemistry | Notes |
|---|---|---|
| `big/strand_400x300` | copy strand, uniform | partnerships dominate |
| `big/gradient_400x300` | copy strand, `RATE_MIN 0.1` at x=0 to `RATE_MAX 1.0` at x=200 | fewer pairs, less diversity |
| `big/classic_400x300` | classic, original replicase, `ENCOUNTER_TURN 0` | same families as small runs, 10x as many |

Analysis outputs: `~/smol/bigscore/` (`real.txt`, `nullscore.txt`, `blocks.txt`,
`blocks2.txt`; drift null in `bigscore/null/`).

## Side-on worlds (Oct 20)

All `VERTICAL 1`, `STRANDCOPY 1`, `DECAY_CHEM 10` (`DECAY_CHEM_SEED 1`), seed planted on the
ground row. Weights are written to `decay_weights.tsv` in each run.

| Run | Shape and light | Outcome |
|---|---|---|
| `vert/strand_120x60` | 120 wide x 60 tall, light by height 0.2 -> 1.0 | turf ~4 high, died at 547k |
| `vert/strand_400x12` | 400 x 12, light by height 0.2 -> 1.0 | running (20M); fragile winners |
| `vert/shade_400x12` | 400 x 12, `SHADE 0.5`, no height field | died at 379k |
| `vert/shade_1200x12` | as above, 1200 wide | started by mistake, stopped at ~3.3M |
| `vert/slab_400x12x5` | 400 x 12 x 5 slices (`DEPTH 5`), `SHADE 0.5`, seed in middle slice | running (20M) |

## Partnership isolation and follow-ups (Oct 21-24)

| Where (aida `~/smol/runs/`) | What | Results in repo |
|---|---|---|
| `partners/` | Batch 1: four pairs alone, vs copier, with parasite; mutation off/on; 10 seeds; 200k | `experiments/partners/agg1.txt`, `lineage_agg.txt` |
| `b2/` | Batch 2: invasion from rarity, reverse invasion, bottlenecks, rescue, accessibility | `experiments/partners/agg2.txt` |
| `dep/` | Removal tests on evolved descendants (from batch-1 mutation-on runs at 200k) | `experiments/partners/dep_agg.txt` |
| `replay/` | strand0 seeds rerun with a checkpoint every 10k (identical to the originals) | used by batch 2 |

## Published comparison and novelty (Oct 20-24)

| Where | What | Results in repo |
|---|---|---|
| `~/published/bio/` | The authors' archive (8 runs to 2M) | `experiments/published/` |
| `pub/` | Published configs with old rules, 5 seeds each, 2M | `experiments/published/assay_pubrerun.txt` |
| `pr0x/` | 24 baseline survivors continued 1M to 2.5M | `experiments/published/assay_pr0x.txt` |
| `pr0log/` | The same 24 baseline runs rerun to 1M with `OFFSPRINGLOG` | `experiments/novelty/` |
| `~/smol/nulls_pr0placed/` | Drift null calibrated from placed births | `experiments/novelty/wn_null_placed.txt` |

Big-run analyses to 10M: `experiments/big/` (copy mixes, lengths, assays, classic transition).

## Processive replicators (Oct 24-)

| Where (aida `~/smol/runs/`) | What | Results in repo |
|---|---|---|
| (scan only) | `transition_scan.py` over big classic and 87 existing small/published runs | `experiments/processive/screen/` |
| `transplant/` | E1, E2, M1, M2, ORD into an empty world and six 2.5M communities; 3 reps; 300k | `experiments/processive/tp_agg.txt`, `scan_tp/` |
| `big/classic_rep*` | three big classic replicates to 10M (running) | |
| `pr0xx/`, `long/` | small worlds to 10M: 18 continued survivors, 20 fresh (running) | |

## Paper runs continued in smol (Nov 15)

The 8 published runs at 2M, converted with `tools/pub_to_ckpt.py` (old partner and turn rules),
continued to 3M with `--reseed N` on chasca, 8 in parallel. Runs are in `/tmp/sm_time/smolpub/`
on chasca (not archived).
- Wall time 2M -> 3M: 35-44 min per run (C++ dev: 2.6 s/step, about 30 days per 1M).
- All 8 survived. Molecules at 3M: 7,920 / 9,516 / 10,957 / 11,788 / 5,154 / 9,757 / 11,195 /
  7,693 (runs 1-8).
- In every run, the most abundant species at 2M had gone by 3M and something else led.
- These runs lacked `COMPAT_DELETION` (found the same day, section 17). A rerun with it is in
  `/tmp/sm_time/smolpub2/`.
