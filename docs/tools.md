# Machines, code and tools

## Where the code is

| What | Where | Notes |
|---|---|---|
| Rust simulator `smol` (source of truth) | chasca `~/projects/smol` | git, **local only, no remote** |
| Deployed copy + all runs | aida `~/smol` (runs in `~/smol/runs`) | not a git checkout; updated by tar (below) |
| C++ fork of Stringmol + analysis scripts | aida `~/stringmol` | github.com/davidmold/stringmol, branch `fast` |
| Upstream contributions | chasca `~/projects/stringmol-uoy` | fork of uoy-research/stringmol |
| Python env for analysis | aida `~/smenv` (`~/smenv/bin/python`) | rapidfuzz, numpy |
| Design files for seeds | `designs/` in the repo | cross pair, G+M, cell, evolved seeds, strand seed |

Commits use `git -c user.name=davidmold -c user.email=davidmold@gmail.com`.

## Deploying a change to aida

```
cd ~/projects/smol && tar czf /tmp/smolsrc.tgz src Cargo.toml && scp /tmp/smolsrc.tgz aida:/tmp/
ssh aida 'cd ~/smol && tar xzf /tmp/smolsrc.tgz && touch src/*.rs src/bin/*.rs && ~/.cargo/bin/cargo build --release'
```

`cargo` is not on PATH in non-interactive ssh. Running simulations keep their old binary;
new switches default off, so rebuilding mid-batch is safe.

Before committing any chemistry change: build the previous commit in a worktree
(`git worktree add /tmp/smolhead HEAD`), run both binaries on a classic config, a gradient
config and the latest special config, and `cmp` `popdy001.dat`, `species.tsv`, `copies.tsv`.

## Binaries (`target/release/`)

| Command | Purpose |
|---|---|
| `smol run CONF` | Run in the current directory |
| `smol resume CKPT [--reseed S] [--nsteps N]` | Continue exactly from a checkpoint; `--reseed` for replay experiments |
| `smol-react PAIRS [MAXSTEPS] [--compat] [--skip N]` plus `--twosided --trace --strand --translate [seed] --stop --bind-t T` | Reaction probe: react listed molecule pairs, report products and steps |
| `smol-pairprob [BINDLEN] [BIND_T]` | Binding probabilities between molecules |
| `smol-codon encode STRING [SEED]` / `decode GENE [SEED]` / `table [SEED]` | Translate between machine strings and letter genes |

Launching a long run on aida so it survives the ssh session:
`cd RUNDIR && (setsid nohup ~/smol/target/release/smol run run.conf > nohup.out 2>&1 < /dev/null &)`.
Don't use `pkill -f` with a pattern that matches your own ssh command line. To stop one run,
find its pid via `readlink /proc/<pid>/cwd`.

## Outputs of a run

`popdy001.dat` (step,species,count), `species.tsv` (id, birth step, parent active, parent
passive, sequence), `copies.tsv` (step, active, passive, product, count, side), `grid.snap`
(zlib frames of species ids, bound flag in bit 31), `checkpoint_N.txt.gz`, `log.txt` or
`nohup.out` (one line per 1000 steps with event counters), and with `LIVEEVERY` `live.grid` +
`live.json`. Vertical runs with `DECAY_CHEM` also write `decay_weights.tsv`.

## Analysis scripts

| Script | Where | What |
|---|---|---|
| `window_novelty.py RUNS_DIR [--glob seed*] [--window 100000] [--thresh 0.25]` | aida `~/stringmol` | The main novelty metric: new sequence families per window (species >= 1% and 20 molecules) |
| `neutral_null.py RUN OUT [--start 300000]` | aida `~/stringmol` | Drift null: Moran process with the run's own N, turnover and new-species rate; writes files the metric reads. Existing nulls in `~/smol/nulls_pr0`, `nulls_strand0`, `nulls_strandTS` |
| `families.py`, `summarise.py` | aida `~/stringmol` | Family clustering per run; summary across runs |
| `probe_families.py pairs/analyse` | aida `~/stringmol` | Are families functional? (builds probe pair lists, analyses `smol-react` output) |
| `tools/partners.py RUN...` | repo; aida `~/smol/partners.py` | Borrower partnerships: first seen, established, ancestry, share |
| `tools/replicators.py RUN... [--window] [--top]` | repo; aida `~/smol/replicators.py` | Dominant replicators across runs |
| `evo_report.sh`, `strand_report.sh` | aida `~/smol` | Batch reports used for the Oct 18 results |
| `viz/build.py RUN OUT.html EVERY` | repo / aida `~/smol/viz` | Animated grid viewer for one run |

Quick scripts from 20 Oct, saved in `tools/` (paths inside are hard-coded to aida; edit before
reuse):

| Script | What |
|---|---|
| `tools/who_copies_whom.py RUN_DIR [WINDOW]` | Last-window copy relations: self vs cross, copier vs code-less, mutual pairs, hosts |
| `tools/block_novelty.py` | Per-100x100-block novelty from checkpoints for `big/classic_400x300`, plus small runs scored the same way |
| `tools/big_summary.py` | Species, effective species, mean length over time and by band, for the three big runs |
| `tools/vertical_profile.py RUN_DIR` | Occupancy and top species by height in a side-on run |
| `tools/vertical_strata.py RUN_DIR...` | Decay weight and share of copiers by depth from the column top (handles DEPTH slabs) |

## Dashboard and live viewer

`web/server.py` runs on aida as the user service `smol-dash`
(`systemctl --user restart smol-dash` after changing it). Port 8095:
`/` lists conditions and runs with their descriptions from `runs/conditions.json`;
`/live?run=<condition>/<run>` shows the live grid (colour by species, brighter when bound,
click to highlight, species list with roles, who-copies-whom, population chart; refreshes
every 3 s). Viewer HTML is cached in `~/smol/viewcache`.

## Data and disk

- The finished small-run batches had their `popdy001.dat`, `species.tsv` and `copies.tsv`
  gzipped on 20 Oct (110 GB -> 29 GB). All the analysis scripts above read either form
  (`ropen`/`rexists` helpers in each script; pre-change copies kept on aida as `.bak`).
- The big runs write roughly 0.2-0.8 GB per million steps each; 100M steps fits in the free
  space (~239 GB) but `copies.tsv` is the part to segment or trim if needed.
- Checkpoints are what replays and per-block scoring need; consider keeping more, spaced out.

## Upstream (Simon Hickinbotham's Stringmol)

- PR #50 (alignment-score fix + regression test) and PR #51 (RNGFILE fix: restarts restore the
  random generator, + restart test) are merged into uoy-research/stringmol `dev` (15 and 19 Oct;
  Simon's note on #51, 15 Nov: apologies for the delay).
- PR #52 (merged 15 Nov): `perf-agentsalign`, rebased onto the new `dev` (dropping `cf715fc`):
  removes the second, discarded alignment per binding attempt and the `swlist` cache that never
  stored anything. Regression, restart and Catch tests pass; popdy identical to `dev` for 20,000
  steps. Time saved: 12% over steps 0-10k, 38% over 10k-20k (v2.0.3 config, seed 4357, two runs
  each on chasca).
- Next: `perf-small` (BLOSUM row-lookup cache, species hash index, single-pass species counts),
  rebased onto `dev` 109be57 as local branch `perf-small-rebase` (3 commits, cf715fc dropped);
  regression, restart and Catch tests pass. Opened as PR #53 (15 Nov): 6% saved over steps 0-10k,
  3% over 10k-20k, popdy identical. The PR also explains why the agent list needs rewriting next.
- Not yet ported: replacing the agent list (95% of the old code's run time).
- Restarting the paper's runs (RSOS 8:210441; data on aida ~/published/bio/biopaperdata/raw/run1-8):
  every run has out1_STEP.conf + RNGstate_STEP.txt every 20k steps to 2,000,000, so `dev` (with
  #51) restarts from them and reloads the archived generator state. Speed from run2 at 2M
  (9,078 molecules, 125x100), chasca, 40 steps each: `dev` 2.6 s/step, #52 2.6, perf-small 2.6,
  identical output. ~33k steps/day, ~30 days per extra million steps per run. #52 and perf-small
  don't help at this population; the agent-list rewrite is what would.

- `tools/pub_to_ckpt.py OUT1_CONF OUT_DIR [--nsteps N] [--set K V] [--new-rules]` turns a C++
  snapshot (e.g. the paper's `out1_2000000.conf`) into a smol checkpoint + species.tsv: every free
  molecule and every bound pair mid-reaction, with both strings and all four pointers (C++ `irwf`
  maps one-to-one onto smol's Pointers). Resume with `--reseed S` (the C++ RNG can't carry over).
  Defaults to the old partner/turn rules, as the paper's v0.2.3.4 used.

## Interventions on running worlds

`tools/ckpt_edit.py RUN_DIR CKPT OUT_DIR [--remove SEQ] [--remove-species ID] [--cull F]
[--region X0,X1,Y0,Y1] [--add SEQ N] [--set KEY VALUE] [--rng-seed N]` writes an edited
checkpoint plus `species.tsv` to OUT_DIR; then `cd OUT_DIR && smol resume CKPT [--reseed S]`.
A no-op edit resumes byte-identically to the original. `--put X,Y SEQ` replaces whatever is in that cell with one
free molecule of SEQ (no clearing around it). `--live-species` copies only the species
present after the edit (pair with `--set SPECIESFORGET 1`), instead of the recipient's whole list. Used for rarity invasion, bottlenecks,
removal and rescue, and for continuing runs with new settings (e.g. `--set NSTEPS`).

## Scanning for processive species

`tools/transition_scan.py RUN [--every N] [--start S] [--min M] [--cap C] [--sample K] [--cache F] [--json F]`:
probes species in a self-reaction and reports processive share, detection and establishment and
cumulative copies. By default it probes the most abundant species (>= M molecules, up to C);
`--sample K` probes the species of K random molecules instead, an unbiased share estimate when
most molecules are rare variants (column `cover` then means the fraction of the sample probed).

## Archiving to S3

`tools/archive_s3.sh [--delete] RUNPATH...` (run on chasca): streams each `aida:~/smol/runs/RUNPATH`
as a tar to `s3://lantic-smol-archive/runs/RUNPATH.tar` (Glacier Instant Retrieval), records the
SHA-256 and size, checks the object's size in S3, and with `--delete` then removes the aida copy.
It skips directories a running `smol` or scan is using, and objects already uploaded.

## Stopping and resuming long runs

- `tools/trim_logs.py RUN_DIR STEP [--dry]`: before resuming a GZIPLOGS run from
  `checkpoint_STEP`, drop the logged data after STEP (whole trailing gzip members, including one
  cut off by a crash; the tail is kept as `*.trimmed_tail`), so steps are not logged twice.
  Trim + resume reproduces an uninterrupted run's popdy, species and copies logs exactly.
- `tools/ckpt_setkey.py IN OUT KEY VALUE`: set a config key in a checkpoint (e.g.
  `SPECIESFORGET 1`); resume from OUT in the same directory.
- `tools/fold_species.sh RUN_DIR`: fold a plain `species.tsv` (copied in by `ckpt_edit.py`) into
  `species.tsv.gz` as its first member, to save disk.
- `smol resume` reads both `species.tsv` and `species.tsv.gz` (since 25 Oct).

## Experiment scripts

- `experiments/partners/mkconf.py`: configs with blocks of one, two or three sequences in a
  checkerboard. `summ.py RUN [NAME=SEQ...]`: shares, placed/attempted exact copies, share of
  dead molecules that were copied. `batch1.sh`: partnership batch 1. `pairs.txt`: the pairs.
- `experiments/queue2.sh`: establishment replays, published-matched runs, pr0 extensions.
