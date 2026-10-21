# Next steps

## Now (agreed 20 Oct evening)

1. **Partnership isolation** (`experiments/partners/`):
   - Batch 1 (running): four obligate pairs (`experiments/partners/pairs.txt`), each member
     alone, the pair alone, the pair beside the plain copier, the pair with its parasite;
     mutation off and on; 10 seeds; 200k steps on 125x100.
   - Batch 2 (to design from batch 1): spread from rarity into resident copier worlds and into
     evolved communities of different ages; bottlenecks; remove one partner and add it back
     (`tools/ckpt_edit.py`).
   - Key question: does selection preserve a reproducible partnership across generations, or
     does the population merely contain compatible molecules that keep replacing one another?
2. **Measurement**: `OFFSPRINGLOG` and `CKPTKEEPEVERY` are built. Still to do: a functional assay
   panel on the reaction probe (self-copying, partner copying, parasite rejection, reaction
   duration, dependence on other species) and partnership-level complexity.
3. **Published comparison**: matched runs and pr0 extensions queued (`experiments/queue2.sh`);
   get the authors' archive (`biopaperdata.tgz`, CC BY, needs a browser download) for their
   configs and evolved strings, then test their reported behaviours in our probe.
4. **Replays at establishment** (queued): strand0 seeds rerun with a checkpoint every 10k
   steps up to 100k past partnership establishment, then resumed with new seeds to see whether
   establishment is reachable, invasive and persistent.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): runs stopped, resumable from
  checkpoints. Ideas kept: copies into own/adjacent column, growth from within, sediment.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
