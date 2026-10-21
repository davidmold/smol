# Open questions and candidate next steps

Nothing here is agreed work; it's the list of threads left open as of 20 Oct.

## Watching the running experiments

- **Do the partnerships go anywhere?** In `big/strand_400x300`, watch for a third partner, a
  defence against borrowers or parasites, or molecules getting longer. Re-run
  `tools/who_copies_whom.py` every few million steps. If the pairs stay fixed for several
  million steps, that is a stronger negative.
- **Score the big runs again** at later millions with the per-block method
  (`tools/block_novelty.py`), to see whether novelty per area falls over time as it does in small
  worlds. Keeping more checkpoints (spaced out) would let this cover the whole run.
- **Gradient run:** does a lineage specialise to the cold edge, and what happens at boundaries?

## Side-on worlds (towards trees)

The missing ingredient is a link between a molecule and what its copies do for its lineage.
Candidates, simplest first:
1. Copies go only into the copier's own column or the next one, so a lineage owns its column.
2. Growth from within (like a cambium): a new copy is inserted at the copier's position and
   pushes the column up, leaving products behind as the living layer rises.
3. Sediment: a decayed molecule leaves an inert block for a while (niche construction).
4. Test why fragile molecules win (hypothesis: fast turnover frees cells at the lit surface).

## Other ideas raised

- Plain 3D worlds: more neighbours means more mixing; a thin slab (e.g. 400x300x4) is the cheap
  test of whether extra connectivity helps partnerships or invaders.
- Template switching (fusion of two molecules' code during copying).
- Parallel tile updates for much larger worlds (changes the update rule, so a new baseline).
- A pairing/partnership metric on the dashboard.
- `TRANSLATE_STOP` 40-seed batch (deprioritised).

## Roadmap tracks not yet done

- Track 1 replay: resume from checkpoints before takeovers or crashes with 20+ new seeds
  (`smol resume CKPT --reseed S`) to see when outcomes become fixed.
- Origin from a random soup (no seed), or a soup of replicase fragments.
- GPU for analysis (all-pairs alignment and reaction probes).

## Housekeeping

- The `smol` repo has no remote: push it somewhere (e.g. a private GitHub repo) for a backup.
- Upstream: open the `perf-agentsalign` and `perf-small` PRs after rebasing; follow PR #51.
- The roadmap document's early sections (dated 16 Oct) could be refreshed from this folder.
