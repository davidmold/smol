# Next steps (24 Oct)

Nothing is running on aida. Completed and written up: offspring logging, partnership batches 1
and 2, the published comparison (archive analysis, old-rules reruns, baseline extended to 2.5M),
establishment replays, and the big runs (stopped at ~10M).

## Now

1. **What replaces the original pairs?** In the mutation-on partnership runs, trace descendants
   of each partner through `species.tsv` template ancestry, build the placed-copy network among
   abundant species, probe evolved descendants, and remove/rescue descendant lineages. Question:
   does mutual dependence persist, become interchangeable among many species, or collapse into
   autonomous copying with parasites?
2. **The big runs, 1M-10M.** Diversity, length, functional assays and partnership share over the
   full runs (their copy logs are attempts only; they predate OFFSPRINGLOG).
3. **Novelty against drift, rechecked** with the null calibrated from placed births, on runs that
   have OFFSPRINGLOG.

## Later

- Does the partnership accumulate further organisation (third-party dependencies, defence)?
  Longer runs from established replays.
- Replace the top-12 weighted probe with observed-reaction percentages (the paper's measure),
  which needs a per-reaction log.

## Paused

- Side-on worlds (VERTICAL, SHADE, DEPTH, DECAY_CHEM): resumable from checkpoints.
- Template switching, parallel tile updates, plain 3D, TRANSLATE_STOP batch.

## Housekeeping

- The `smol` repo has no remote: push it somewhere for a backup.
- Upstream: open `perf-agentsalign` and `perf-small` after rebasing; follow PR #51.
