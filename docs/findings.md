# What we've learned

In roughly the order it was found. "Survival" means still alive at 1M steps unless stated.
Most batches use the same 40 seeds (RANDSEED 18018, 19019, ... 57057) on a 125x100 grid.

## The chemistry in brief

- 33 symbols: 26 letters (no-ops; runs of them are labels) and 7 opcodes `$ % ^ ? } > =`.
  A letter's complement is 13 places on (A<->N, B<->O, ..., R<->E, F<->S). Opcodes are their
  own complements.
- Two neighbours bind with a probability set by a Smith-Waterman alignment of one's complement
  against the other; in practice everything binds through 2-5-symbol matches, and short exact
  matches bind with probability 1. The molecule whose aligned region starts later becomes the
  **active** and runs its code from the alignment start; the other is the **passive**.
- Pointers I (instruction), R (read), W (write), F (flow), each on either string.
  `^X` toggles a pointer's string (A=I, B=R, C=W, otherwise F); `>X` moves a pointer to F;
  `$label` searches; `?` tests; `=` copies one symbol (with mutation); `%` cuts; `}` ends.
- Grid cells hold one molecule; molecules never move; copies go to a random empty neighbour
  or are lost. Decay 0.0005 per step; a bound pair dies together (old rule, kept).
- Full rules: the specification document (see README).

## 1. The baseline over its first 1M steps (Oct 16-17)

- **One sequence transition, then parasite turnover (to 1M).** In every run the 64-symbol seed replicase is replaced by
  one shorter family of 27-32 symbols (between 29k and 259k steps), which holds to the end.
  After that almost every new family is a parasite (1-16 symbols, copied by others, copying
  nobody). A second replicator family is rare and late (731k, 884k, 971k in 3 of 6 runs,
  never above 11%).
- **Novelty against drift: the original conclusion does not survive the corrected null.**
  Rechecked 24 Oct (`experiments/novelty/`) on the 24 surviving baseline runs, rerun to 1M with
  placement logging (identical runs: same seeds, byte-identical population files). The null
  calibrated from *placed* new species has about a third of the old null's rate of new
  sequences, and then real runs produce **more** new families than drift, not fewer (new
  families per window from 400k-1M: real 4.7-6.7, placed-calibrated null 2.2-3.1, old
  registered-births null 5.1-11.9). Families present: real 12.6-14.7, placed null 13.7-14.3. So
  the drift comparison depends entirely on calibrating the supply of new sequences, and with the
  correct calibration selection *raises* family turnover in these runs (consistent with
  replicator-parasite churn) rather than suppressing it. The old text follows for the record.
  Sequence-family novelty per 100k-step window, against a
  neutral drift null with the same population, turnover and rate of new sequences:
  new families fall 6.4 -> 2.8 in the real runs while the null rises 3.2 -> 11.5; families
  present fall 14.1 -> 10.3 against 11.8 -> 27.7. Families = edit distance <= 0.25.
- **Families are functional.** Reaction-probe profiles overlap 0.59-0.86 within a family,
  0.15-0.38 between families.
- **Survival is close to a coin flip**, 50-60% whatever the rules (baseline `pr0`: 24 of 40).
  Any comparison needs 40 or more seeds.
- **Replicators live on a knife edge over interruption.** A 30-symbol copy takes ~90 steps.
  Independent decay (each partner rolls separately, so reactions break twice as often) killed
  16 of 17 seeds.
- **Binding is not where parasites win.** Replicators copy parasites almost as readily as kin
  (33% vs 38% exact-copy rate). The template has no moves, so no arms race can start.
- **A mechanism can be real and not matter.** The old code's partner-choice bias and turn rule
  caused a dramatic collapse in one seed but made no difference to survival across 114 seeds.
- **Classify by function, not shape**: roles come from the copy log and the reaction probe.

## 2. Strand displacement (Oct 18): harmful as built

A free molecule may break into a reaction (`DISPLACE`). Survival fell as the rate rose:
15-16, 11, 8-10, 0-1, 0 of 40 at rates 0.01, 0.03, 0.1, 0.3, 1.0 (baseline 24). Survivors were
no more novel. `DISPLACE_RESTORE` (undo what the abandoned reaction wrote) separated the
causes: copies were written onto the copier's own tail, so break-ins left it damaged. Restore
lifted survival at 300k from 0 to 8 of 20 at rate 0.3 and 0 to 3 of 20 at 1.0. Damage dominates
at moderate rates, interruption at high rates.

## 3. Binding rule (Oct 18): no clear effect

`BIND_T 3` (P = 1 - exp(-score/3), so longer matches bind more strongly): 17 of 40 survived
against the baseline's 24, within the run-to-run noise at 40 seeds, and replicators still
shrank (`bindT3`, run on chasca). `BINDLEN` lowers binding about equally for every pair.

## 4. Two-sided reactions (Oct 18): deadlock

With `TWOSIDED` both partners run code. Two replicators each copy the other's string onto
their own tail, so each string grows exactly as fast as it is read; reactions run to the
2,000-symbol limit with no product and the grid dies out. Same root cause as displacement.

## 5. Copy strand (Oct 18): robust but more convergent

`STRANDCOPY`: copies are built on a separate nascent strand; the template must be the partner;
`%` releases the strand. New seed `RERERERE^BQQ$=?>%}` (18 symbols).
- `strand0`: 27 of 40 survived. `strandTS` (with TWOSIDED, deadlock gone): 40 of 40.
- Evolution became *more* convergent: last-window novelty 0.6 and 0.2 new families against
  drift nulls of 3.7 and 3.1. Replicators shrank to a median 13 (strand0) and 9 (strandTS).
- `strand0` is the baseline for every copy-strand experiment.

## 6. Evolved seeds and the cross pair (Oct 18)

Seeding with replicators evolved in the baseline: `evo57` (code-first, 26 symbols) 34 of 40,
`evo21` (letter prefix, 28 symbols) 24 of 40, against the baseline's 24. A designed
**cross-replicating pair** (A copies only B, B copies only A; `designs/cross.txt`) survived
**40 of 40** and kept its length (median 38 vs 39 at the start), with both types in the top
three at 1M in 32 of 40. Its novelty was lower than the baseline's (2.4 vs 4.7 new families).
The only design so far that held its length.

## 7. Spontaneous partnerships (Oct 19-20)

Obligate mutualisms evolve on their own in the copy-strand chemistry (cf. Barricelli's
symbiogenesis).
- In small `strand0` runs: **borrowers** with no `=` of their own, which use `^A` to run the
  partner's code. Present in 27 of 27 survivors, established at a median of ~160k steps,
  arising from a truncated fragment plus a `^B`->`^A` mutation. Tool: `tools/partners.py`.
- In `big/strand_400x300`: **complementary-prefix twins**, e.g. `RE^BQQ$=?>G}` <-> `ER^BQQ$=G?>G}`:
  zero self-copies, ~80-90k copies each way per 20k steps. ER-ER can't bind; RE-RE gets stuck.
- By 0.9M steps in the same run: 65% of exact copies happen inside two-way pairs, only 3% are
  self-copies. Top pair `S^BO$=J?>G}` <-> `F^O>` (~2.1M copies each way, no self-copies):
  `^O` moves the flow pointer to the partner and `>` jumps into its code. Nearly half the
  copy-strand population lacks `=`, and most of those are partners, not passengers.
- The pair together is ~15 symbols, no longer than the seed: a cheaper division of labour,
  not added complexity.

## 8. Gene and machine: a von Neumann constructor (Oct 19)

`TRANSLATE`: with `^C` toggling the unused W pointer, `=` reads letter pairs through a balanced
676->33 codon table and writes the translated symbol (opcodes pass through). Seed: gene G
(100 letters, does nothing as code) + machine M (42 symbols) that copies G blind and
translates G into a new M. Tool: `smol-codon` (encode/decode/table). Design files
`designs/gm.py`, `gm.txt`, `gm.conf`.
- `gm` (G+M alone): survived 40 of 40, G length stable (94-102 letters).
- `gmc` (G+M with 25 plain copiers 36 cells away): G+M survived in 0 of 40. Copiers win on
  speed and space.
- Seed 28028 crash: a `^S$S` parasite exploited mid-code binding plus opcode passthrough in
  translation; a `^F`->`^L` resistance mutation then evolved. `TRANSLATE_STOP` (a non-pair is a
  stop codon) closes the passthrough; a 40-seed run with it has not been done.

## 9. Damage and "cells" (Oct 19)

`DAMAGE`: each symbol of every molecule is replaced at random with a small chance per step
(an environmental hazard). The idea: in biology the machine's robustness, not its speed, is
its advantage. Then a "cell": a cell-cycle machine M (`designs/cell.py`, 44 symbols) that holds
its gene G while copying and translating it, with `DAMAGE_SHIELD_BOUND` protecting molecules
inside a reaction.
- Plain copiers tolerate damage up to ~3e-4; shielded cells die at 1e-4.
- The cell lost to copiers at every damage rate tried (`dmg/`, `cell/`, `cellhi/`).
- David's steer at this point: having copiers is fine; the goal is open-endedness, not making a
  chosen design win. Fewer, bigger, longer runs with a live viewer.

## 10. Big worlds (Oct 20)

One run each on 400x300 (120,000 cells, 9.6x the area), still running.
- **Sturdier, not more inventive.** Classic at 1M: 8,992 species alive, 365 effective species
  (inverse Simpson), top species 1.5%, against 200-900, 10-40 and 10-30% in small runs. No crash
  so far, against ~40% extinction in small worlds. But the families are the ones small runs
  produced.
- **Novelty per area is unchanged.** Scored per 100x100 block (about one small world),
  counting a family as new only if unlike anything abundant anywhere in the grid 100k steps
  earlier: 1.8-4.0 new families per block per 100k, against 2.1-3.1 per small run scored the
  same way. Whole grid ~30 per 100k, about 12 small worlds' worth. Only 800k-1.2M could be
  scored per block (positions come from checkpoints, of which 5 are kept).
- **The novelty metric is not scale-free.** At 60k-100k molecules no species in the drift null
  reaches the 1% bar, so the whole-grid null is empty after 300k. Big worlds need per-block
  scoring.
- **Chemistries diverge.** Classic stays with self-copiers (53% of copies; 3% go to code-less
  molecules; mean length ~28). Copy strand goes to partnerships (65% of copies in two-way pairs;
  mean length ~10). The gradient run has fewer pairs (41%), more self-copying (18%) and about a
  third of the uniform run's effective diversity; copy code is least common in the warm middle.
- The classic chemistry needs the original 65-symbol replicase; the short copy-strand seed
  dies within 15k steps without `STRANDCOPY`.
- Taken together: the "bigger worlds" leg of the clean negative, supporting the
  expensive-ratchet view (one run per condition).

### 10b. The big runs to ~10M steps (analysed 24 Oct; `experiments/big/`)

Copy logs in these runs are attempts (they predate `OFFSPRINGLOG`). Measures: `copy_mix.py`
(exact copies by kind, 20k-step window before each 0.5M), `length_traj.py`, `assay.py` (top 12
species reacted pairwise every 1M).

- **Copy strand, uniform:** the partnership ecology persists for all 10M steps. Mutual copying is
  52-71% of exact copies throughout; one-way copying (templates copied but not copying back) falls
  from ~30% to ~15%; templates without copy code fall from ~50% to ~25-30%; self-copying rises
  from ~3% to 20-29% after 2M. Among the top species, parasitic cross-pairs fall from 0.51 at 1M
  to 0-0.12 from 6M. Lengths stay ~8-10 (median).
- **Copy strand, gradient:** similar but more parasitic: mutual 55-68%, one-way ~30%, self 5-15%;
  median length drifts down from 11 to 5-6.
- **Classic: a late transition to longer, processive replicators.** Until ~3M the classic world is
  as before (self-copying 55-62%, one-way ~31-35%, median length ~30). Between 3M and 4.5M copying
  collapses and recovers, and from 4.5M on: self-copying 94-96%, one-way 4-6%, mutual ~0; median
  length rises from ~30 to 60-75 (90th percentile 100-260); the population is extremely diverse
  (top 12 species are 3-4% of molecules). The new replicators (e.g.
  `EKM>BHCOM^$EOC>C%M$==C=?H>PFYPBXCOPB$E>`, ~40 symbols) **do not finish their reactions**: in the
  probe one binding makes copy after copy of the partner (46 exact copies in 5,000 steps, about
  one per 110 steps) and is still running at the step limit. This processive (amplifying) mode of
  reproduction did not appear in any small run or in the published runs, and it is the first
  sustained increase in molecule length we have seen. One run.
- **How it arose** (`mix_classic_100k.txt`, `top_classic.txt`: top species every 100-500k,
  self-reacted in the probe):
  1. To ~3.0M the world is ~26-symbol self-replicators (`??>C$BU^BQ>C$=?W>F$B%R}XHO`): one copy
     per ~100-step reaction.
  2. 3.1-4.0M, a crisis: a family ending `...^B>C%T$==?>` (doubled `=`, often with tandem repeats,
     31-63 symbols) takes over, together with the fragment `%T$==?>`. None of them makes an exact
     copy of itself in the probe. Exact copying falls ~20-fold (2.1M to ~0.1M per 20k steps) and
     98-99% of what remains is one-way: templates copied by species they don't copy back.
  3. From ~4.15M, molecules from that family with a tripled copy (`...?Q>C%FTB$===?>^B...`,
     68-74 symbols) appear that copy processively (26-30 exact copies per 5,000-step reaction). By
     4.5-5.0M they and relatives (42-74 symbols) are among the top species, self-copying recovers
     to 0.8-0.95 of exact copies, and one-way copying drops to 2-6%.
  So the processive replicators came out of the parasitic family that had displaced the old
  replicators, not from the old replicators directly. This is the first observed case of a
  transition after the initial one, but it is a single run.

## 11. Side-on worlds: towards trees (Oct 20, exploratory)

`VERTICAL`: x is horizontal, y is height with the ground at the bottom; a copy can only go on
the ground or on top of an occupied cell; when a molecule decays, the column above falls by
one. Light makes reactions faster higher up (`RATE_MIN`/`RATE_MAX` by height) or comes from
exposure (`SHADE`: rate x0.5 per molecule stacked above). `DECAY_CHEM 10`: each symbol has a
random fixed stability weight (0.07-4.5 with seed 1; tough: J G ? C > =, fragile: F U } ^) and a
molecule decays at DECAY times the mean weight of its symbols. `DEPTH n`: n slices side by side.
All copy-strand chemistry.
- Without falling, column bases decayed and left a floating canopy. With falling, a 120x60
  world stayed a turf ~4 high, never used the light, and died at 547k.
- 400x12 (light spans reachable height): fills up, survives; the winner `SFSF^B$=?>}` is
  *fragile* (decay weight 1.78).
- With shading: in one seed of three the surface was stable junk (`RER`, weight ~0.45, capping
  155 of 400 columns) shading the copiers underneath, which keep copying it on top of
  themselves: a smothering lid, not a tree. Other seeds uniform. 400x12 shading run died at 379k.
- 5-slice slab with shading: spreads through all slices; one of three test seeds died; fragile
  molecules (weight 1.5-1.7) dominate again.
- Diagnosis: nothing yet rewards a molecule for what its copies do. A copier gains nothing from
  where its copies land and will copy whatever is next to it, including the lid.

## 12. The published runs, analysed with our tools (Oct 20)

The authors' archive has full snapshots every 20k steps of their 8 surviving 2M-step runs, and
their configs. Runs 1-4 use exactly our baseline config (125x100, 121 seeds in a block, DECAY
0.0005, MUTATE 0.0002). Details and numbers: `experiments/published/`.
- **Our reaction probe reproduces their reported reaction times** on their own evolved strings
  (run 2, self-reaction of the top species): 184, 263 and 225 steps at 600k, 680k and 1.5M
  against their 183, 262 and 224. Their parasitic share reaching ~0 at 1.5M in run 2 also shows.
  So the chemistry matches at the reaction level.
- **Over the first 1M, their runs and ours look the same** (medians across runs): 5-8k
  molecules, median length falling from 64 to ~25-30, top species ~26-35 symbols, 80-90% of
  same-species pairs self-replicating, 60-75% of productive cross-species pairs parasitic,
  reaction times 130-170 steps.
- **The change is in their second million, and it is in interactions, not length.** Parasitic
  share falls from ~0.7 to 0.31-0.55 (0.31 at 2M), mutual copying rises from ~0.3 to 0.45-0.69,
  reaction time rises from ~140 to 180-205 steps, and population from ~7k to ~10k. Median and top
  lengths do not rise (22-30).
- **Our runs past 1M (21 Oct).** The 24 baseline survivors continued to 2.5M (`runs/pr0x`):
  18 reached it, 6 died after 1M. The published configs rerun with the old rules (`runs/pub`,
  5 seeds each): 5 of 15 reached 2M (box 1/5, random pairs 3/5, slot 1/5; they had 8 of 20).
  Both show a **weaker version** of the published shift (medians across survivors):

  | | parasitic share | mutual share | reaction time | N |
  |---|---|---|---|---|
  | published, 1.0-1.1M -> 1.9-2.0M | 0.71-0.79 -> 0.31-0.56 | 0.21-0.29 -> 0.43-0.69 | 167-188 -> 162-184 (peak 205) | 8k -> 10k |
  | our baseline rules, 1.1M -> 2.4-2.5M | 0.70-0.78 -> 0.60-0.66 | 0.22-0.30 -> 0.34-0.40 | 136-184, no trend | 7.5k -> 10.7k |
  | published configs, old rules, 1.0-1.2M -> 1.9-2.0M | 0.66-0.78 -> 0.61-0.62 (0.44-0.47 at 1.5M, 1.8M) | 0.23-0.34 -> 0.38-0.39 | ~150 -> 178-206 | 9-10k -> 10.3k |

  Mutual copying rises and parasitism falls in all three, but less in ours, and reaction time
  rises only with the old rules. Length does not rise in any of them. Small numbers of
  surviving runs (5-18), so this is a direction, not a measured difference.
- Caveat: "parasitic" and "mutual" here are weighted by abundance and binding probability over
  the top 12 species (`tools/assay.py`); the paper's percentage is over observed reactions.

## 13. Partnership isolation (Oct 21-22)

Four obligate pairs from the big runs (`experiments/partners/pairs.txt`; in the probe neither
member binds its own kind and each copies the other), copy-strand chemistry, 125x100.
Batch 1 (`agg1.txt`, 256 runs, 200k steps) and batch 2 (`agg2.txt`, 464 runs).

- **Members alone die out** (0 of 16 runs); the pair alone fills the grid (40/40).
- **Mutation off, the pair is a stable ecology.** Beside an equal block of the plain copier, three
  pairs exclude it completely and g1 coexists (copier 15%). With a parasite the pair holds ~30%
  of the grid against ~70% parasite in all 40 runs.
- **Mutation on (200k), the original pairs become a minority** (4-27% of molecules), replaced by
  short fragments (s1-s3) or newly evolved ~15-symbol self-copiers (g1, and where the copier was
  present). 117 of 120 worlds survive. Whether descendants keep the partnership is not yet
  measured.
- **Invasion of a resident copier world from rarity is chemistry-specific.** s1 and g1 never
  establish (0/40): the copier can't bind s1, and copier-g1 complexes stall. s2 and s3 take
  over from 20 molecules (10/10, 8/10), the copier going extinct, but the invasion often crashes
  the whole world first (from 2 molecules: 5/10 and 8/10 worlds died). In the probe s2A is copied
  by the copier without copying back, and copier-s3A reactions never end, so copiers are
  parasitised or locked up. The copier never invades a pair world (0/40).
- **Bottlenecks:** cut to 1%, pair worlds survive 35/40 against the copier's 10/10; cut to 0.1%,
  4/40 against 1/10. Obligate pairs are only slightly more fragile.
- **Rescue:** with one partner removed, adding 20 back after 0, 2k or 5k steps restores the pair
  in 60/60 runs; after 10k in 5/20; after 20k in 0/20 (the remaining partner has died out).
- **Accessibility vs introduction** in evolving communities (strand0 replays, MUTATE 0.0002,
  150k steps after the checkpoint). A partnership counts as established if, in some 10k window,
  at least 30% of placed exact copies happen in mutual cross-species pairs (counts at a 10%
  threshold are in agg2.txt and are higher):

  | seed (orig. establishment) | checkpoint | replays establishing | own pair introduced | foreign s1 introduced |
  |---|---|---|---|---|
  | 18018 (110k) | 30k / 60k / 90k | 5/10, 6/10, 10/10 | 2/3, 2/3, 3/3 | 3/3, 3/3, 3/3 |
  | 20020 (70k) | 20k / 40k / 60k | 1/10, 4/10, 10/10 | 3/3, 3/3, 3/3 | 2/3, 2/3, 3/3 |
  | 23023 (100k) | 30k / 60k / 90k | 3/10, 2/10, 7/10 | 3/3, 3/3, 3/3 | 3/3, 3/3, 3/3 |

  Provisional reading: supplying pairs makes substantial mutual copying easier to obtain in the
  three tested communities, and spontaneous establishment becomes likelier as the community ages.
  That supports an accessibility limitation; it does not exclude establishment or persistence
  barriers, because the criterion counts any reciprocal pair reaching the threshold in one 10k
  window, not the introduced partnership persisting (the introduced sequences were still present
  at the end in 42 of 54 intro runs, but presence is not continued mutual dependence). 3 seeds.
- **What replaces the original pairs (24 Oct; `lineage_agg.txt`, `dep_agg.txt`).** Each species
  is assigned to the founder it descends from through template ancestry (`lineage.py`).
  - In the mutation-on runs at 200k, **the partnership persists through descendants**: in 90 of
    117 surviving runs at least half of all placed exact copies are made between the A and B
    lineages in both directions; same-species self-copies are ~0 for s1-s3. The molecules I had
    called "other" are mostly the borrower's own descendants: mutants that lost the jump (e.g.
    `F^OQ`, `F^OR`, `F^O` from `F^O>`), copied by the A lineage and making nothing (44-86% of the
    B lineage). The working core is carried by 1-5 effective species pairs. Exceptions: g1 is
    mixed (some autonomous copying), s2 with its parasite is dominated by the parasite lineage,
    and where the plain copier was present it persists alongside (15-36%).
  - **Removal tests on the evolved communities** (from 200k, mutation on, 50k more steps):

    | pair | unchanged | B lineage removed | A lineage removed | B removed, 20 of an evolved B added |
    |---|---|---|---|---|
    | s1 | 10/10 | 0/10 | 0/10 | 10/10 |
    | s2 | 10/10 | 7/10 | 0/10 | 9/10 |
    | s3 | 10/10 | 2/10 | 0/10 | 10/10 |
    | g1 | 6/8 | 4/8 | 0/8 | 5/5 |

    The B lineage is entirely dependent on A. The A lineage still depends on B in s1 and mostly
    in s3; in s2 and g1 it often survives B's removal, either by self-copying variants
    (`S^OF^B$=?L>}` in g1, `>$$^BQQ$=?>H}` in s2) or by a **new partnership inside the A
    lineage** (s3: `SE^BQ$=I?>G}` <-> `FS^BQ$=I?>G}`, complementary prefixes, zero self-copies).
    An evolved B descendant restores the partnership in 34 of 35 runs.
  - So under mutation the organisation is maintained by descendants rather than by the
    original sequences, with a fringe of passengers shed by the borrower lineage; the copier
    lineage sometimes carries standing variation that can replace its partner, either by
    becoming autonomous or by re-forming a partnership with itself.
- **Against the three possibilities** (hard to reach, can't invade, can't persist): for
  the tested pairs form self-maintaining, mutually dependent ecologies; mutual copying is
  provisionally easier to obtain when supplied (accessibility); against a resident autonomous
  copier, invasion depends on the chemistry of the encounter; under mutation the *original
  sequences* become a minority, but whether the *partnership* persists through their descendants
  is not yet established (next: ancestry plus removal tests on descendants). Whether
  partnerships accumulate further organisation is the open question.

## 14. Processive replicators: detection, screen, transplants, long small runs (Oct 24-25)

**Definition** (`tools/transition_scan.py`). A species is *processive* if, in a probe
self-reaction (mutation off, 5,000-step limit, `smol-react --self`), it binds itself and makes at
least 3 exact copies of itself in that one reaction; reaching the step limit alone does not
count. *Detection*: the first sampled step with a processive species of >= 20 molecules.
*Establishment*: processive species hold >= 10% of molecules for >= 3 consecutive samples. Every
species with >= 5 molecules is probed (up to 4,000 per sample), not just the most abundant; the
scan also reports cumulative placed (or attempted) copies, so emergence can be set against
reproductive activity as well as elapsed steps.

- **Validation on the big classic run** (every 500k, species >= 20 molecules): detection at
  3.0M (2 processive species, 0.1% of molecules: a `...%T$==?>` family member, i.e. the precursor
  was already present, rare, before the crisis), 0.7% at 4.0M, establishment at 4.5M (21%), and
  59-64% of all molecules from 6M on.
- **Screen of existing histories** (every 100k, species >= 5 molecules): no processive species
  in any of 87 runs: our 40 small baseline runs (0-1M), the 24 survivors continued to 2.5M, the
  15 old-rules reruns of the published configs (to 2M) and the 8 published runs (to 2M). All of
  these stop before the 3M at which the big world first showed them, so this is not evidence of
  a size requirement (time confound).
- **Transplants** (`experiments/processive/`, 123 runs of 300k steps, mutation on, 3 replicates;
  20 molecules in a cleared 5x5 patch; `tp_agg.txt`). Sequences: early processive E1, E2 (from
  5.0M; 42 and 50 symbols; 45 and 34 copies per reaction), mature M1, M2 (9.8M; 40 and 61
  symbols; 46 and 34), and an ordinary pre-transition replicator ORD (2.8M; 26 symbols; 1 copy).

  | recipient | E1 | E2 | M1 | M2 | ORD | resident only |
  |---|---|---|---|---|---|---|
  | empty 125x100 world | 3/3 (0.84) | 3/3 (0.80) | 3/3 (0.75) | 3/3 (0.87) | 0/3 | - |
  | 5 evolved communities (seeds 29029, 36036, 40040, 47047, 56056) | 15/15 | 14/15 | 15/15 | 15/15 | 0/15 | 0/15 |
  | evolved community seed21021 | 1/3 | 0/3 | 3/3 (0.42) | 1/3 | 0/3 | 0/3 |

  Cells: runs with processive species >= 10% of molecules at 300k (median final share). The
  communities are small classic worlds evolved for 2.5M steps (from `pr0x`). Processive
  replicators established in 64 of 72 community introductions, usually within 50-150k steps,
  reaching 42-88% of molecules; no control ever became processive. One community (seed21021)
  largely resisted (5/12). Descendants were tracked as processive species by the scan, not as
  the introduced sequence.
- **Reading:** once they exist, processive replicators need neither a big world nor an empty
  one to spread: they invade most evolved small communities from 20 molecules. So the open
  question is whether they can *arise* in small worlds and how often in big ones (below). A concrete hypothesis from the original
  history: the parasitic crisis through which they arose may usually destroy small populations.
- **Bug found and fixed:** runs resumed from an edited checkpoint with `GZIPLOGS` had two species
  files (the copied `species.tsv` and `species.tsv.gz` for new species), and the scan read only
  the first. The first transplant summary (no establishment in communities, "boom then bust")
  was an artefact of this and was withdrawn the same evening. `ckpt_edit.py` now writes the gz
  file when `GZIPLOGS` is on, and the analysis tools read both files.

- **Small worlds to 10M (25-26 Oct; `experiments/processive/scan_long/`).** Scored every 100k
  (species >= 5 molecules), from 2.5M for the survivors. All 38 runs are finished.

  | set | reached 10M | extinct (step) |
  |---|---|---|
  | `pr0xx`: 18 survivors continued 2.5M -> 10M | 13 | 5 (5.19, 5.53, 5.72, 6.08, 6.55M) |
  | `long`: 20 fresh runs 0 -> 10M | 8 | 12 (6 before 1.5M; 2.1, 3.1, 3.9, 6.0, 7.0M) |

  **No processive species reached detection (>= 20 molecules) in any run.** Processive species
  were seen below that, briefly, in three runs: long seed66066 at 3.3M (10 species, 0.7% of
  molecules), pr0xx seed44044 at 3.1M (1 species) and pr0xx seed21021 at 10M (3-7 species,
  0.3%, largest 16 molecules). The first two were gone at the next sample (100k later); the
  seed21021 sighting is the run's last observation, so its fate is unknown. That run is being
  continued to 13M with checkpoints every 50k (26 Oct). A finer scan of 9.0-10M (1,000-molecule
  samples every 20k, `scan_21021/`) shows the lineage appearing between 9.96M and 9.98M, with a
  33-molecule species at 9.98M: by the definition above that is a **detection** (>= 20
  molecules) in a small world, missed by the 100k-step scan; it is at 0.6-0.7%, far from
  establishment.
- **The small-world variants** (`seqs_sizes.json`). seed66066 at 3.3M: 38 symbols, 46 copies
  per self-reaction, and the same layout as the big-world processive replicators (a `>C%`
  copy block); 17 + 14 + 2 molecules of three variants. seed21021 at 10M: 95 symbols, 24 copies
  (16 and 8 molecules of two variants). seed44044 at 3.1M: one species of 462 symbols making
  5 copies (5 molecules).
- **A blind spot in the scan, and the check.** The default scan only probes species with >= 5
  molecules. In 7 of the 21 runs that reached 10M (long 64064, 66066, 68068; pr0xx 21021, 22022,
  29029, 44044) these cover only 0-20% of molecules for long stretches: almost every molecule is
  a unique variant. The big classic world went the same way during its crisis (cover 0.15-0.3 at
  3-4.5M), and the big replicates are mostly in that state too (cover 0.08-0.4). So a new option,
  `transition_scan.py --sample K`, probes the species of K molecules drawn at random (weighted by
  abundance), which estimates the processive share without regard to rarity
  (`experiments/processive/scan_sample/`). On the original big run it gives 0 at 3.0M (below
  the sample's resolution), 1.1% at 4.0M, 30% at 4.5M and 80% at 6.0M. That matches the
  abundance scan on the timing of the transition but not the size: the abundance scan gave
  ~21% at 4.5M and 59-64% from 6M, consistent with it missing rare processive variants. So the
  two agree on when, not on how much. With K = 300 every 100k over 2.5-10M in the 7 low-cover small runs it finds
  processive molecules only at the three sightings above (3, 4 and 1 sampled molecules: 1.0%
  at 3.3M in 66066, 1.3% at 10M in 21021, 0.3% at 3.1M in 44044, so these shares are rough) and
  none at the other 525 of 528 samples; a scan probing up to 3,000 species of any size every 500k
  (`scan_min1/`) agrees. What this rules out: a missed *established* population. With 300 draws
  and no hits, the 95% upper bound on the share at that sample is ~1% (~0.3% for the 1,000-draw
  big-world samples), assuming the probe classifies correctly; brief appearances between samples
  are not ruled out. The probe now stops on a failed `smol-react` call instead of recording a
  negative; none of the ~190k cached classifications came from a failed probe.
- **Reading, provisionally.** Processive variants do arise in some small worlds and (in two of
  three cases) vanish. One hypothesis is that the barrier is getting a foothold rather than
  persisting: the transplants show that *selected big-world genotypes*, supplied as 20 adjacent
  molecules in a cleared 5x5 patch, spread in most evolved small communities. But that protocol
  supplies partners and space, and the naturally arising small-world variants may be poorer
  genotypes or a poor fit to their community, so the transplants don't yet say why they fail.
  Next: introduction size (1, 2 adjacent, 5, 20, with matched clearing) and the small-world
  variants themselves alongside the big-world ones. The two small cohorts (20 fresh runs; 18
  survivors continued from 2.5M) are differently selected and only 21 runs reached 10M, so they
  are kept separate in any frequency estimate.
- **Activity, not just steps.** The original big run logged only attempted copies (products,
  placed or not): 6.3 x 10^8 by its first detection at 3.0M and 1.2 x 10^9 by establishment.
  The fresh small runs that reached 10M made 1.8-3.9 x 10^8 attempts (the continued survivors
  1.2-3.9 x 10^8 over 2.5-10M alone), about a third to two thirds of the big world's activity
  before detection (`scan_long/attempts.txt`). So small worlds have been watched for three times
  the steps but somewhat less copying, and one big world is a sample of one: "not in small
  worlds" is suggestive, not yet a rate comparison. (Placed products are 14-29% of attempts in
  these runs; one run, pr0xx seed36036, places 94%.)
  Extinction is the commoner outcome (17 of 38 so far, about half of them late, after 2M).
  One late extinction looked at in detail (pr0xx seed47047, 5.19M) was not preceded by a
  parasitic crisis: exact copying held at ~1 x 10^5 per 250k, and one-way copying *fell*
  (0.36 -> 0.12 of copies) while self-copying rose (0.55 -> 0.71) just before the end.
- **Big replicates so far (26 Oct; `scan_big/`, `scan_sample/`).** None of the three 400x300
  classic replicates has processive species yet: none in the abundance scan (>= 20 molecules) or
  in random samples of 1,000 molecules every 500k, to 4.0M (seed19019), 4.5M (seed20020) and
  3.5M (seed21021) at the last scan (26 Oct; the runs themselves are further on); the original
  run had 1.1% at 4.0M and 30% at 4.5M. Three replicates without a transition would mean the
  original result was not reproduced in three runs, not that it is rare: even at a 50% chance
  per run, three misses happen 12.5% of the time. seed20020's history
  differs: 6.4 x 10^6 species ever by 2.9M against 7.9 x 10^7 for seed19019 by 2.1M, at
  comparable copying (4.8 x 10^8 attempts by 2.5M), and its cover fell from 0.77 to 0.15 over
  3-4.5M, like the original crisis. Running on to 10M.
- **Runs lost to memory, and the fix (25 Oct).** The simulator kept every sequence ever made
  (including never-placed products) as a key in memory, so a small run grew to ~15 GB by 6-9M
  steps and a big one to ~13-19 GB by 2M. Five runs were killed by the out-of-memory killer
  (pr0xx 21021, 24024, 29029, 44044 and big seed21021); one more (long seed66066) was stopped
  at a checkpoint to protect the rest. All were resumed and the small ones finished. The species table is now keyed by a 128-bit hash of the
  sequence (outputs byte-identical to before on a test run), and the new switch
  `SPECIESFORGET` lets a run forget sequences absent from the world at a checkpoint (an
  extinct sequence that reappears gets a new id; the dynamics are unchanged, verified). The
  big replicates continue with it on (seed20020 from 4.9M, 19019 from 2.2M, 21021 from 1.7M);
  the small runs without it. All stopped runs resume from
  their last checkpoint after trimming their logs back to it (`tools/trim_logs.py`; a trimmed
  resume reproduces the uninterrupted logs exactly, apart from `lives.tsv`, whose per-molecule
  bookkeeping is not checkpointed).
- **A second resume bug found and fixed:** `smol resume` read only a plain `species.tsv`, so a
  GZIPLOGS run resumed with its species only in `species.tsv.gz` started with an empty table
  and gave reappearing sequences new ids. It now reads both. No reported run was affected
  (every earlier resume had a plain copy of its species file).

## Caveats on the measurements (review, 20 Oct)

- **Copy attempts are not offspring.** Until `OFFSPRINGLOG`, `copies.tsv` counted a product
  (and `species.tsv` registered its species) before checking for room to place it. In a full
  copy-strand world only ~4% of attempts are placed; in a filling classic world ~39%. All copy
  counts and partnership shares in sections 7 and 10 are attempts. The drift null's rate of new
  species is calibrated from `species.tsv`, which includes never-placed products, so it is too
  high.
- **The novelty metric is narrow.** It compares each window with the previous one only, so a
  family that returns counts as new; a small sequence change with a large functional effect stays
  inside a family; the null fixes population size, has no space and uses its own mutation
  process, so real-vs-null differences are not purely selection.
- **The big runs are one run per condition.** Blocks are spatial sub-samples, not replicates.
- **"Expensive ratchet" is not yet decomposed.** Useful innovations might be hard to reach,
  reachable but unable to invade, or able to invade but unable to persist. These call for
  different experiments (replay, introduction at low abundance, persistence after invasion).
- **Run length.** Hickinbotham, Stepney & Hogeweg (2021) report replicator length and reaction
  time rising, and parasitic reactions falling to near zero by 1.5M, in the second half of 2M-step
  runs (12 of 20 runs went extinct). Most of our runs stop at 1M.

## The pattern across everything

Short and fast wins on length whenever nothing stops it: changes to binding, interruption,
reaction symmetry, copying mechanics, translation, damage and space all left that pressure in
place, and only the designed cross pair held its length. But length is not the whole story:
interactions keep changing after lengths settle (parasitism down, mutual copying up), and the
one new organisation that evolved by itself, obligate partnerships, forms self-maintaining
ecologies. Whether it builds on itself is the current question.
