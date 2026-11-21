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
  surviving runs (5-18), so this is a direction, not a measured difference. Rechecked with
  COMPAT_DELETION in section 29: still no significant difference, still a weaker shift.
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
  establishment. **It then vanished:** 0.2% at 10.00M and none in any of the 102 samples from
  10.02M to 12.04M (`scan_21021/post10M.txt`). So the small-world lineages so far arise, reach a
  few dozen molecules at most, and are gone within ~50-100k steps.
- **Introduction size (26 Oct; all 270 runs; `experiments/processive/sizes.py`,
  `sizes_agg.txt`).** Clear the same 5x5 patch in an
  evolved community (pr0x at 2.5M), add N molecules of one sequence, run 200k steps with
  mutation; established = processive share >= 10% at the end. Cells: established / runs.

  | recipient | sequence | N = 1 | N = 2 (adjacent) | N = 5 | N = 20 |
  |---|---|---|---|---|---|
  | seed36036 (permissive) | E1 (big world, 42 symbols) | 2/12 | 10/12 | 6/6 | 3/3 |
  | | S66 (small world, 38) | 0/12 | 6/12 | 6/6 | 1/3 (3/3 persist) |
  | | S21 (small world, 95) | 0/12 | 5/12 | 4/6 | 3/3 |
  | | S44 (small world, 462) | 0/12 | 0/12 | 0/6 | 0/3 |
  | seed21021 (resistant) | E1 | 0/12 | 1/12 | 0/6 | 1/3 |
  | | S66 | 0/12 | 0/12 | 0/6 | 0/3 |
  | | S21 | 0/12 | 1/12 | 0/6 | 1/3 |
  | | S44 | 0/12 | 0/12 | 0/6 | 0/3 |

  "Established" here means a processive share >= 10% at the 200k endpoint, a looser test than
  the scan's (>= 10% for three consecutive samples); neither means takeover. S66 at N = 20
  persisted in all three permissive runs but crossed 10% in only one. Clearing alone (NONE)
  never produced processive replicators (0/6). Outcomes were mostly all or nothing: a single E1
  either took over (~87% of molecules within 40k steps; 2 runs) or was gone within 20k (10 runs).

  **Design caveat.** Every treatment clears the same 5x5 patch, so a single molecule starts
  with all eight neighbouring cells empty and must wait for residents to grow back into
  contact, while an adjacent pair has a reaction partner at once. This is a fair comparison of
  colonisation from different founding groups, but not of how a natural mutant fares: a mutant
  usually appears next to the molecules that made it. And every introduction uses the same
  location, so differences between recipients are partly differences between two local
  neighbourhoods.

  What the table supports, under this protocol:
  - **Founding group:** an adjacent pair clearly helps (E1 2/12 -> 10/12; S66 0/12 -> 6/12; S21
    0/12 -> 5/12), consistent with an encounter limit. For the small-world variants a pair
    succeeds about half the time, not usually.
  - **Genotype:** S66 and S21 can establish in a small world; S44 (462 symbols) never did. This
    does not show that length is the reason, or that S66 and S21 are as competitive as E1.
  - **Recipient:** a large difference between the two snapshots tested. seed21021 admitted 4 of
    132 introductions (4 sequences x 33; the 3 clearing-only controls aside), all E1 or S21 at
    N = 2 or 20. seed21021 at 2.5M is an ancestor of the 9.98M community in which S21's lineage
    arose and was lost, not the same ecological state.
  So these results narrow the possibilities but do not yet say why the natural seed21021 lineage
  failed in its actual surroundings. Next: introductions without clearing (replacing a resident),
  at several locations, into checkpoints near S21's emergence as well as at 2.5M, with the first
  encounters of each newcomer logged.
- **Introductions without clearing (27 Oct; `experiments/processive/natural.py`,
  `natural_agg.txt`, `natural_families.txt`).** One processive molecule replaces the resident of
  a randomly chosen occupied cell, or two replace that cell and its right-hand neighbour; nothing
  else is touched (a replaced resident's partner is freed). 16 locations per recipient, shared by
  both sequences and both N, so treatments are paired. The newcomers and their descendants are
  traced (`TRACESEQ`/`TRACEAT`). 200k steps, mutation on. Established = processive share >= 10%
  at 200k *and* the top processive species belongs to the introduced family (by its copy block:
  E1's `>C%TB$=E===?>`, S21's `>C$====?`). Cells: N = 1 / N = 2, out of 16 each.

  | recipient | E1 | S21 |
  |---|---|---|
  | seed36036 at 2.5M | 1 / 5 | 0 / 9 |
  | seed21021 at 2.5M | 3 / 4 | 2 / 3 |
  | seed21021 at 9.9M | 1 / 2 | 0 / 0 |
  | seed21021 at 10.0M | 0 / 3 | 0 / 0 |

  - **seed21021 at 2.5M is not a resistant community.** Over 16 random locations without
    clearing it admitted 12 of 64 introductions, against 15 of 64 for seed36036. The cleared
    single-location result above (4 of 132) was largely a property of that protocol and place.
  - **A second molecule helps in seed36036** (E1 1 -> 5, S21 0 -> 9) but hardly in seed21021
    at 2.5M (3 -> 4, 2 -> 3). Without clearing, E1 singletons did no better than from a cleared
    patch (1/16 here, 2/12 there), and E1 pairs worse (5/16 against 10/12).
  - **S21 never established in its own late community:** 0 of 64 at 9.9M and 10.0M, against
    6 of 64 for E1 there. Introduced singly, 13 and 12 of 16 S21 founders (9.9M, 10.0M) were
    first bound as the *passive* partner by a resident (the resident runs its code on the
    newcomer), and 8-9 of 16 had their own sequence altered by a reaction before dying. In
    seed21021 at 2.5M most founders' first reaction was instead as the *active* partner (10-11
    of 16), and few were altered (2 of 16).
  - **Native establishment in two replays.** In 2 of the 32 E1 runs from 9.9M the introduced
    lineage died out but a native S21-family processive population reached 16% and 31% by 10.1M
    (e.g. `TDNP>BO^$O%HK?>C$====?>MZ$I>`, 28 symbols, 80 self-copies per reaction): the family
    that produced S21 arose again from the same state and, this time, took hold. None arose from
    10.0M (0 of 64).
  - **Clean replays** (`replays.py`, `replays_summary.txt`; 24 reseeded 200k-step replays from
    each checkpoint, nothing introduced). From 9.9M, 1 of 24 established a native processive
    population: `Q>BO^$O%H>C$===?O>A$$T$TQ>` (26 symbols), detected at 9.94M with 42 molecules
    (0.5%), 10% by 10.04M, 36% by 10.1M. Two more replays had 1-3 processive species briefly;
    the other 21 none. From 10.0M, where the original lineage was present (16 + 8 molecules),
    0 of 24 established; it was detected again in one replay (10.04M) and lost, and one ended at
    1.9%. The clean (unperturbed) rate from 9.9M is therefore 1 of 24; the other two native
    successes came in E1-introduction runs, which perturb the community, so they are not
    equivalent trials (3 of 88 counts them all). From 10.0M: 0 of 24 clean, 0 of 88 in all. The
    clean success still settles the possibility question: a small world can make the processive
    transition spontaneously, without a big world or an introduction (from one evolved state). So the original loss
    at 10.02M was the usual outcome, not bad luck: once this lineage is in this state it is
    almost always lost, while a fresh start from 9.9M occasionally succeeds.
  - **The winners are short.** All three native successes are 26-32-symbol members of the
    S21 family with one copy block (`...>BO^$O%H...>C$===?...>`), matching the second half of
    S21 itself (`...Q>BO^$P%HIL>C$====?O>AMY$Z>...`). S21, the variant that arose in the
    original run and failed (95 symbols), carries two such blocks.
  - **The short form does better where S21 failed** (28 Oct; W26 =
    `Q>BO^$O%H>C$===?O>A$$T$TQ>`, the 26-symbol replay winner, introduced exactly as above at
    the same 16 locations with the same random streams, so each W26 run is paired with an S21
    and an E1 run). Established, N = 1 / N = 2 of 16: from 9.9M W26 1 / 2, from 10.0M W26 0 / 3,
    against S21 0 / 0 at both; 6 of 64 for W26 against 0 of 64 for S21 (one-sided Fisher
    p = 0.014) and 6 of 64 introduced-E1 successes. All six W26 populations are W26-like
    (26-29 symbols, e.g. `Q>BAOW^Q$O%H>C$===?O>A$$=$TQ>`). The successes are at different
    places: only one location (9.9M, pair, rep 5) worked for both W26 and E1. W26 founders are
    treated much like S21's: first bound as the passive partner by a resident in 14-15 of 16
    singles, and 6-12 per group altered. So the late community uses newcomers as templates
    whatever they are; what differs is that W26 sometimes comes through it and S21 never did.
    This is a genotype difference; it does not yet say whether length, the second copy block or
    susceptibility to resident reactions is the cause (and "inert" padding would not isolate
    length cleanly, since every symbol also affects binding and search labels). Still rare (about 1 in 10), so the
    native transition there is a long-odds event either way.
  - **How founders failed** (all 512 founders; `natural_agg.txt`). Few decayed without ever
    binding (0-3 per group). A common end is the first reaction: the founder binds a resident
    and both die in that reaction with no product, since under the old pair-decay rule a pair
    decays together (up to 16 of 32 founders in a group, e.g. S21 pairs in seed21021 at 2.5M).
    Others react repeatedly without a placed descendant (14 of 16 S21 singletons in seed36036),
    or have descendants that then die out; lineages that went extinct did so after a median of
    ~0.9-12k steps (by group). Products lost for lack of room are common (the worlds are full).
  - Checked by hand: the fates of the 15 founders of seed36036 / E1 / N = 1 read from their
    traces match the scorer's counts. A traced lineage being alive is not the same as a
    processive population: one S21 lineage in seed36036 spread to the trace limit while no
    processive species remained (its descendants had lost processivity).
  - Caveats: one 200k-step window per run, an endpoint test, 16 locations per recipient, one
    replicate set; the trace's exact-copy flag in these runs compared current strands, and the
    scorer recomputes exact copies from species ids (the simulator now does the same).
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
- **Reading, provisionally.** Processive variants do arise in some small worlds and (in all
  three cases) vanish. One hypothesis is that the barrier is getting a foothold rather than
  persisting: the transplants show that *selected big-world genotypes*, supplied as 20 adjacent
  molecules in a cleared 5x5 patch, spread in most evolved small communities. But that protocol
  supplies partners and space, and the naturally arising small-world variants may be poorer
  genotypes or a poor fit to their community, so the transplants alone don't say why they fail.
  The introduction-size results above now point to starting number (encounter) as a major
  barrier, with genotype and community as further ones. The two small cohorts (20 fresh runs; 18
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
- **Big replicates: final (29 Oct; `scan_big/`, `scan_sample/`).** None of the three 400x300
  classic replicates (seeds 19019, 20020, 21021) produced processive replicators by 10M steps:
  none in the abundance scan (species >= 20 molecules) or in random samples of 1,000 molecules,
  every 500k from 0.5M to 10M (20 samples each). The original run had detection at 3.0M,
  1.1% at 4.0M, 30% at 4.5M and ~80% from 6M. Nor is it a matter of activity: by 10M the
  replicates had made 5.0-6.3 x 10^9 attempted copies, 4-5 times the original's 1.2 x 10^9 at
  its establishment (by 3.0M they had 6.0-10.8 x 10^8 against its 6.3 x 10^8). So the original
  transition was **not reproduced in three runs**. That alone does not make it rare: even at a
  50% chance per run, three misses happen 12.5% of the time; but together with the 0 of 38
  small worlds and the 1-in-24 clean native success from the one small-world state that came
  closest (seed21021 at 9.9M), the evidence now points to a chancy event rather than a routine one at
  this scale. The replicates' histories differ from each other: seed20020 made far fewer new
  species (2.3 x 10^8 ids by 10M, against 5.0 and 5.6 x 10^8) at similar copying, and its
  cover fell from 0.77 to 0.15 over 3-4.5M, as the original's did during its crisis, without a
  transition following.
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

## 15. Disturbance: does clustered disturbance help innovations establish? (Oct 29)

**Question.** Can a genotype-blind ecological change make spontaneous innovations establish
more often? Hypothesis: occasional local openings let new reproductive strategies expand
alongside compatible neighbours before resident interactions eliminate them.

**Design** (`experiments/disturb/`). From seed21021 at 9.9M (where 1 of 24 clean replays made
the processive transition), 40 reseeded 200k-step runs per treatment, chemistry unchanged,
nothing introduced: *base*; *clustered*, every 1,000 steps remove every molecule in one random
7x7 patch; *dispersed*, every 1,000 steps remove whatever is in 49 random cells. Disturbance for
the first 150k steps, then off. Rep k shares its random stream across treatments and its
disturbance seed across the two disturbed ones. Scored with `--sample 1000` every 10k;
established = >= 10% processive for >= 3 consecutive samples. Established populations were
continued to 400k without disturbance. One schedule, chosen in advance; no sweep.

| treatment | runs | any processive seen | detected (>= 20 mol.) | established | extinct | mean molecules | placed copies (median) | removed by disturbance |
|---|---|---|---|---|---|---|---|---|
| base | 40 | 5 | 3 | 2 | 0 | 11,374 | 1.21 x 10^6 | 0 |
| clustered | 40 | 5 | 2 | 1 | 0 | 11,394 | 1.22 x 10^6 | 6,718 |
| dispersed | 40 | 6 | 2 | 0 | 0 | 11,357 | 1.21 x 10^6 | 6,703 |

(Fisher, established: clustered vs base p = 1.0; dispersed vs base p = 0.49; clustered vs
dispersed p = 1.0. `results/disturb_agg.txt`.)

- **Decision: neither disturbance improved establishment, so this hypothesis is closed** (by
  the rule set beforehand). Disturbance changed neither how often processive variants appeared
  (5-6 runs of 40 in every treatment) nor how often they established.
- **Turnover did not stand in for potential:** occupancy and placed copies were the same in all
  treatments; nothing went extinct.
- **The baseline itself is a better estimate now:** 2 of 40 clean replays from 9.9M
  established, with the earlier 1 of 24, 3 of 64 (~5%).
- **Lasting, not blooms:** all three established populations (base rep3, base rep18,
  clustered rep35) were still there at 400k with disturbance off: 50%, 43% and 19% of
  molecules (the last falling from a 48% peak). All are short S21-family replicators again
  (29-38 symbols, e.g. `$=P>BO?^$O%HK?>C$====?M>MZ$I>`). One more base run (rep6) ended at
  10.8% without meeting the three-sample rule.
- **Scope:** this tested one modest dose. 49 cells every 1,000 steps removed ~6,700 molecules
  over 150k steps, under 1% of the run's deaths (it places ~1.2 million products). "No effect"
  holds for this dose; a heavier dose would be a different question, and the rule set
  beforehand was to close rather than sweep.

## 16. What first differs between a lineage that establishes and one that disappears? (Oct 30)

From existing histories only (`experiments/whyfail/`): all 168 runs from the two seed21021 states
(120 disturbance runs and 24 replays from 9.9M, 24 replays from 10.0M). `lineages.py` probes every
species that reached 3 molecules and groups processive species into lineages by parent links (a
processive species whose parents are not processive founds a new lineage). `compare.py` compares
lineages over the 5,000 steps after each first reached 5 molecules, from `copies.tsv` (who
copied whom) and `grid.snap` (species map every 500 steps).

- **Origins.** 118 lineages got a foothold (>= 5 molecules) from 9.9M and 53 from 10.0M. Each is
  a separate origin: processive forms arise again and again from non-processive residents of
  the same family, which carry the copy core (`...BO^$O%HK?>C$====?M>MZ$I>` and relatives) behind
  long repetitive prefixes; a shortening makes them processive. No founder sequence occurred
  twice in 118, so the histories cannot compare one genotype in two places.
- **Outcomes.** From 9.9M, 3 lineages established, 3 more grew large (peak >= 500 molecules,
  e.g. disturb base rep18, which the molecule-sample scan counted as established; the lineage
  count misses its rarest variants) and 112 failed. From 10.0M (the 95-symbol S21 lineage
  present) all 53 failed.
- **Space did not differ:** 3-4% of neighbouring cells were empty for every group, and
  successful lineages were not more clustered.
- **Not obviously a better genotype:** failed founders were similar in length (median 40 vs 30)
  and made more self-copies in the probe (median 54 vs 35).
- **Later mutation is not the first difference:** in all three established lineages the form
  that came to dominate arose 42-48k steps after the foothold, but the lineages had already
  started growing (x2.0-2.4 in the first 2,000 steps; failures x0.8).
- **The first difference is who copies them.** In the first 5,000 steps successful and large
  lineages had more of their reactions with kin (0.48-0.65 vs 0.19), and when residents ran
  their code on a member, the residents were mostly relatives (non-processive residents carrying
  the founder's copy core), and 70% of those copies were new members. In the 15 most-copied failed
  lineages most such copying was by unrelated residents (2,043 of 2,683 attempts), and 15% gave members
  (`results/same_geno.txt`).
- **The same signal in the introduction runs** (`intro_neigh.py`, `intro_motif.py`; recipient
  cells from the 9.9M/10.0M checkpoints). W26's 6 successes had a larger share of neighbouring
  cells holding S21-family residents (carrying `>C$==`): 0.40-0.88 against a mean of 0.52 for
  its 58 failures (one-sided Mann-Whitney p = 0.008; pairs only p = 0.039). Among W26 pairs,
  4 of 13 established where >= 70% of neighbours were family, 1 of 10 where < 50%. **E1 at the
  same cells shows the opposite** (3 of 10 at < 50%, 1 of 13 at >= 70%; p = 0.97), so the
  effect is not that some places are generally better: it depends on the newcomer's fit to its
  neighbours. Close copies of W26's own motif did not help (none next to any success).
- Caveats: few successes (6 natural, 6 introduced); the family definition (`>C$==`) was chosen
  after looking, out of three tried; a lineage counts only species that reached 3 molecules.

**Prediction.** A new S21-family processive replicator establishes much more often when it
arrives among residents of its own broad family, which copy it faithfully into new members,
than among unrelated residents, which use it as a template and alter it; for an unrelated
replicator (E1) the same neighbourhoods should give no advantage, or the reverse.

## 17. Reaction-level equivalence with the C++ code (Nov 15)

Details: `experiments/equivalence/README.md`.

- **Test:** 23,040 reactions between the paper's own evolved molecules (top 12 species, 8 runs,
  every 100k steps, every ordered pair), about 25 million steps. Each reaction was run in the C++
  spatial code and in smol, with matching random draws, and compared step by step.
- **Current dev:** every reaction matches exactly, with mutation off and at 1%: binding
  probability, choice of active molecule, every pointer at every step, products and outcome.
  The exceptions are 17-18 reactions where a copy fills the 2,000-symbol buffer and the C++
  code writes past its end.
- **The paper's version (0.2.3.4):** the same, except 0.7% of reactions. These all hit the
  uninitialised alignment score at a `?` with the read pointer at the end of a string (fixed in
  #50), so the paper's binary was build-dependent there.
- **One rule difference found:** in the C++ code a copy-deletion also skips the next
  instruction. Now a switch, `COMPAT_DELETION`, which `pub_to_ckpt.py` turns on.
- **Not yet tested:** world-level rules (update order, partner choice, placement), which have
  only been compared statistically.

## 18. Do parasites hold the population down? (Nov 16)

The paper says evolved resistance to parasitism "explains how the population size is larger
than seen before". The premise is that parasites hold the population down. This tests it
directly.

**Method** (`experiments/parasites/`):
- Take the 8 published runs at 1M, converted with `pub_to_ckpt.py` (old rules).
- Classify every species present by what it did in the next 20k steps:
  - R: made exact copies.
  - P (parasite): never copied anything but was copied exactly.
  - N: neither.
- Parasites were 26-62% of molecules; N was 5-7% (23% in run 1).
- Treatments: remove all P; remove all N; and, as controls, remove the same number of random
  molecules (cullP, cullN), since freeing space alone lets the population regrow.
- 3 reseeds each, 100k steps, mutation on.

**Result:**
- **Removing parasites releases the population, but only for a while.**

  | Steps after the edit | 2k | 10k | 20k | 50k | 100k |
  |---|---|---|---|---|---|
  | noP minus cullP, median of 8 runs | +1,631 | +3,087 | +1,325 | -6 | +323 |
  | Runs with a positive difference | 8/8 | 8/8 | 7/8 | 4/8 | 6/8 |

  At 10k, run 3 has 11,818 molecules against 3,912 unedited.
- **Removing N does nothing** (median -430 at 10k).

**Reading:** parasites are what holds the population below a lone replicator's capacity. The
paper's premise survives a direct test. But a one-off removal is undone within about 50k steps,
presumably as parasites re-arise by mutation (not yet measured). So a lasting rise needs
something that keeps parasite load down, which is consistent with their resistance story but
does not by itself show resistance is the mechanism. Caveats:
- Only 3 seeds per treatment; the populations swing widely.
- The classification is behavioural, over smol's first 20k steps.
- In our own runs the population rose with little fall in our parasitism measure (section 12).

## 19. Time-shift swaps of replicators and parasites (Nov 16)

Details: `experiments/parasites/` (`swap.py`, `swap.sh`, `pswap.sh`, `results_swap.txt`).

**Method:**
- For each published run, take the 1M and 2M worlds (converted, old rules), classified as in
  section 18.
- Replace one class (replicators R, or parasites P) with molecules drawn by abundance from the
  same run at the other time. The control replaces that class with fresh draws from its own
  time, so the disturbance is the same.
- 8 runs x 4 swaps x 5 seeds, 100k steps, mutation on.
- Parasite share = molecules whose sequence was classed P at 1M or 2M.

**Results** (median differences from the matched control over 8 runs):

| Swap | +10k population | +100k population | +10k parasite share | +100k parasite share |
|---|---|---|---|---|
| 2M replicators into 1M worlds | +2,909 (8/8 up) | +1,745 (8/8) | -0.20 | -0.08 |
| 1M replicators into 2M worlds | -2,564 (0/8 up) | -713 (0/8) | -0.01 | -0.02 |
| 2M parasites into 1M worlds | +1,294 (5/8 up) | -117 (3/8) | -0.32 | -0.03 |
| 1M parasites into 2M worlds | +1,019 (8/8 up) | +491 (6/8) | -0.31 | -0.02 |

**Reading:**
- **Replicators improve, and the improvement lasts.** Later replicators raise the population
  of an earlier world for the whole 100k steps. Earlier replicators lower the population of a
  later world for the whole 100k steps.
- **Parasites do not improve.**
  - Moved in either direction, they do worse: their share falls by about 0.3 against the
    replicators of another time.
  - The effect fades as the share returns, presumably through fresh parasites mutated from the
    resident replicators.
  - This fits the paper's own claim that new parasites arise from mutated replicators rather
    than from evolving parasite lineages: parasites are matched to the replicators of their own
    time, not getting better.
- **The replicators' advantage is not mainly keeping parasites out.**
  - Against the older 1M parasites, the later replicators do cut the parasite share (-0.20),
    which is resistance to foreign, older parasites.
  - But in their own 2M worlds, swapping back to 1M replicators lowers the population with no
    change in parasite share.
  - Parasite share at 2M is still 31-40% in most runs (7% in run 4).
  - So most of the later replicators' advantage is tolerance (losing less to the parasites they
    have) or productivity, not resistance. "Evolved resistance explains the larger population"
    overstates it.

**Caveats:**
- 5 seeds per swap; the populations swing widely.
- The classification is behavioural, over smol's first 20k steps from each snapshot.
- These are smol runs under the old rules, not the C++ code.
- A swap breaks all reactions in the swapped class; the controls share that.

## 20. Do replicators learn to reject their own parasite mutants? (Nov 16)

Details: `experiments/parasites/mutassay.py`, results `mutassay_paper.txt`.

**Method:**
- In the 8 published runs, every 100k steps, take the top replicator: the most abundant
  species that copies itself in the probe.
- Make all of its single-symbol substitutions and deletions (about 800-2,000 mutants).
- React the replicator with each mutant, mutation off.
- Mutants that can no longer copy themselves are its would-be parasites (55-75% of mutants).

**Measures:**
- Binding: P(bind) with would-be parasites against functional mutants.
- Copying: binding-weighted share of each kind that gets copied in an encounter.
- disc = 1 - copyP/copyF.

**Result:**
- **No exclusion at binding:** P(bind) is the same for would-be parasites and functional mutants
  (about 0.62-0.65) throughout, in every run.
- **Discrimination is low and mostly flat:**
  - Median disc is 0.24-0.31 from 200k to 2M, with no trend.
  - Would-be parasites are copied in about 65-73% of encounters against about 93-96% for
    functional kin.
  - The baseline about 0.25 is present from 200k; it is presumably structural (mutations that
    break copying also break what the checks read), not an evolved improvement.
- **Two runs do improve in the second million:**
  - Run 2: disc rises from about 0.2-0.3 to 0.5-0.7 from 1.1M.
  - Run 4: disc rises from about 0.25 to 0.55-0.65 from 1.3M.
  - Run 4 is also the run whose parasite share fell most (26% at 1M to 7% at 2M). Run 2's share
    did not fall (35% to 32%).
- **More checks are not better checks:** the number of '?' tests in the top replicator rises
  from 2 to 3-5 in runs 4, 7 and 8, but discrimination rises only in run 4. In run 7 it falls
  as the checks accumulate.

**Reading:** the paper's claim holds in a minority of runs, not as a general trend:
- Explicit discrimination against near-self parasites does evolve in some runs (clearly in run 4,
  where it matches a large fall in parasites).
- Most runs keep copying their own parasite mutants about two-thirds as readily as their kin for
  the whole 2M steps.
- "Increasingly complex" mechanisms (more '?') appear without becoming more discriminating.

Caveats:
- One replicator per snapshot.
- Single-step mutants only.
- smol's reaction probe; the reactions match the C++ code (section 17).

## 21. Tolerance fingerprints and "evolution evolves", 1M vs 2M (Nov 16)

Details: `experiments/parasites/mechanisms.py`, results `mechanisms_paper.txt`. For each
published run at 1M and 2M: the top replicator S against the world's 10 commonest parasites in
the probe, the world snapshot, and the 100k-step control runs that follow.

| Measure | 1M | 2M | Runs that rose |
|---|---|---|---|
| S self-copy time (steps per copy) | 178 | 153 | 1/6 (it fell in 5) |
| Parasite time per copy / kin time per copy | 0.73 | 0.77 | 4/5 |
| Copies per self-binding | 1 | 1 | 0/6 |
| Two-way copying in S-parasite encounters | 0 | 0 | 0/5 |
| Replicator clustering (neighbours, vs random) | 1.34 | 1.24 | 2/8 |
| Copy-operator executions per product symbol | 2.01 | 2.03 | 3/6 |
| Non-exact products that are >2 edits from both reactants | 0.80 | 0.86 | 5/8 |

**Tolerance:**
- Two of the candidate mechanisms are ruled out: no run makes more than one copy per binding,
  and replicators are never copied back by parasites.
- Clustering, if anything, falls. Caveat: the index's ceiling falls as the replicator share
  rises.
- Parasites keep most of their speed advantage (copied in about 75% of a kin copy's time at both
  times).
- What does change is the replicators' own speed: self-copy time falls by about 14% in 5 of 6
  runs. So the clearest difference is simply faster replicators. That is efficiency, which the
  paper set aside on the basis of Figure 5 (well-mixed, no parasites); with parasites present, a
  modest speed gain may matter more than Figure 5 suggests.
- Not yet measured: parasites trapped in long unproductive reactions, which the paper's t3
  description suggests.

**Evolution evolves:**
- Effective mutation exposure is already doubled at 1M in every run with a top replicator
  (self-scan, 2.0 copy operations per product symbol) and does not rise further.
- Of new strings made by reactions, about 80-86% differ from both reactants by more than 2
  edits: emergent rearrangements (partial copies, duplications, mis-cleaves) are by far the main
  source of new sequences, not point mutations, as the paper says. Their share rises only
  slightly (5/8 runs).
- So the phenomena are real, but they are established early. There is no evidence here of
  evolvability itself increasing over the second million.

## 22. Runs with parasites suppressed from the start (Nov 16)

Details: `experiments/parasites/suppress.py`, results `suppress_*.txt`.

**Method:**
- 30 runs from the seed replicators (the paper's box config, old rules, COMPAT_DELETION), 2M
  steps each, in 5k-step chunks.
- After each chunk: P arm removes every species that was copied exactly but never copied
  anything and is at least 5k steps old. The cull arm removes the same number at random. The
  none arm only restarts.
- 10 seeds each.

**Survival to 2M:**
- none: 4/10.
- P: 8/10 (deaths at 0.69M and 1.29M).
- cull: 3/10, plus one still running.

**How much parasitism the P arm removes:**
- Established parasite share at chunk ends after 1M (medians): none 0.30, P 0.08 (regrown
  within 5k steps), cull 0.16.
- Parasites younger than 5k steps are never removed or counted. So the P arm reduces parasitism
  a lot but does not remove it, and fresh mutant parasites, the kind that matter most (sections
  19-20), are still present.

**Results** (top replicator, medians across surviving runs):
- **Shortening:** length falls from about 60 to 25-30 by 300k in every arm, with or without
  parasites.
- **Self-scan:** 2 copy operations per product symbol appears without parasites. At 2M, 6 of 8
  P-arm runs are self-scanning (`suppress_expo.txt`), as are the none runs. Self-scan is not
  driven by persistent parasites.
- **Checks and discrimination:**
  - The '?' count is 2 in every arm from 200-300k. These are largely the copy loop's own test,
    so the count is a weak measure of discrimination.
  - Mutant discrimination is about 0.15-0.25 in the P arm and 0.2-0.35 in the none arm:
    slightly higher with parasites late on.
- **Second million:**
  - Without established parasites, replicators keep getting shorter (about 20) and faster
    (self-reaction about 120-140 steps).
  - With parasites (only 3 surviving none runs), they stay longer (30-33), slower (about
    210-250 steps), and gain a third or fourth '?' in some runs.

**Reading:**
- The main changes the paper attributes to parasitism (shortening, self-scan, the basic check
  structure) evolve just the same when established parasites are kept below about 8%.
- Parasitism is associated only with some late extra length, slowness and checks, in 3 runs.
- Without parasites the system ends simpler and faster, and survives far more often.
- Caveat: fresh parasites remain in the P arm. A stronger suppression is needed before
  concluding that parasitism plays no part in self-scan.

## 23. PARASITEFREE: no string is placed unless it can copy itself (Nov 16)

Details: `experiments/parasites/pfree_*.txt`.

**Method:**
- 20 runs from the seed replicators (box config, old rules, COMPAT_DELETION), 2M steps:
  10 with PARASITEFREE 1 (pf), 10 identical controls (ctl).
- PARASITEFREE blocks about 340-870 products per 1000 steps once the grid is full.

**Survival and population:**
- pf: 10/10 survive and hold a full grid (about 12,200 molecules) from 100k on.
- ctl: 5/10 survive (deaths at 73k-306k; two survivors still running at writing), 6-11k molecules.

**Results** (top replicator, medians):

| | pf | ctl |
|---|---|---|
| Length at 200k | 60 | 31 |
| Length at 1M | 42 | 28 |
| Length at 2M | 38 | 28 |
| '?' tests | 1 throughout | 2 from 200k |
| Self-reaction steps | 215 -> 105-135, falling steadily | 170-230 |
| Self-scan (copy ops per product symbol) | 1.00 in 28 of 29 snapshots | 2.0 in every survivor from 300k |

- **Self-scan:** the one pf exception is run 210 at 2M (2.00).
- **Shape of the replicators:**
  - pf replicators keep much of the seed's long non-coding head and shorten only slowly.
  - ctl replicators lose it by 200k and duplicate the copy-loop motif (`$=?>...^B$=?>`), which is
    the self-scan.
- The mutant-discrimination numbers are not comparable: in pf the top replicator rarely copies
  even functional mutants when bound to them, because it is usually the passive partner.

**Reading:**
- With parasites impossible, the changes the paper attributes to parasitism mostly do not
  happen: the fast collapse to short replicators, self-scan (and with it the doubled mutation
  exposure), and the second loop test.
- The replicators instead get steadily faster while staying long.
- Section 22's tentative conclusion (self-scan evolves without parasites) was wrong: there, fresh
  parasites were never removed.

**Confound:** PARASITEFREE also blocks every product that cannot copy itself, including the
non-viable intermediates of the mutational cascades the paper shows leading to the first
self-scanning replicator (their figure 6). Self-scan may be absent because parasites are absent,
or because the route to it is blocked. The test that would separate the two: let such products
exist but never be copied.

## 24. Competition: self-scanners against their non-scanning contemporaries (Nov 16)

Details: `experiments/parasites/compete.py`, `compete_results.txt`, `compete_pairs.tsv`,
`compete_probe.txt`.

**Method:**
- The chemistry is not changed; only the starting mix is.
- For 6 published runs, take the first 20k snapshot where a self-scanning replicator (D) is
  among the top replicators, with the commonest non-scanning replicator (A) of the same
  snapshot.
- Take up to 5 of that snapshot's parasites (non-self-copiers that A or D copies).
- Fresh 125x100 world: 1000 A + 1000 D at random cells, with or without 1000 of the parasites;
  mutation off or on.
- 10 seeds per condition, 30k steps.

**Result** (share of D among A + D at 30k, median of 10 seeds):

| Run | mutation off, no parasites | mutation off, parasites | mutation on, no parasites | mutation on, parasites |
|---|---|---|---|---|
| 1 | 0.72 | 1.00 | 0.43 | 1.00 |
| 2 | 0.42 | 1.00 | 0.87 | 1.00 |
| 3 | 0.25 | 1.00 | 0.76 | 1.00 |
| 7 | 0.74 | 1.00 | 0.79 | 0.95 |
| 4 | 0.82 | 0.39 | 0.68 | 0.54 |
| 5 | 0.99 | 0.00 | 0.99 | 0.61 |

- **In 4 of 6 runs (1, 2, 3, 7) self-scan is an anti-parasite adaptation:**
  - With parasites, D takes over completely within 10-20k steps, in every seed.
  - Without parasites and with mutation off, it is neutral or worse (0.25-0.74). In run 3 it
    loses.
  - With mutation on and no added parasites, D gains in runs 2 and 3 (0.87, 0.76), as parasites
    arise from both competitors.
  - In the probe, those parasites are copied readily by A and hardly at all by D (exposure, i.e.
    P(bind) x copied: A 0.60-0.91, D 0-0.12), and D lets go of them almost at once (5-26 steps
    bound against 97-163 for A).
- **In runs 4 and 5 it is reversed:** D wins without parasites and loses with them (run 5: 0.00
  in 8 of 8 surviving seeds). D still resists those parasites better than A in the probe
  (exposure 0.16-0.25 against 0.45-0.57), and neither parasite set copies A, so the pairwise probe
  does not explain the reversal. It may be a spatial or community effect. Unexplained for now.

**Reading:**
- This is the strongest support yet for the paper's claim, and it comes without changing the
  rules. Where self-scan first appears, it is usually favoured because of parasites: it lets a
  replicator refuse and release the parasites its competitor keeps copying.
- The advantage is not universal. In 2 of 6 runs the same contrast favours the other replicator.
- Together with sections 19-23: self-scan and the early restructuring are parasite-driven, while
  later improvement is mostly tolerance and the parasites do not themselves improve.

## 25. Are parasites like defective interfering particles? Genealogy and regeneration (Nov 17)

Details: `experiments/parasites/dip.py`, `dip_regen_agg.py`, results `dip_genealogy_ctl.txt`,
`dip_regen_population.txt`, `dip_regen_genealogy.txt`.

**Genealogy** (5 surviving normal runs from section 23, 2M steps each):
- Method:
  - A successful parasite is a species copied exactly at least 20 times that never made a copy.
  - Its template parent is the nearer (edit distance) of its two parent species in species.tsv.
  - Its lineage is followed back through parasite parents to the first non-parasite ancestor.
- **Only 42-58% of parasite species are made directly from a non-parasite** (mostly a
  replicator); the rest descend from other parasites, up to hundreds of generations deep.
- Weighted by how often they are copied, **only 28-44% of parasitism comes from parasites made
  directly from replicators.**
- The largest parasite clades have 200-520 species. Parasite lineages last a median 150k-440k
  steps, against 23k-85k for replicator species.
- So the paper's "new parasite species continually arise from mutated replicators, rather than
  from evolving parasite lineages" does not hold in these runs: parasite lineages are real and
  long-lived.
- **The lineages do not get better.** P(bind) x copied by the current top replicator, for
  parasites made directly from a replicator against those descended from parasites: 0.28 vs 0.32,
  0.26 vs 0.26, 0.19 vs 0.12, 0.39 vs 0.46, 0.40 vs 0.39. No consistent difference.

**Regeneration** (3 parasite-free worlds at 2M, continued 200k steps, 5 seeds each):

| Steps after the switch | +0 | +10k | +20k | +100k | +200k |
|---|---|---|---|---|---|
| Parasites allowed, mutation on (population) | 12,235 | 11,968 | 9,817 | 8,464 | 10,362 |
| Parasites allowed, mutation off | 12,276 | | | | 12,213 |
| Still parasite-free | 12,285 | | | | 12,226 |

- With mutation on, the first successful parasite species appears within 1,000-2,000 steps in
  all 15 runs. 2 of 15 go extinct.
- With mutation off or parasites still blocked, nothing changes.
- Within 200k steps the new parasites already form lineages: 41-66% of species, and 12-67% of
  parasitism (median about 40%), descend from parasites. Exploitation by lineage depth is again
  about equal (deeper slightly higher in 9 of 15 runs, lower in 3; mostly by 0.01-0.1).

**Reading:**
- Like defective interfering particles, parasites arise de novo from the replicators almost
  immediately once mutation can make them, they need the replicators, and they suppress the
  population.
- Unlike the "regenerated each time" picture (and the paper's claim), about half of them
  propagate as parasite lineages that persist for hundreds of thousands of steps.
- Those lineages drift rather than escalate: no clear improvement in exploitation along them.
- The best description is persistent, non-escalating parasite lineages, continually topped up
  from the replicators.

## 26. Why self-scan lost in competition runs 4 and 5 (parasites adapt to the scanner; decoys protect A)

The competition (section 24) had D (self-scanner) beat A with parasites present in runs 1, 2, 3 and 7,
but lose in runs 4 and 5. Scripts and results in `experiments/parasites/`.

- **Run 5:** one fast D-specialist, `$$YHOB`. A and D bind it equally (0.80); A lets go after 8 steps
  without copying it, D copies it in 37 steps (its own copy takes 211). In the first 1,000 steps D made
  5,070 copies of it and 841 of itself.
- **Snapshot pools (`ratio.py`):** same A and D against 31 pools taken from snapshots with different
  D:A ratios (10 seeds each, mutation off). Pools from snapshots where D's family had taken over beat D
  in 11 of 13 cases; pools from before D existed beat D in 2 of 6. The original run 4 and 5 pools came
  from snapshots where D already outnumbered A by 5.7x and 4.2x (other runs 0.7-2.1x).
- **Evolved pools (`evolve.py`, causal):** parasites evolved for 150k steps in worlds of D alone or
  A alone (mutation on, 3 seeds per run), then used in the competition. D-evolved pools beat D in
  35 of 35 (median D share 0.00-0.01); A-evolved pools mostly favour D (median 0.95-0.99), beating it in
  7 of 28. 4 of 18 A-only worlds went extinct, 0 of 18 D-only.
- **Prediction (`predict.py`):** each host reacted with each parasite over 10 seeds (`smol-react --seed`,
  added for this). "The host with the fastest-growing parasite (P(bind) x success x self-copy steps /
  parasite copy steps) loses" gets the winner right in 85 of 93 pools.
- **Run 4 (`knockout.py`, `dropone.py`):** D's specialist `>$$OZHOB` matters, but the decisive members
  are *decoys*: long relatives of A that A copies but that cannot copy themselves. They bind A's
  parasites for 160-360 steps without copying them and release D's parasite within 10-45 steps (one
  also copies it). Without the decoy, the original run 4 pool goes from D share 0.39 to 1.00; without
  both decoys, the 60k pool goes from 0.00 to 1.00 (10/10 seeds). Removing members that copy other
  parasites ("alternative hosts") did not explain it.
- **Reading:** self-scan's advantage is real but temporary: it beats parasites adapted to other
  replicators and loses once parasites adapt to it, an arms race rather than steady improvement.
  Decoys are a second tolerance mechanism (besides faster self-copy); how common they are is untested.

## 27. The stepping-stone confound: PARASITEFREE 2 and the genealogy of the first self-scanner (Nov 17)

Section 23's confound: PARASITEFREE 1 blocks every product that cannot copy itself, which also
blocks the non-viable intermediates on the way to self-scan. (The test proposed there, "let such
products exist but never be copied", cannot work: a molecule that cannot copy itself only has
descendants if it is copied.) Details: `experiments/parasites/scanpath.py`, `scanpath_ctl.txt`,
`scanpath_pf2.txt`, `pfree_complexity_pf2.txt`, `pfree_expo_pf2.txt`.

**Genealogy (controls):** the first self-scanner (self-copier, copy ops per product symbol >= 1.5,
>= 20 molecules) traced back through species.tsv to the nearest non-scanning self-copier. In all 6
control runs where self-scan appeared, the path runs through 1-8 species that cannot copy
themselves, typically a long mutant with a duplicated copy loop that replicators copied (up to
hundreds of times); the self-scanner is a truncated copy of it.

**PARASITEFREE 2:** a product that cannot copy itself is placed unless it is an exact copy of
either reactant. Non-copying mutants appear (and are re-made by mutation) but never multiply.
Same 10 seeds and set-up as section 23, 2M steps.

| | ctl | pf (1) | pf2 |
|---|---|---|---|
| Survive to 2M | 5/10 | 10/10 | 9/10 (one died at 1.97M) |
| Population, median at 1M | ~10,100 | ~12,100 | ~12,200 |
| Self-scan appears (>= 20 molecules) | 6/10, 21k-182k | 1 snapshot, at 2M | 10/10, 49k-288k |
| Top replicator self-scans, 300k / 1M / 2M | every survivor | 0 / 0 / 1 | 5/10, 3/10, 5/9 |
| Top replicator length at 2M | 28 | 38 | 17 |
| '?' tests | 2 | 1 | 2 |
| Self-reaction steps at 2M | ~180 | ~120 | ~92 |

- In pf2 the paths to the first self-scanner have 1-18 non-copying intermediates, none of them ever
  copied exactly: their numbers come from recurrent mutation.
- pf2 blocks 400-1,500 exact copies per 1,000 steps: fresh non-copying mutants still waste copier
  time, but the population stays at the parasite-free level.

**Reading:** PARASITEFREE 1's result was mostly the blocked route. Shortening, the second '?' and
self-scan come from the mutational route through non-copying intermediates, not from parasite
pressure. Multiplying parasites supply selection: with them, self-scan wins and stays (every
control survivor; competition, section 24); without them it comes and goes (about half of pf2
snapshots), as in the competition without parasites, where D is roughly neutral.

## 28. How later replicators cope: stalling on parasites, not decoys (Nov 17)

Details: `experiments/parasites/decoy.py`, `decoy_survey.txt`, `swapprobe.py`,
`swapprobe_results.txt`.

**Method:** the 8 published runs at 100k and every 200k to 2M (86 snapshots). The 30 commonest
species of each snapshot are reacted pairwise (smol-react --compat, one trial). Replicator: makes
an exact copy of its partner or itself in some pair. Parasite: never does, but is copied exactly by
a replicator. (A self-copy test alone misses runs where replication is by partners, e.g. run 6
after 400k and run 1 after 1M.) Each parasite's expected bound time (parasite x partner abundance x
P(bind) x steps) is split into: copied; held >= 50 steps with no copy, by a replicator ("stalled")
or by a non-replicator ("decoy"); released sooner.

**Decoys (findings 26) are rare:** decoy molecules are 0% of the population in the median
snapshot and >= 2% in 10 of 86; decoy holding is a median 1-4% of parasite-bound time and does not
track population (Spearman -0.05).

**Stalling grows:**

| Period | Median pop | Copied | Stalled by replicators | Held by decoys |
|---|---|---|---|---|
| <= 600k (32 snapshots) | 5,754 | 0.65 | 0.08 | 0.04 |
| 0.8-1.2M (24) | 7,350 | 0.70 | 0.18 | 0.01 |
| 1.4-2M (30) | 9,462 | 0.30 | 0.59 | 0.02 |

- Stalls produce no product at all (not botched copies).
- 1M -> 2M: stalling rises and copying falls in 8 of 8 runs.
- Spearman with population: stalled +0.47, copied -0.44, parasite share of the top 30 -0.71.

**Swap (findings 20):** against the 2M parasites, the 1M top replicator copies (0.75-0.97) where
the 2M one holds them in runs 1, 3 and 5; in run 4 the 2M replicator does not bind them at all.
These are four of the five largest population drops when 1M replicators replaced 2M ones
(-6,929, -4,114, -3,309, -1,947 at +10k). Runs 2 and 8, where both hold, had the two smallest drops
(-613, -326). Run 7 (both hold, drop -3,181) is the exception; run 6 (both copy, -1,720) is middling.

**Reading:** what findings 20 called tolerance is mostly resistance by stalling on the established
parasites, not discrimination of own mutants (fresh mutants are still copied 65-73% of the time,
findings 21), which is why parasite share stays up: new parasites keep arising, and parasite
lineages adapt to the new replicators (findings 26).

**Caveats:** one trial per pair; only the 30 commonest species; the swap probe uses only the top
replicator. A direct intervention on stalling has not been done.

## 29. World-level statistical check against the paper, with COMPAT_DELETION (Nov 17)

Details: `experiments/published/cmp3.py`, `cmp3_results.txt`, `assay_pubdel.txt`,
`pooled_assay_pubdel.txt`. Runs on aida in `~/smol/runs/pubdel/`.

**Method:** the paper's three configs (box, random pairs, slot) with the old rules and
COMPAT_DELETION, 8 seeds each, 2M steps. Seeds 1101-5505 have the same starting worlds as the
October rerun (section 12); 6606-8808 are new (the random-pairs ones reuse seed 1101's layout).
`tools/assay.py` every 100k, compared over runs that survived to 2M, averaged over 1.0-1.2M and
1.9-2.0M; Mann-Whitney U against the paper's 8 runs.

| Survivors | paper (8/20) | Oct rerun, old rules (5/15) | rerun + COMPAT_DELETION (7/24) |
|---|---|---|---|
| Population, 1.9-2.0M | 10,130 | 9,690 | 9,999 |
| Parasitic share, 1.9-2.0M | 0.47 | 0.50 | 0.72 |
| Mutual share, 1.9-2.0M | 0.53 | 0.50 | 0.28 |
| Reaction time, 1.9-2.0M | 175 | 178 | 158 |
| Change in parasitic share, 1.0-1.2M -> 1.9-2.0M | -0.19 | -0.08 | -0.04 |

- Survival by config: box 3/8, random pairs 1/8, slot 3/8.
- No difference from the paper is significant (all p >= 0.06, 12 tests); ranges overlap widely
  (paper's parasitic share at 2M spans 0.17-0.89).
- Both smol reruns show a weaker second-million shift than the paper; COMPAT_DELETION does not
  close it (if anything the shift is smaller).
- Identical starting worlds diverge at once (chaos), so seeds 1101-5505 are not a paired test.

**Reading:** smol is statistically indistinguishable from the paper at n = 7-8 survivors, but a
modest difference in the second million cannot be excluded; it is not the deletion rule. Settling
it would need ~100 smol runs and fresh C++ runs (the paper's 8 runs are all there is). The
parasitism results (sections 18-28) do not depend on it: they use the paper's own worlds converted
to smol, or compare smol conditions with each other.

## 30. Dose-response: graded parasite amplification (Nov 18)

Details: `experiments/parasites/dose.sh`, `dose_summary.py`, `dose_results.txt`, `dose_*.txt`.

**Method:** PARASITEFREE 2 with `PFPROB` q: every product is placed except that an exact copy of a
non-self-copier is blocked with probability q. Mutants that cannot copy themselves always appear,
so the route to self-scan (section 27) is open at every level; only parasite amplification varies.
q = 0.25, 0.5, 0.75 run fresh; q = 0 is the section 23 controls, q = 1 the section 27 pf2 runs.
Same seeds 201-210 and set-up throughout, 2M steps. Measured share of exact copies let through:
0.75, 0.50, 0.25 (as set).

| q | 0 | 0.25 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|
| Survive to 2M | 5/10 | 6/10 | 7/10 | 8/10 | 9/10 |
| Self-scan appears (>= 20 molecules) | 6/10 | 8/10 | 8/10 | 10/10 | 10/10 |
| Median first appearance | 122k | 125k | 114k | 76k | 83k |
| Top replicator self-scans, per-run share of snapshots (median) | 1.00 (n=4) | 0.80 (7) | 0.50 (10) | 0.40 (10) | 0.63 (10) |
| Self-reaction steps at 2M (median) | 165 | 154 | 118 | 104 | 92 |
| Length at 2M | 24 | 23 | 20 | 20 | 17 |
| Population at 2M | 9,910 | 10,402 | 10,291 | 12,086 | 12,269 |

- Spearman rho(q, per-run self-scan share) = -0.39 over 41 runs with >= 3 scored snapshots
  (t = -2.62).
- The fall is mostly between q = 0 and 0.5; from 0.5 to 1 self-scan holds the top spot about half
  the time.
- Self-reaction time, length, survival and population all move steadily with q.

**Reading:** parasite amplification selects for self-scan in a graded way, and self-scan costs
copying speed; without the pressure it is roughly neutral (comes and goes). Together with section
27: the route makes self-scan, multiplying parasites make it win.

**Caveats:** at q = 0 only 4 runs live long enough to score (survivors may be biased towards
self-scan); one trend test on 41 runs; snapshots scored every 200k.

## 31. Time-shift matrix: parasites against host worlds from every other time (Nov 20)

Details: `experiments/parasites/timeshift.sh`, `timeshift_agg.py`, `timeshift_results.txt`. Runs on aida
in `~/smol/runs/parasites/timeshift/`.

**Method:** the 8 published runs at 200k, 400k, ..., 2M (10 times), converted and classified as
in section 18 (20k window). For every host time H and parasite time P of the same run, the H
world's parasites (class P) were replaced by parasites drawn by abundance from the P world
(`swap.py ... P`); H = P (own parasites redrawn) is the control. 3 seeds, 20k steps; 2,400 runs,
no extinctions. d = relative change in population at +10k against the control: d > 0 means the
foreign parasites are less harmful than the world's own.

| Lag (parasite time - host time) | -1.8M | -1.0M | -0.4M | -0.2M | +0.2M | +0.4M | +0.6M | +1.0M | +1.4M |
|---|---|---|---|---|---|---|---|---|---|
| median d | +0.25 | +0.12 | +0.11 | +0.08 | -0.01 | +0.04 | +0.13 | +0.36 | +0.59 |
| cells with d > 0 | 8/8 | 36/40 | 50/64 | 53/72 | 34/72 | 37/64 | 40/56 | 32/40 | 23/24 |

(+20k gives the same picture.)

- **Own-time parasites are the most harmful.** Parasites from the past are less harmful in 300 of
  360 cells (median d +0.11); from the future in 253 of 360 (+0.17). Only parasites from just
  after the host's time (+0.2M, +0.4M) are about as harmful as its own.
- **Not escalation:** under an arms race, future parasites would be more harmful than the host's
  own. They never are on median; the further in the future, the less harmful (Spearman +0.40).
- **Not cycling either:** old parasites do not regain their effect at long lags. Past
  parasites are less harmful by a roughly constant margin however old (Spearman +0.10).
- **Two phases.** Before about 1M everything changes fast (large d off the diagonal). From 1.2M
  on, parasites from 1.2-2M are interchangeable in late worlds (median d +0.02, > 0 in 99/160),
  while parasites from 1M or earlier stay less harmful there (+0.12, > 0 in 173/200).
- Runs differ: in runs 1, 3 and 5 past parasites are the less harmful; in 4, 6, 7 and 8 future
  ones.

**Reading:** parasites are adapted to their own time, and replicators move away from them for
good, but parasites do not get generally better: directional turnover (Red Queen without
recurrence), fast in the first million and nearly settled in the second. This fits the fixed-early
self-scan (section 21), the arms race of section 26 (parasites adapt to whatever replicator is
common) and the stalling of section 28, and it does not fit escalating complexity.

**Caveats:** d is relative to population, so small early worlds give large values; classes come
from one 20k window; only parasites are swapped (the replicator half was done for 1M/2M only).

## 32. QSEARCH: a "password" label test, and whether it brings escalation (Nov 21)

Details: `experiments/parasites/qsearch_summary.py`, `qsearch_summary.txt`, `timeshift_smol.sh`,
`ts_ctl.results.txt`, `ts_qsearch.results.txt`, `mutassay_qsearch.txt`. Runs on aida in
`~/smol/runs/qsearch/`, `ts_ctl/`, `ts_qsearch/`.

**Why:** section 31 found turnover, not escalation. One explanation: recognition is too narrow.
A `?` label test (2+ letters) compares its label only with the partner at the read pointer, where a
deletion-mutant parasite usually matches. Label tests do appear in the published runs (molecules
carrying one: ~0 at 100k, 15-36% by 2M in runs 1, 3, 4, 8), so evolution reaches for them.

**Switch:** `QSEARCH 1`: a `?` label test passes if its label matches anywhere in the strand under
the read pointer (best local match, same probability rule): a check that the partner carries a
tag, which tag and check can change together. With QSEARCH 0 a run is identical to before
(checked to 60k); with 1 it diverges once a label test first runs (9k in seed 201).

**Runs:** the section 23 control set-up and seeds 201-210 plus QSEARCH 1, 2M steps. Time-shift
matrices (as section 31, from smol's own checkpoints, 2 seeds, 20k steps) for these runs and for
the 10 controls.

| | ctl | QSEARCH |
|---|---|---|
| Survive to 2M | 5/10 | 7/10 |
| Molecules with a label test, median, 1.2-1.4M | 0.011-0.041 | 0.119-0.163 |
| Top species carries one | 0 of 5 throughout | 0-2 of 8 |
| Population at 2M (median) | 10,628 | 11,275 |
| Parasite share at 2M (median) | 0.27 | 0.26 |
| Own-mutant discrimination after 1.5M | ~0.25-0.32 | 0.33-0.46 (3 runs) |
| Time-shift: past parasites, median d (> 0) | +0.13 (86/105) | +0.08 (279/354) |
| Time-shift: future parasites, median d (> 0) | +0.05 (88/129) | +0.13 (255/354) |

- Escalation would make future parasites more harmful (d < 0). With QSEARCH they are less
  harmful, the more so the further ahead (+0.03 at +0.2M rising to +0.43 at +1.6M), and past
  parasites are less harmful by a near-constant margin: the same turnover as sections 31 and the
  smol controls. From about 1M the late block is nearly flat (d +0.0 to +0.15).
- Password checks spread somewhat more but stay in minor species; discrimination of own mutants
  rises a little late in the 3 runs assayed.

**Reading:** giving recognition a much larger space did not produce escalation within 2M steps.
The checks are taken up only weakly, probably because deletion-mutant parasites keep the tag and
each check costs time. So the ceiling is not (only) the narrowness of recognition.

**Caveats:** 10 runs, 2 seeds per matrix cell; the control matrix rests on the 4-5 surviving
control runs; the mutant assay covered 3 QSEARCH runs.

**Housekeeping:** the two matrices first filled aida's disk (323 GB: every swap run kept a copy of
the run's large species list). With David's OK the copies of species.tsv, copies.tsv and grid.snap
in the swap-run folders were deleted (logs, popdy, checkpoints and classes kept); both time-shift
scripts now delete them as each run finishes.

## 33. Gene and machine: does separating machinery from template bring escalation? (Nov 21)

Details: `experiments/parasites/classify_gm.py`, `gm_summary.py`, `gm_summary.txt`,
`ts_gm.results.txt`, `timeshift_extra.py`, `timeshift_extra_results.txt`. Runs on aida in
`~/smol/runs/gm/`, `ts_gm/`.

**Set-up:** the section 8 design (gene G, 100 letters, inert as code; machine M, 42 symbols, copies G
blind and translates it into a new M), with STRANDCOPY, TRANSLATE and TRANSLATE_STOP, the old
partner and turn rules and COMPAT_DELETION (as the controls), 10 seeds, 2M steps.

**Classes (`classify_gm.py`):** over a 20k window, R = makes exact copies (machines); species that
are copied but never copy are a working gene (G) if one of the three commonest machines, bound to
it, copies it and translates it into a machine that, bound to it again, copies it and makes itself;
otherwise a gene-parasite (P). A first version required genes to be all letters; evolved genes
carry the odd non-letter outside the coded part, so most working genes were classed P and the first
swap set was invalid (kept in `~/smol/runs/ts_gm_invalid/`, not deleted). The tables below are
from the corrected classification.

| Step | 200k | 600k | 1M | 1.4M | 2M |
|---|---|---|---|---|---|
| Runs alive | 10 | 10 | 10 | 10 | 10 |
| Population (median) | 9,061 | 9,336 | 9,312 | 8,720 | 9,166 |
| Machines / working genes / gene-parasites | .27/.31/.36 | .27/.44/.23 | .27/.43/.22 | .26/.40/.28 | .26/.28/.43 |
| Commonest machine length | 42 | 42 | 42 | 42 | 42 |

- All 10 runs survive (controls 5/10). The machine barely changes: 42 symbols and the same core
  in every run (e.g. run 301 at 2M: `WGXWEHHR=XEK^B?WGEW}$=?>%^F$LUBO>B^C$=?>P`).

**Time-shift matrix** (as section 32, swapping gene-parasites only; 1,000 cells, 2 seeds):

| | gene and machine | controls | QSEARCH |
|---|---|---|---|
| d by lag, all hosts: Spearman (t) | -0.15 (-4.6) | -0.12 (-1.8) | +0.10 (+2.6) |
| past / future parasites, median d | +0.030 / -0.025 | +0.128 / +0.052 | +0.083 / +0.129 |
| hosts and parasites both from 0.6M: Spearman | +0.03 | -0.18 | +0.12 |
| parasites from 0.2-0.4M vs other times, median d | +0.06 to +0.07 | +0.08 to +0.26 | +0.10 to +0.23 |
| parasites from 0.6M on vs other times, median d | -0.03 to +0.01 | +0.01 to +0.31 | +0.04 to +0.15 |

- Future gene-parasites are more harmful than a world's own (the escalation sign), but only
  because parasites from 0.2-0.4M are weaker: gene-parasites improved once, by about 0.6M.
- From 0.6M on, parasites from any time are as harmful as any other (d ~ 0, no lag trend): no
  local adaptation in time, no turnover, no further escalation.
- In the ordinary-replicator systems parasites from every other time are less harmful than a
  world's own (turnover).

**Reading:** separating machine from template removes the turnover rather than producing an arms
race. The host side (the machine) does not evolve, so gene-parasites improve once and then have
nothing to track; the system is stable (10/10 survive) and static. An arms race needs both sides to
keep changing; here the machine is frozen, probably because any change to it must come through a
translated gene mutation and the blind copier has no recognition step to improve.

**Caveats:** 10 runs, 2 seeds per cell; the 2M host row is odd (old parasites more harmful,
d -0.06 to -0.18), unexplained.

## 34. Why is the gene-and-machine machine frozen? (Nov 21)

Details: `experiments/parasites/gm_scan.py`, `gm_scan2.py`, `gm_scan3.py` (seed gene and machine,
`gm_scan*_seed.txt`), `gm_linkage.py`, `gm_linkage_results.txt`.

Section 33: the commonest machine stays at 42 symbols with the same core in all 10 runs. Three
candidate causes, tested:

**1. A fragile genetic code? No.** Every single-letter substitution in the coded part of the seed
gene (2,100), translated by the seed machine, then tested for a closed loop:
- synonymous 1.6% (the balanced random table keeps the symbol for 2.9% of single-letter codon
  changes, and the letter/opcode kind for 66.5%); a different working machine 38%; broken 60%.
  For comparison, 25-45% of an ordinary replicator's single mutants still copy themselves.
- But by position: the label head (`WGEWLHHRJXQQ`) tolerates 92-98% of substitutions; the working
  core (`^B?WGEW}$=?>%^F$LUBO>B^C$=?>%}`) tolerates 0-6% at most positions (exceptions: `}`,
  and a few symbols near the end).
- Single insertions and deletions shift the reading frame: 95% broken.
- In-frame duplications inside the core work in 19-60% of cases (2-30 letters), including
  duplications of the copy loop (e.g. `...}$=?$=?>%...`), the kind of change that built self-scan in
  ordinary replicators. So structural variation of the core is reachable.

**2. No genotype-phenotype linkage (a better machine helps all genes)? No.** Over 1-2M, each of
the 20 commonest machines copies the gene that encodes it 23-55 times (run medians) more often
than that gene's share of all gene copies would give: products are placed next to their makers,
so space links a machine to its own gene.

**3. Remaining explanation: the machine sits at an optimum.** Viable variants exist and are
linked to their genes, yet none spreads, so they are neutral (head edits, which do drift) or
worse (a longer machine translates more slowly). A blind copier has no recognition step that a
change could improve against gene-parasites (the label search is passed by any gene-parasite that
keeps the label). This is the "short and fast wins" pressure again, now on the machine.

**Reading:** a more robust codon table would not unfreeze the machine; the test of whether
machine change can be selected would be a competition, as in section 24: a duplicated-loop
machine against the seed machine, with and without gene-parasites.

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
