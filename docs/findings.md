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

## 1. The baseline consolidates over 1M steps (Oct 16-17)

- **One step forward, then churn.** In every run the 64-symbol seed replicase is replaced by
  one shorter family of 27-32 symbols (between 29k and 259k steps), which holds to the end.
  After that almost every new family is a parasite (1-16 symbols, copied by others, copying
  nobody). A second replicator family is rare and late (731k, 884k, 971k in 3 of 6 runs,
  never above 11%).
- **Selection suppresses novelty.** Sequence-family novelty per 100k-step window, against a
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

  Partnerships are reachable: they form spontaneously in most replays, more reliably closer to
  the original establishment time. When supplied, a pair establishes almost always, even early
  and even when foreign. So in these communities the limit on *when* partnerships appear is
  accessibility (the right variants arising), not invasion; and a community a few tens of
  thousands of steps older makes spontaneous establishment more likely. Caveats: 3 seeds;
  "established" counts any mutual pair, not only the introduced one (the introduced sequences
  were still present at the end in 46 of 54 intro runs).
- **Against the three possibilities** (hard to reach, can't invade, can't persist): for
  partnerships themselves, they are reachable and invade evolving communities easily; against a
  resident autonomous copier, invasion depends on the specific chemistry of the encounter; and
  under mutation the original pairs persist only as a minority. The next transition beyond
  pairs has not been observed, so these experiments don't yet say what blocks it.

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

Within 1M steps, short and fast wins whenever nothing stops it. Changes to binding, interruption, reaction
symmetry, copying mechanics, translation, damage and space all left that pressure in place.
The only design that held its length was the designed cross pair; the only new interaction
that evolved by itself (partnerships) is itself a way of getting shorter.
