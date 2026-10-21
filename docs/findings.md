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

## 1. The baseline is not open-ended (Oct 16-17)

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

## The pattern across everything

Short and fast wins whenever nothing stops it. Changes to binding, interruption, reaction
symmetry, copying mechanics, translation, damage and space all left that pressure in place.
The only design that held its length was the designed cross pair; the only new interaction
that evolved by itself (partnerships) is itself a way of getting shorter.
