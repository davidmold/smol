# Parasitism in Stringmol: what the evidence shows

Nov 17, 2026 · David Mold

> Snapshot of the Claude Doc of the same name, exported 18 Nov 2026. The doc is the working version; this copy is re-exported when it changes. Details behind every number are in `findings.md` (sections 17–30) and `experiments/`.

## Summary

Parasitism is real and limiting, but the evidence supports only part of the story in Hickinbotham, Stepney & Hogeweg (2021, *R. Soc. Open Sci.* 8: 210441). Parasites hold the population down, and they are what make self-scan win. But self-scan itself is built by a mutational route through parasite-like intermediates, later gains come mostly from replicators stalling on established parasites rather than from better discrimination (which stays weak), and parasite lineages persist and adapt to new defences: an arms race, not steady improvement.

| Claim in the paper | Verdict | Key evidence |
| --- | --- | --- |
| Parasites limit the population | Supported, by intervention | Removing parasites raised the population in 8/8 runs, by a median of 3,087 molecules at +10k steps |
| Evolved resistance explains the larger population later on | Partly | Later replicators stall on established parasites (held without copying: 8% → 59% of their bound time), which tracks the swap results. But they still copy their own fresh mutants, so parasite share mostly holds steady from 1M to 2M. |
| Replicators evolve explicit, increasingly complex self/non-self discrimination | Partly | It clearly evolves in 2 of 8 runs; in the rest, own-mutant parasites are copied about 2/3 as readily as kin, for 2M steps |
| Replicators slow replication (self-scan) to blunt parasites | Supported, and graded | Self-scanners beat their non-scanning contemporaries when parasites are present in 4 of 6 competition runs. In the other 2, the parasites had already adapted to the scanner. Its hold on the top spot rises with how far parasites can multiply (dose-response). |
| New parasites arise from mutated replicators, not from parasite lineages | Not supported | About half of parasite species, and 56–72% of parasitism, descend from other parasites; lineages last 150k–440k steps |
| Evolution evolves: higher mutation rates and emergent mutation operators | Real, but fixed early | Mutation exposure doubles with self-scan; about 80–86% of new strings are rearrangements; neither keeps rising after 1M |
| Parasitism drives complex replicators and ecosystems | Partly | Self-scan, shortening and a second `?` all arise without parasite lineages (PARASITEFREE 2, 10 of 10 runs). Parasite pressure decides whether self-scan dominates, in proportion to how far parasites can multiply. Parasites don't get better at exploiting in general, but they do adapt to a new defence. |

## Setup

All experiments ran in smol, a Rust reimplementation of spatial Stringmol, using the paper's own worlds and rules. It's fast enough to make these experiments practical: about 1,000–2,000 steps per second on a full 125×100 grid.

- **Equivalence with the C++ code.** 23,040 reactions between the paper's own molecules were traced step by step with matching random draws. Against current `dev`, every reaction is identical except where a string hits the 2,000-symbol limit (a C++ bug, fixed in uoy-research/stringmol#54). Against v0.2.3.4, the only differences come from the uninitialised alignment score fixed in #50.
- **Rule switches.** Runs use the paper's old partner and turn rules (`COMPAT_PARTNER`, `ENCOUNTER_TURN`), plus `COMPAT_DELETION`, the C++ copy-deletion quirk the equivalence test found.
- **Data.** The 8 published runs (snapshots every 20k steps to 2M) were converted to smol worlds with `pub_to_ckpt.py`. New runs started from the paper's seed replicator in its "box" configuration.
- **Classes.** A species is classed by what it does over a window of steps:
  - *replicator*: made an exact copy of something;
  - *parasite*: was copied exactly but never made a copy;
  - *other*: neither.
- **World level.** World-level rules (update order, partner choice, product placement) haven't been matched draw for draw, but a statistical check against the paper's runs found no significant difference (section 7, item 4).

Details, scripts and raw results are in the smol repo: `docs/findings.md` sections 17–30 and `experiments/parasites/`.

## 1. Parasites hold the population down

Removing parasites from the paper's worlds at 1M steps released the population in all 8 runs. Removing the same number of random molecules didn't, and nor did removing the non-copying, never-copied molecules.

- **Method:** parasites were 26–62% of molecules at 1M. Each removal was compared with removing the same number of random molecules, as a control for freed space. 3 seeds, 100k steps, mutation on.
- **Result:** population gain from removing parasites over the matched random removal (median of 8 runs; runs with a gain in brackets):

| Steps after removal | 2k | 10k | 20k | 50k | 100k |
| --- | --- | --- | --- | --- | --- |
| Extra molecules | +1,631 (8/8) | +3,087 (8/8) | +1,325 (7/8) | −6 (4/8) | +323 (6/8) |

- **Removing other non-copiers:** no gain (−430 at 10k, positive in 2/8 runs).
- **The effect is temporary:** it fades within about 50k steps as parasites return. In run 3 the population went from about 3,900 to 11,800 within 10k steps.

## 2. Replicators improve; parasites are matched to their time

Later replicators lift an earlier world for good, but parasites do worse whenever they're moved to another time.

- **Method:** one class (replicators or parasites) of each run's 1M or 2M world was replaced with molecules drawn from the same run at the other time. Each control replaced that class with fresh draws from its own time, so the disruption was the same. 5 seeds, 100k steps.
- **Result** (medians of 8 runs, as differences from the control):

| Swap | Population +10k | Population +100k | Parasite share +10k |
| --- | --- | --- | --- |
| 2M replicators into 1M worlds | +2,909 (8/8 up) | +1,745 (8/8 up) | −0.20 |
| 1M replicators into 2M worlds | −2,564 (0/8 up) | −713 (0/8 up) | −0.01 |
| 2M parasites into 1M worlds | +1,294 (5/8 up) | −117 (3/8 up) | −0.32 |
| 1M parasites into 2M worlds | +1,019 (8/8 up) | +491 (6/8 up) | −0.31 |

- **Replicators:** later ones are better in a lasting way.
- **Parasites:** they lose about a third of their share against replicators of another time, earlier or later, and recover as new parasites appear. They don't get better over time; they track the replicators of their own era.

## 3. Resistance by stalling, not by discrimination

Later replicators cope with the parasites already around by binding them and letting go without copying them. They still copy their own fresh mutants, so new parasites keep coming and parasite share stays up. (This section first concluded "mostly tolerance"; the stalling result below revises that.)

- **Two ways to cope:**
  - *Resistance* lowers how often parasites get copied, so their share falls.
  - *Tolerance* leaves parasite share unchanged but lowers the cost per encounter.
- **Parasite share mostly holds steady from 1M to 2M.** The clear exception is run 4, where it falls from 26% to 7%.
- **Swapping back costs population, not share:** putting 1M replicators back into 2M worlds lowered the population, with no change in parasite share.
- **Discrimination against their own mutants is weak.** For the top replicator every 100k steps, all its single-symbol mutants were made and reacted with it:
  - would-be parasites (mutants that can't copy themselves) bind exactly as often as working mutants (0.62–0.65);
  - they're copied in 65–73% of encounters, against 93–96% for working kin;
  - this 0.24–0.31 discrimination is flat from 200k to 2M.
- **Two runs do improve:** run 2 (to 0.5–0.7 from 1.1M) and run 4 (to 0.55–0.65 from 1.3M). Run 4 is also where parasites collapsed.
- **More checks aren't better checks:** extra `?` tests appear in runs 7 and 8 without better discrimination.
- **Other candidat****e mechanisms, 1M to 2M:**
  - *ruled out:* several copies per binding (still 1), and two-way copying (none);
  - *not supported:* tighter spatial clustering (1.34 to 1.24);
  - *what does change:* self-copy time falls from 178 to 153 steps (in 5 of 6 runs), while parasites keep most of their speed advantage (parasite/kin time ratio 0.73 to 0.77).
- **Stalling explains most of the swap result.** Across the 8 published runs, replicators increasingly bind established parasites, hold them for 50+ steps and release them without copying. That rises from 8% of parasite-bound time (up to 600k) to 59% (1.4M–2M), while copying falls from 65% to 30%. From 1M to 2M it rose in all 8 runs. Against the 2M parasites, the 1M top replicator copies where the 2M one holds in runs 1, 3 and 5, and the 2M replicator doesn't bind them at all in run 4. Those are four of the five biggest population drops when 1M replicators were swapped back in. The two runs with no difference (2 and 8) had the two smallest drops; run 7 is the exception. Details in section 7, item 5.

## 4. Does parasitism drive the complexity?

Yes, but not the way the paper describes. The route to self-scan runs through parasite-like intermediates whether or not parasites can multiply; multiplying parasites supply the selection that makes it win, in proportion to how far they can multiply. Each test below had a catch, now resolved.

**Suppression.** From the paper's 200k worlds, established parasites were removed every few thousand steps; a matched control had random molecules removed instead.

- 8 of 10 suppressed runs survived to 2M, against 4 of 10 controls.
- Established parasite share was 0.08, against 0.30 in controls.
- Catch: the removal rule spared parasites younger than 5k steps, so fresh parasites stayed. This arm can't say what happens with no parasites at all.

**PARASITEFREE.** smol places a new molecule only if it can copy itself exactly when bound to a copy of itself. From 200k:

- all 10 runs survived, at about 12,200 molecules; 5 of 10 controls survived;
- the top replicator stayed longer (60/42/38 symbols at 200k/1M/2M, against 31/28/28);
- it kept fewer `?` checks (1, against 2);
- it got much faster at copying itself (215 steps falling to 105–135, against 170–230);
- self-scan was absent in 28 of 29 snapshots, but present in every surviving control.

So when every non-self-copier is blocked, replicators get faster but stay long and never evolve self-scan. At first sight that supports the paper's claim, but the next test shows the block was cutting off the route to self-scan.

- Catch, now tested: the block also removed the stepping stones to self-scan. In every control run where self-scan appeared, it came through 1–8 non-self-copying intermediates. With those allowed but unable to multiply (PARASITEFREE 2), self-scan appears in all 10 runs, but dominates in only about half of them (see section 7, item 2).

**Competition.** A self-scanner (D) and a non-scanner (A) from the same run were mixed 50:50 in a fresh world, with and without parasites, mutation on and off. Share of D at 30k steps:

| Run | No mutation, no parasites | No mutation, parasites | Mutation, no parasites | Mutation, parasites |
| --- | --- | --- | --- | --- |
| 1 | 0.72 | 1.00 | 0.43 | 1.00 |
| 2 | 0.42 | 1.00 | 0.87 | 1.00 |
| 3 | 0.25 | 1.00 | 0.76 | 1.00 |
| 7 | 0.74 | 1.00 | 0.79 | 0.95 |
| 4 | 0.82 | 0.39 | 0.68 | 0.54 |
| 5 | 0.99 | 0.00 | 0.99 | 0.61 |

- **Runs 1, 2, 3 and 7:** self-scan wins outright when parasites are present, and is unremarkable without them. This is the cleanest support for the paper's claim that we have.
- **Why:** A copies parasites in 60–91% of encounters, D in 0–12%. D stays bound to a parasite for 5–26 steps, against 97–163 for A.
- **Runs 4 and 5 go the other way****, and that's now explained****:** their parasite pools came from snapshots where D already dominated, so they held parasites adapted to D. Parasites evolved on D in a fresh world beat D every time (35 of 35 pools). See section 7, item 1.

**Dose-response.** With the route to self-scan left open and only parasite multiplication varied, self-scan's hold on the top spot falls from 1.00 to about 0.5 as multiplication is cut, while replicators get faster and shorter and survive more often (section 7, item 3).

## 5. Where parasites come from: like defective interfering particles, but with lineages

The paper says new parasites "continually arise from mutated replicators, rather than from evolving parasite lineages". The first half holds; the second doesn't.

**Genealogy.** Each parasite species was traced back through its recorded parents to the nearest non-parasite.

- 42–58% of parasite species are one step from a replicator (depth 1), so new parasites do keep arising from replicators.
- Weighted by how often they actually got copied, depth 1 accounts for only 28–44%. **56–72% of parasitism comes from parasites descended from parasites.**
- Parasite clades reach 200–520 species, and lineages last 150k–440k steps, against 23k–85k for replicator species.
- Deep parasites aren't better exploiters than fresh ones (copied at about the same rate), so lineages persist rather than improve in general. They do adapt to whichever replicator is common, though: parasites evolved on a self-scanner beat it in 35 of 35 tests (section 7, item 1).

**Regeneration.** Three parasite-free worlds (from the PARASITEFREE runs) were restarted with the block switched off, 5 seeds each.

- With mutation on, parasites reappear within 1–2k steps.
- Population falls from 12,235 to 9,817 after 20k steps and 8,464 after 100k, then partly recovers to 10,362 at 200k. 2 of 15 runs went extinct.
- With mutation off, or with the block left on, population stays at about 12,200.
- By 200k, new parasite lineages already make up 41–66% of parasite species.

**Reading.** Like defective interfering particles in viruses, parasites are cheap deletion-type mutants that a working replicator is constantly producing and copying. Unlike them, they then go on to form long-lived lineages of their own.

## 6. "Evolution itself evolves"

There's only modest support for this claim in our runs.

- **Mutation exposure barely moves:** the average number of copying passes a sequence gets exposed to goes from 2.01 to 2.03 between 1M and 2M. That's because self-scan is already everywhere by 1M.
- **Self-scan is early:** it appears between 80k and 200k steps in most runs (run 5 at 580k), so most of the change happens well before the period the paper focuses on.
- **Rearrangements creep up:** among inexact copies, the share that are rearrangements (rather than point changes) rises from 0.80 to 0.86.

So the way variation is produced does shift, but mostly early, and only slightly after that. And because self-scan dominates only where parasites can multiply (section 7, item 3), the doubled mutation exposure is itself partly a product of parasite pressure.

## 7. Follow-up experiments

1. **Runs 4 and 5 in the competition (closed).** Self-scan wins against parasites adapted to other replicators, and loses to parasites adapted to itself.
   - **Causal test.** Parasites were evolved in fresh worlds holding only D or only A (150k steps, mutation on, 3 seeds, all six runs), then used in the A-vs-D competition. D-evolved pools beat D in 35 of 35 cases (median D share 0.00–0.01). A-evolved pools mostly favoured D (median 0.95–0.99), beating it in 7 of 28. As a side result, 4 of 18 A-only worlds died out from their own parasites, against none of 18 D-only worlds.
   - **Snapshot pools agree.** Pools from snapshots where D's family had taken over beat D in 11 of 13 cases. The original pools for runs 4 and 5 came from snapshots where D already outnumbered A by 5.7× and 4.2×.
   - **Run 5:** one fast D specialist, `$$YHOB`. D copies it 5.7 times for each copy of itself; A binds it but lets go in 8 steps.
   - **Run 4:** D's specialist `>$$OZHOB`, plus *decoys*: long relatives of A that A keeps copying but that can't copy themselves. They bind A's parasites for 160–360 steps without copying them, and release D's parasite within 10–45 steps (one also copies it). Removing the decoys flips run 4 completely: D's share goes from 0.39 to 1.00 in the 100k pool, and from 0.00 to 1.00 (all 10 seeds) in the 60k pool.
   - **Prediction rule.** Across all 93 pools, "the host with the fastest-growing parasite loses" gets the winner right in 85. The misses involve third-party molecules like these decoys, which pairwise measures can't see.

   Scripts and results: `ratio.py`, `evolve.py`, `predict.py`, `knockout.py`, `dropone.py` and their `*_results.txt`.

2. **Remove the stepping-stone confound (closed).** The route to self-scan runs through parasites, but parasite lineages aren't needed to make it.
   - **Design note.** The test first proposed here (place non-copiers but never copy them) couldn't work. A molecule that can't copy itself only has descendants if it is copied, so blocking that blocks the route too.
   - **Genealogy.** In the 6 control runs where self-scan appeared, the first self-scanner's ancestry runs through 1–8 non-self-copying intermediates. Typically a long mutant with a duplicated copy loop was copied by replicators, and the first self-scanner is a truncated copy of it.
   - **PARASITEFREE 2.** Every product is placed, but an exact copy of a molecule that can't copy itself never is. So such mutants appear but can't multiply. Same 10 seeds as before, 2M steps:
     - self-scan appears in all 10 runs (49k–288k), against 6 of 10 controls and none in PARASITEFREE 1 until 2M;
     - it gets there by the same route, with no intermediate ever copied exactly;
     - the top replicator shortens (17 symbols at 2M, against 28 in controls and 38 in PARASITEFREE 1), gains a second `?`, and copies itself fast (about 92 steps);
     - population stays near the parasite-free level (about 12,200; 1 run died at 1.97M);
     - **but** self-scan is the top replicator in only 5/10, 3/10 and 5/9 runs at 300k, 1M and 2M, against every surviving control.
   - **Reading.** PARASITEFREE 1's result came mostly from the blocked route. The structural changes (shortening, the second `?`, self-scan itself) come from the mutational route through non-copying intermediates, not from parasite pressure. What multiplying parasites add is selection: they make self-scan win and stay. This matches the competition, where D is roughly neutral without parasites and wins outright with them.

   Scripts and results: `scanpath.py`, `scanpath_ctl.txt`, `scanpath_pf2.txt`, `pfree_complexity_pf2.txt`, `pfree_expo_pf2.txt`.

3. **Dose-response (done: graded, levelling off at half).** PARASITEFREE 2 with a blocking probability q: every mutant is placed, so the route to self-scan stays open, and each exact copy of a non-self-copier is blocked with probability q. That varies only how far parasites can multiply. q = 0 is the controls and q = 1 is PARASITEFREE 2; same 10 seeds throughout, 2M steps.

   | Blocking q | 0 | 0.25 | 0.5 | 0.75 | 1 |
   | --- | --- | --- | --- | --- | --- |
   | Survive to 2M | 5/10 | 6/10 | 7/10 | 8/10 | 9/10 |
   | Runs where self-scan appears | 6/10 | 8/10 | 8/10 | 10/10 | 10/10 |
   | Top replicator self-scans (median share of each run's snapshots) | 1.00 | 0.80 | 0.50 | 0.40 | 0.63 |
   | Self-copy steps at 2M | 165 | 154 | 118 | 104 | 92 |
   | Length at 2M | 24 | 23 | 20 | 20 | 17 |

   - Self-scan's hold falls steadily as parasite amplification is cut, mostly between q = 0 and 0.5, then levels off at about half the time. Rank correlation over 41 runs: −0.39 (t = −2.6).
   - At the same time replicators get faster and shorter, and survival rises, all in steps with q.
   - **Reading.** Multiplying parasites select for self-scan in proportion to how far they can multiply, and self-scan costs copying speed. Without that pressure it is roughly neutral: it appears, and holds the top spot about half the time.
   - **Caveats:** at q = 0, only 4 runs live long enough to score, and survivors may be biased towards self-scan. Snapshots within a run aren't independent, so the test uses one value per run.

   Scripts and results: `dose.sh`, `dose_summary.py`, `dose_results.txt`; switch `PFPROB`.

4. **Does smol match the C++ code at world level? (statistical check done; a draw-for-draw replay isn't worth it.)** Reactions are identical. A whole-world replay would need smol to copy the C++ update order, list handling and bugs, and chaos would end agreement at the first mismatch, so we did a statistical check instead. The paper's three configs were rerun in smol with every compat switch on (8 seeds each, 2M steps) and compared with the paper's 8 surviving runs and our October rerun:
   - survival 7/24 against the paper's 8/20 (and 5/15 in October);
   - no significant difference in population, self-replication, mutual copying, parasitic share, reaction time or top length at 1.0–1.2M or 1.9–2.0M (every p ≥ 0.06 across 12 tests, survivors only);
   - **but** both smol reruns show a weaker second-million shift: parasitic share falls by a median of 0.04 (deletion fix on) and 0.08 (October) within each run, against 0.19 in the paper. Mutual copying at 2M is 0.28 against 0.53 (p = 0.12).

   So smol is indistinguishable from the paper at this sample size. The deletion rule doesn't explain the weaker shift, which may be noise, but a modest difference can't be ruled out with 7–8 survivors per side. Settling it would take about 100 smol runs plus fresh C++ runs, since the paper's 8 are all there is. None of our parasitism results depends on it: they use the paper's own worlds or compare smol conditions with each other. Details: `experiments/published/cmp3.py`, `cmp3_results.txt`.

5. **How tolerance works (closed: it is mostly stalling, and decoys are rare).** Survey of the 8 published runs at 100k, then every 200k to 2M (86 snapshots), with the 30 commonest species in each reacted pairwise. Replicator: makes an exact copy of something. Parasite: is copied exactly but never copies.
   - **Decoys are rare.** Non-replicators that hold parasites make up 0% of molecules in the median snapshot and 2% or more in only 10 of 86. They take a median 1–4% of parasite-bound time, and that share doesn't track population (rank correlation −0.05). The run 4 decoys were an unusual case.
   - **Replicators stall instead.** Replicators binding parasites, holding them 50+ steps and releasing them with no product: 8% of parasite-bound time up to 600k, 18% at 0.8–1.2M, 59% at 1.4–2M, while copying falls from 65% to 30%. Holding rose and copying fell from 1M to 2M in all 8 runs, and holding tracks population (rank correlation +0.47).
   - **It explains most of the swap result** (section 3): the runs where the 2M replicator holds or ignores parasites that the 1M one copies (1, 3, 4, 5) are four of the five biggest population drops.
   - **Reading.** Later replicators resist the established parasites by stalling on them, not by discriminating their own mutants. Fresh mutants are still copied 65–73% of the time, so parasites keep arising, and their lineages adapt (item 1).

   Caveats: one reaction trial per pair; only the 30 commonest species; the swap comparison uses only the top replicator. Not yet tested: a direct intervention on stalling. Scripts and results: `decoy.py`, `decoy_survey.txt`, `swapprobe.py`, `swapprobe_results.txt`.

**Still open.**

- A direct test of stalling: change how long replicators hold parasites and see whether population follows.
- Whether smol differs slightly from the C++ code over the second million (about 100 runs per side).
- The exact shape of the dose-response curve (more seeds per level).
