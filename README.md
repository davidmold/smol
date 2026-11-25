# smol

[![DOI](https://zenodo.org/badge/1459239374.svg)](https://doi.org/10.5281/zenodo.24443699)

A fast reimplementation, in Rust, of the spatial **Stringmol** automata chemistry
(Hickinbotham, Stepney et al., University of York;
[uoy-research/stringmol](https://github.com/uoy-research/stringmol)). It was written from a
description of the model rather than translated from the C++ code, and is used for long and
large evolution experiments: a full 125 x 100 world runs at roughly 200 steps per second on one
core, and much larger worlds are practical.

It reads Stringmol-style configs (`GRIDX`, `GRIDY`, `DECAY`, `MUTATE`, `RANDSEED`, `AGENT`,
`GRIDPOS`, ...), plus switches for the variations studied here (`docs/switches.md`).

## What it has been used for

The main study so far tests the account of parasitism in Hickinbotham, Stepney & Hogeweg (2021):
**[`docs/parasitism.md`](docs/parasitism.md)** summarises what the evidence shows, with
`docs/findings.md` and `experiments/` holding every script and raw result behind it. In brief:
parasites are real and limiting, but in the published setup they settle into a turnover rather
than an escalating arms race. A gene-and-machine variant whose machine checks a second tag on its
gene first looked like a continuing arms race, but a 110-run confirmation with its tests fixed in
advance did not support it (`docs/findings.md` section 36). All the gene-and-machine variants
show changing host-parasite compatibility. Contests of evolved before-and-after pairs suggest that
turnover is mostly hitchhiking: the machine's translation-start label changes with little effect
of its own, carried by machine-gene pairs that are better overall. Lineage replays in the earlier worlds
support this: evolved backgrounds helped rare founders establish whichever label they carried.
Any part played by parasites is small and unresolved (section 37).

## How it relates to Stringmol

**The same chemistry, not the same program.**

- Matches: the opcode set and execution, binding by Smith-Waterman alignment, copying with
  point mutation, decay, a toroidal grid with Moore-neighbourhood reactions, products placed
  in an empty neighbouring cell. Our reaction probe reproduces the reaction times reported in
  Hickinbotham, Stepney & Hogeweg (2021, *R. Soc. Open Sci.* 8: 210441) for strings from their
  archived runs (184, 263, 225 steps against their 183, 262, 224), and over the first million
  steps runs of their configurations look statistically like theirs.
- Differs: a different random generator and update details, so it does not reproduce a C++
  run step for step and cannot continue one exactly. Two behaviours of the C++ code are
  switches, off by default: `COMPAT_PARTNER 1` (first free neighbour as partner) and
  `ENCOUNTER_TURN 1` (a failed binding attempt uses up the partner's turn). Alignment of an
  empty string scores 0 (as in uoy-research/stringmol#50).
- In their second million steps the published runs show parasitism falling and mutual copying
  rising; reruns of their configurations here show the same direction but weaker, even with
  the old-rule switches (few surviving runs, so not yet a measured difference; `docs/findings.md`
  section 12).
- It does not read the C++ code's `out1_*.conf` snapshots or write its log formats; its outputs
  are `popdy001.dat` (step, species, count), `species.tsv`, `copies.tsv` and gzipped checkpoints.

## Build and run

```
cargo build --release
mkdir run && cd run
../target/release/smol run ../designs/base.conf      # the baseline 125 x 100 world
../target/release/smol resume checkpoint_500000.txt.gz [--reseed S] [--nsteps N]
../target/release/smol-react pairs.txt 5000 --self  # reaction probe
```

## Documentation

`docs/README.md` is the place to start: the research question, where things stand, and links
to the findings (`docs/findings.md`), experiment index, switches and tools. The notes are a
working lab book, so they mention the machines and storage the runs used.

## Licence

GPL-3.0-or-later (see `LICENSE`). smol is a separate program, not derived from Stringmol's code.
(Stringmol's repository ships the GPL v2 licence text, while its source headers say GPL v3 or
later.)
