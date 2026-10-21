# Stringmol / smol: start here

Status as of **20 October 2026**. This folder is the one place to pick the project up from.
It summarises everything; the long design documents on claude.ai are linked at the bottom
and are not needed to follow what has been done.

## The question

Can a chemistry of self-copying strings keep producing new kinds of replicator, or does it
always settle into one ecology and stop? We use Stringmol (Hickinbotham, Stepney et al.,
University of York) because molecules bind by sequence and copy each other through explicit
opcodes, so parasites and arms races can arise without being designed in.

Working assumption (from the roadmap): **the ratchet is expensive, not just undiscovered.**
The useful experiments measure what makes a population escape a stuck state, rather than
waiting for complexity to appear.

## Where things stand, in one paragraph

Baseline Stringmol is not open-ended: it makes one replicator transition (the 64-symbol seed
shrinks to ~30 symbols), then settles into replicator-plus-parasite churn, with *less*
novelty over time than a neutral-drift null. Every chemistry change tried since (strand
displacement, two-sided reactions, a binding rule, copying onto a separate strand,
translation from gene to machine, environmental damage, cells) has left the same pressure in
place: **short and fast wins whenever nothing stops it.** A 10x bigger world is sturdier and
more diverse but invents new families at the same rate per unit area. The one genuinely new
interaction is in the copy-strand chemistry: **obligate partnerships** where a tiny borrower
(e.g. `F^O>`) runs its partner's copy code. They are still "short and fast" in disguise.
The current exploratory thread is side-on (VERTICAL) worlds with support, light and
chemistry-dependent decay, aiming at tree-like structure; so far fragile molecules win and
the only stratification seen was a lid of stable junk smothering the copiers.

## Files in this folder

| File | What it holds |
|---|---|
| [findings.md](findings.md) | Everything learned, experiment by experiment, with the numbers |
| [experiments.md](experiments.md) | Every run condition: where it is on aida, how it was configured, its outcome |
| [switches.md](switches.md) | Every config key of the simulator, its default and what it does |
| [tools.md](tools.md) | Machines, code locations, binaries, analysis scripts, dashboard, deploying, data handling |
| [next.md](next.md) | Open questions and candidate next steps |

## Machines and code at a glance

- **chasca** (WSL, David's PC): the Rust simulator `smol` in `~/projects/smol` (git, local
  only: **no remote, so no off-machine backup**). Edit and test here.
- **aida** (`ssh aida`, Ubuntu, 12 threads, 2x RTX 5060 Ti): runs everything. `~/smol` is a
  deployed copy of the Rust code plus all runs in `~/smol/runs`. `~/stringmol` is the C++ fork
  and the Python analysis scripts (families, novelty, drift null).
- **Dashboard / live viewer:** http://192.168.1.129:8095 (runs list) and
  `http://192.168.1.129:8095/live?run=<condition>/<run>` (live grid while a run is going).

## Running right now (on aida)

| Run | What | Progress (20 Oct, ~20:15) |
|---|---|---|
| `big/classic_400x300` | classic chemistry, 120,000 cells | 1.66M steps, 107k molecules |
| `big/strand_400x300` | copy-strand chemistry, uniform | 1.34M steps, 118k molecules |
| `big/gradient_400x300` | copy strand + reaction-rate gradient | 1.68M steps, 113k molecules |
| `vert/strand_400x12` | side-on world, light by height | 10.3M steps |
| `vert/slab_400x12x5` | side-on slab, 5 slices deep, shading | 0.53M steps |

The big runs are set for 100M steps (2-3 weeks); the vertical ones for 20M.

## The long-form documents (claude.ai)

These are living documents; this folder summarises them as of the date above.

- [Stringmol roadmap](https://claude.ai/artifact/UrRNvz3KaX2CGM7Jh71URX): purpose, tracks,
  dated results. Its "Where things stand" and checklists date from 16 Oct and are partly stale;
  the "What we've learned" section is current to 20 Oct.
- [Stringmol in Rust: specification](https://claude.ai/artifact/9r83UUx8AEroMoQtz1zF5u): the
  chemistry rule by rule, the old code's quirks and what was kept or fixed.
- [Granular copy: design note](https://claude.ai/artifact/3F4TKKy17b99AWYe2NfLBD):
  displacement, two-sided reactions, copy strand, with results.
- [Giving Stringmol a constructor](https://claude.ai/artifact/UQsfTrA3EH7NmMspE2NtdN): the
  von Neumann gene + machine design (TRANSLATE).
