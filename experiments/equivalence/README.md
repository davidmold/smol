# Reaction-level equivalence with the C++ Stringmol (15 Nov)

Question: given the same two molecules, does smol do exactly what Stringmol does?

## Method

- **Molecules:** the 12 most abundant species in each of the 8 published runs (Hickinbotham,
  Stepney & Hogeweg 2021) at every 100k steps from 100k to 2M (`mkpairs.py`, run on aida). That
  gives 160 sets.
- **Pairs:** within each set, every ordered pair reacts, each molecule taking a turn as
  initiator. That is 23,040 reactions, about 25 million execution steps.
- **C++ side:** `spreact_*.cpp` does what the spatial simulation does:
  1. compute the binding probability and choose the active molecule (`AgentsAlign` +
     `ReactionSetupExecution` in dev; `get_sw` + `set_exec` in 0.2.3.4);
  2. place the pair on a grid;
  3. call the spatial step function (`ReactionExecuteOpcode_Spatial` / `spatial_exec_step`)
     until the reaction ends.
- **smol side:** `smol-react --compat --steptrace` does the same.
- **Output:** both print the same text, so the two outputs can be compared byte for byte:
  - the binding probability and which molecule is active;
  - before every step, the opcode and all four pointers (string and position), plus both
    string lengths;
  - every product;
  - how the reaction ended, and both final strings.
- **Random numbers:** both use MT19937 seeded with 4357, reseeded for every pair, so random
  draws (the soft `$` search, the `?` label test, mutation) match draw for draw. The C++ probe
  undoes the draws a cleave spends placing its product on the grid, since smol-react has no grid.
- **Mutation:** each comparison ran with mutation off and with MUTATE 0.01.

## Builds

```
# dev (uoy-research/stringmol upstream/dev 109be57), from a worktree at dev/
g++ -O2 -std=c++11 -w -Idev/src -o spreact spreact_dev.cpp \
    $(ls dev/src/*.cpp dev/src/*.c | grep -v "stringmol.cpp\|webapi_") -lm
# 0.2.3.4 (tag), from a worktree at v0234/
g++ -O2 -std=c++11 -w -fpermissive -Iv0234/src -o spreact spreact_v0234.cpp \
    $(ls v0234/src/*.cpp v0234/src/*.c | grep -v "stringmol.cpp\|webapi_\|mt19937-2.c\|tests.cpp\|smspatial.cpp") -lm
./cmpg.sh TAG ./spreact eq.conf "" PAIRFILE                 # mutation off
./cmpg.sh TAG ./spreact eqm.conf "--mutate 0.01" PAIRFILE   # MUTATE 0.01 in the conf
```

## Results

Per-set results are in `results_*.txt`. Every differing reaction is classified in
`differences.txt`.

| C++ version | mutation | identical reactions | differing reactions | cause |
|---|---|---|---|---|
| dev | off | 23,023 / 23,040 | 17 | C++ writes past its 2,000-symbol buffer |
| dev | 0.01 | 23,022 / 23,040 | 18 | same |
| 0.2.3.4 | off | 22,885 / 23,040 | 155 | uninitialised alignment score (fixed in #50) |
| 0.2.3.4 | 0.01 | 22,889 / 23,040 | 151 | same |

- **Against current dev, every reaction matches step for step, except at the length limit.**
  - What matches: the binding probability, the choice of active molecule, every pointer at
    every step, every product and every outcome. With 1% mutation the insertions,
    substitutions and deletions also match draw for draw.
  - The exceptions are reactions where a runaway copy fills the 2,000-symbol buffer. The traces
    still agree there. Only the final string differs: the C++ code writes one symbol past its
    buffer, losing the terminating zero, so the string it reports runs into whatever memory
    follows. smol stops at 2,000.
- **Against the paper's version (0.2.3.4), the only differences come from one bug.**
  - Every one of them is a `?` executed while the read pointer is at the end of a string.
  - In 0.2.3.4 that alignment reads an uninitialised score, so whether the jump fires depends
    on what the build left on the stack. This is the bug fixed in uoy-research/stringmol#50.
  - The build here jumps. smol and dev treat the empty alignment as no match. The paper's own
    binary could have gone either way.
- **One rule difference found and closed.** In both C++ versions, a deletion during copying
  also skips the next instruction (`OpcodeCopy` advances I on deletion, then again at the end).
  smol's simulator did not do this; only `smol-react --compat` did. It is now a switch,
  `COMPAT_DELETION` (default 0, so earlier smol runs are unchanged). `tools/pub_to_ckpt.py` turns
  it on along with `COMPAT_PARTNER` and `ENCOUNTER_TURN`.

## Not covered here

This compares reactions, not whole worlds. World-level rules still differ in implementation,
and have only been compared statistically so far:

- the order molecules act in each step;
- partner choice (`COMPAT_PARTNER`) and failed-encounter turns (`ENCOUNTER_TURN`);
- decay (a bound pair decays together in both);
- where products are placed.

The decisive test would replay a whole C++ run draw for draw.
