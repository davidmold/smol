# Config keys of `smol`

Config files use the old Stringmol format: `KEY value` lines, `%` or `#` comments,
`AGENT sequence count Q` followed by `GRIDPOS x y [z]` lines. Every experiment is a switch
that is **off by default**; with all switches off the output is byte-identical to the
baseline (checked after every change by running old and new binaries on the same configs and
comparing `popdy001.dat`, `species.tsv` and `copies.tsv`). The parsed config is copied into
every checkpoint, so a resumed run keeps its switches.

## World and run

| Key | Default | Meaning |
|---|---|---|
| `GRIDX`, `GRIDY` | 125, 100 | Grid size (a torus unless `VERTICAL`) |
| `DEPTH` | 1 | n > 1: n slices of GRIDX x GRIDY side by side (thin 3D slab, 3x3x3 neighbourhood, slices don't wrap); live view stacks the slices |
| `DECAY` | 0.0005 | Chance per step that a molecule (or a bound pair, together) is destroyed |
| `MUTATE` | required | Per copied symbol: indel with this chance, substitution with the same again |
| `MAXLEN` | 2000 | Longest allowed string; a write past it ends the reaction |
| `RANDSEED` | 1 | Seed for the run's random-number generator (Xoshiro256++) |
| `NSTEPS` | 1,000,000 | Run length; the run also stops on extinction |
| `POPDYEVERY` | 100 | Interval for `popdy001.dat` (step, species, count) |
| `SNAPEVERY` | 500 | Interval for frames in `grid.snap` (viewer) |
| `COPYEVERY` | 1000 | Interval for aggregated rows in `copies.tsv` |
| `CKPTEVERY`, `KEEPCKPT` | 50,000, 3 | Checkpoint interval and how many to keep |
| `LIVEEVERY` | 0 | Write `live.grid` / `live.json` for the live viewer every N steps |
| `CKPTKEEPEVERY` | 0 | Checkpoints at multiples of N are never deleted (for replays) |
| `OFFSPRINGLOG` | 0 | `copies.tsv` gains a 7th column, placed products (column 5 stays all attempts); `lives.tsv` per species and window: step, species, died, mean lifetime, died having made >=1 placed product, placed products made, died having had >=1 placed exact copy, placed exact copies. An exact copy is credited to the template (the partner) if it matches, else to the maker. Doesn't change the run |
| `SPECIESFORGET` | 0 | At each checkpoint forget sequences that are absent from the world and older than the previous checkpoint, to bound memory in long runs. A forgotten sequence that reappears gets a new species id (so `species.tsv` can list a sequence twice). Doesn't change the dynamics. Turn on mid-run with `tools/ckpt_setkey.py` |
| `TRACESEQ`, `TRACEAT`, `TRACEMAX` | off, all cells, 200,000 | Lineage trace: at start or resume, tag the molecules carrying the sequence TRACESEQ (only those in the cells listed in TRACEAT `x,y;x,y`, if given); placed products inherit the tag if they are an exact copy of a tagged parent, or a mutant closer (edit distance) to the tagged parent. `trace.tsv` gets one line per event for tagged molecules, up to TRACEMAX: `tag` (step, slot, born, species, x, y), `bind` (... role A/P, partner slot, born, species, tagged), `end` (... steps bound), `prod` (maker and partner slot/born/species/tagged, product species, exact copy of M(aker)/P(artner)/-, tagged, new slot), `lost` (a product with no room), `die` (... bound at death, age, products made, exact copies). Doesn't change the run |
| `DISTURB_EVERY`, `DISTURB_SIZE`, `DISTURB_MODE`, `DISTURB_UNTIL`, `DISTURB_SEED` | 0 (off), 7, 1, 0, 1 | Disturbance: every DISTURB_EVERY steps (until step DISTURB_UNTIL; 0 = no end) remove every molecule in DISTURB_SIZE^2 cells: one random square patch (MODE 1, clustered) or that many random cells across the world (MODE 2, dispersed). A removed molecule's partner is freed. Cells come from DISTURB_SEED and the step (resume-safe). The log gets `DISTURB to step N: removed K`. Off, outputs are byte-identical to before |

## Baseline rule switches (the old code's behaviour, for comparison)

| Key | Default | Meaning |
|---|---|---|
| `COMPAT_PARTNER` | 0 | 1 = old partner choice: first free neighbour in scan order |
| `ENCOUNTER_TURN` | 0 | 1 = a failed binding attempt also uses up the partner's turn |
| `PARASITEFREE` | 0 | 1 = a product is placed only if it can copy itself: bound to a copy of itself, mutation off, it makes an exact copy within 5000 steps in any of 3 trials (cached per sequence). Other products are discarded and counted as `blocked` in the log. Removes parasites at source; also removes molecules that only copy partners. 2 = only exact copies of such molecules are blocked: mutants that cannot copy themselves are placed (so mutational routes through them stay open) but never multiply. With `PFPROB` q < 1 (mode 2 only), each such exact copy is blocked with probability q (counted as `blocked`; the rest as `pfpassed`), for graded parasite amplification |
| `COMPAT_DELETION` | 0 | 1 = a deletion during copying also skips the next instruction, as the C++ code does (checked draw-for-draw against it) |
| `QSEARCH` | 0 | 1 = a `?` label test (2+ letters) passes if its label matches anywhere in the strand under the read pointer (best local match, same probability rule), not only at the read pointer: a "password" check on the partner. `smol-react --qsearch` |
| `DECAY_INDEPENDENT` | 0 | 1 = partners in a reaction decay separately; a survivor is freed |

## Binding

| Key | Default | Meaning |
|---|---|---|
| `BINDLEN` | 0 | Alignments shorter than n bind in proportion to their length |
| `BIND_T` | 0 | > 0: absolute-score rule P = 1 - exp(-score/T) |

## Reaction mechanics

| Key | Default | Meaning |
|---|---|---|
| `DISPLACE` | 0 | A free molecule may break into a reaction: P = DISPLACE_RATE x P(new) x (1 - P0) x max(0, 1 - copied/DISPLACE_K) |
| `DISPLACE_K`, `DISPLACE_RATE` | 10, 1 | Hardening and overall scale for displacement |
| `DISPLACE_RESTORE` | 0 | On displacement, undo what the abandoned reaction wrote |
| `TWOSIDED` | 0 | The passive also runs its own code, one instruction after the active each step |
| `STRANDCOPY` | 0 | `=` builds the copy on a separate nascent strand; the template must be the partner (own-side `=` does nothing); `%` releases the strand (or cuts at F if empty); an unfinished strand is released when the reaction ends |
| `TRANSLATE` | 0 | With STRANDCOPY: when `^C` has toggled W to the partner side, `=` reads two letters and writes their codon-table symbol (opcodes pass through) |
| `TRANSLATE_SEED` | 1 | Seed of the balanced 676 -> 33 codon table (`smol-codon table`) |
| `TRANSLATE_STOP` | 0 | In translate mode a non-letter pair is a stop codon: the reaction ends |

## Environment

| Key | Default | Meaning |
|---|---|---|
| `DAMAGE` | 0 | Chance per symbol per step that a symbol is replaced by a random one, any molecule |
| `DAMAGE_SHIELD_BOUND` | 0 | Molecules in a reaction are not damaged |
| `RATE_MIN`, `RATE_MAX` | 1, 1 | Reaction-rate field: a molecule's turn happens with probability rate. Across x (RATE_MIN at x=0, RATE_MAX in the middle column, smooth on the torus), or by height when `VERTICAL`. Decay is unaffected |
| `VERTICAL` | 0 | Side-on world: bottom row is the ground, no wrap in y; copies only go on the ground or on an occupied cell; a removed molecule drops its column by one |
| `SHADE` | 0 | With VERTICAL: rate x s for every molecule stacked above |
| `DECAY_CHEM` | 0 | S > 0: each symbol gets weight S^u (u uniform in -1..1, normalised to mean 1); a molecule decays at DECAY x the mean weight of its symbols; table written to `decay_weights.tsv` |
| `DECAY_CHEM_SEED` | 1 | Seed of that table |

## Specified but not built

From the specification document: `COMPCOPY` (copy writes the complement), `LENDECAY`
(length-dependent decay), `DIFFUSE` (molecules move), `MAXRXN` (reaction step limit),
`SOUP` (start from random strings).
