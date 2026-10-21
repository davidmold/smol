# Published comparison: Hickinbotham, Stepney & Hogeweg (2021)

Data: their archive `biopaperdata.tgz` (CC BY; https://doi.org/10.15124/305dfdb6-9483-4c5b-8a01-c030570b9c31),
on aida in `~/published/` (extracted to `~/published/bio/biopaperdata/`). It holds the 8 runs that
survived to 2M steps: full state every 20k steps (`raw/runN/out1_STEP.conf`), their configs
(`script/runN/initrunN.conf`) and R analysis outputs.

Their configs: runs 1-4 "box" 125x100 with 121 seeds in an 11x11 block (the same as our
baseline `pr0` config), runs 5-6 100x125 with 51 seeds in random pairs, runs 7-8 250x50 with a
121 block. DECAY 0.0005, MUTATE 0.0002. Stringmol 0.2.3.4 (old partner and turn rules).

Converted with `tools/pub_convert.py` to `~/published/conv/runN` (popdy001.dat + species.tsv;
free molecules from AGENT lines, bound ones from the two `###` lines of each REACTION block).

- `cmp_pub.txt`, `cmp_pr0.txt`: `tools/length_traj.py` (N, median/p90 length, top species).
- `assay_pub.txt`, `assay_pr0.txt`: `tools/assay.py` per run every 100k (top 12 species reacted
  pairwise in `smol-react`, classic chemistry, mutation off). `pooled_*`: medians across runs.
