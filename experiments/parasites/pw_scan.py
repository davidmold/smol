#!/usr/bin/env python3
"""pw_scan.py [SET...] -> the translation-start password of each run's commonest machine at every
classification window, and the transitions between consecutive windows.

The machine's translate instruction searches its gene for the complement (ROT13) of the label
after `^F$` (`LUBO` in the seed; the seed gene carries YHOB at the end of its 16-letter head).
A transition is a change of the commonest machine's label between two consecutive windows
(t, t + 200k). Prints, per set, the label sequence of each run and a list of transitions."""
import sys, os, re, glob, collections
R = os.path.expanduser("~/smol/runs")
comp = lambda s: "".join(chr((ord(c) - 65 + 13) % 26 + 65) for c in s)
def label(m):
    x = re.search(r"\^F\$([A-Z]+)", m or "")
    return x.group(1) if x else None
def window(w):
    seq = {}
    for l in open(f"{w}/species.tsv"):
        p = l.rstrip("\n").split("\t")
        if p[0].isdigit() and len(p) >= 5: seq.setdefault(int(p[0]), p[-1])
    by = collections.defaultdict(list)
    for l in open(f"{w}/classes.tsv"):
        a, b, c = l.split()
        if int(a) in seq: by[b].append((int(c), seq[int(a)]))
    for k in by: by[k].sort(reverse=True)
    return by
def runs(name):
    ts = f"{R}/ts_{name}"
    for rd in sorted(glob.glob(f"{ts}/s*"), key=lambda x: int(x.split("/s")[-1])):
        ws = sorted(glob.glob(f"{rd}/w*k"), key=lambda x: int(re.search(r"w(\d+)k", x).group(1)))
        out = []
        for w in ws:
            if not os.path.exists(f"{w}/classes.tsv"): continue
            t = int(re.search(r"w(\d+)k", w).group(1)) * 1000; by = window(w)
            M = by["R"][0][1] if by["R"] else None
            nR = sum(c for c, _ in by["R"]); top = by["R"][0][0] / nR if nR else 0
            out.append((t, M, label(M), top, by))
        yield os.path.basename(rd), out
if __name__ == "__main__":
    for name in sys.argv[1:] or ["cf_gmt", "cf_gmq", "cf_gmd"]:
        n = 0; print(f"== {name}")
        for r, out in runs(name):
            labs = [x[2] or "-" for x in out]
            print(f"  {r}: " + " ".join(labs))
            for a, b in zip(out, out[1:]):
                if a[2] and b[2] and a[2] != b[2]:
                    n += 1
                    print(f"     transition {a[0] // 1000}k -> {b[0] // 1000}k: {a[2]} -> {b[2]} (top machine share {a[3]:.2f} -> {b[3]:.2f})")
        print(f"  {n} transitions")
