#!/usr/bin/env python3
"""pw_plot.py SCORES READOUT TRANSITIONS OUT.png -> per-transition view of the round-2 password
contests (no parasites): the whole later pair's advantage and each available label swap, marked by
whether the label change is silent in translation (filled) or moves/weakens the start (hollow) and
whether it changes the label's length (orange name). Panel A: the transitions with both swaps."""
import sys, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
s = json.load(open(sys.argv[1])); r = json.load(open(sys.argv[2])); t = json.load(open(sys.argv[3])); OUT = sys.argv[4]
rows = []
for x in r:
    k = x["k"]; g = lambda c: s.get(f"{k}:{c}", {}).get("p0")
    rows.append(dict(k=k, old=x["old"], new=x["new"], sil=x["old_gene"]["silent"] and x["new_gene"]["silent"],
                     main=g("main"), fwd=g("fwd"), recip=g("recip")))
BL, OR, AQ, TX, T2 = "#257de0", "#f47230", "#00a888", "#191919", "#666666"
def panel(ax, rs, title):
    rs = sorted(rs, key=lambda x: -x["main"])
    ax.axhline(0, color="#bbbbbb", lw=1, zorder=0)
    for i, x in enumerate(rs):
        ax.plot([i], [x["main"] - .5], "o", ms=8, color=BL, zorder=3)
        for key, col, dx, mk in (("fwd", OR, -0.2, "s"), ("recip", AQ, 0.2, "D")):
            if x[key] is not None:
                ax.plot([i + dx], [x[key] - .5], mk, ms=7, mfc=col if x["sil"] else "white", mec=col, mew=2, zorder=4)
        dl = len(x["new"]) - len(x["old"])
        ax.text(i, -0.62, f"t{x['k']:03d} {x['old']}→{x['new']}" + (f" (label {dl:+d})" if dl else ""),
                rotation=60, ha="right", va="top", fontsize=8, color="#c2410c" if dl else TX)
    ax.set_xlim(-0.7, len(rs) - 0.3); ax.set_ylim(-0.6, 0.6); ax.set_xticks([])
    ax.set_ylabel("gene share at 30k minus 0.5", color=T2, fontsize=9)
    ax.set_title(title, loc="left", fontsize=11, color=TX)
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#bbbbbb"); ax.tick_params(colors=T2, labelsize=8); ax.grid(axis="y", color="#eeeeee", lw=0.8)
both = [x for x in rows if x["fwd"] is not None and x["recip"] is not None]
one = [x for x in rows if (x["fwd"] is not None) != (x["recip"] is not None)]
fig, axs = plt.subplots(2, 1, figsize=(13, 11.5), gridspec_kw={"height_ratios": [1, 1.15]})
panel(axs[0], both, f"A. The {len(both)} transitions with both swaps (no parasites)")
panel(axs[1], one, f"B. The {len(one)} transitions with one swap (no parasites)")
h = [Line2D([], [], marker="o", ls="", color=BL, ms=8, label="whole later pair vs earlier pair"),
     Line2D([], [], marker="s", ls="", mfc=OR, mec=OR, ms=7, label="new label on earlier pair (forward swap)"),
     Line2D([], [], marker="D", ls="", mfc=AQ, mec=AQ, ms=7, label="new label on later pair (vs reciprocal swap)"),
     Line2D([], [], marker="s", ls="", mfc="white", mec="#888", mew=2, ms=7, label="hollow: label change moves or weakens the translation start"),
     Line2D([], [], ls="", label="orange name: label changes length")]
fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, fontsize=9)
plt.subplots_adjust(top=0.9, hspace=0.8, bottom=0.17, left=0.07, right=0.99)
plt.savefig(OUT, dpi=150)
