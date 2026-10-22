# -*- coding: utf-8 -*-
"""§3d, final piece: age-specific UNDER-ASCERTAINMENT. The severity (IHR) proxy
lifted the infant share but left adults too high. Here the detection fraction d_a
is anchored to an independent seroprevalence study (Italy; Infect Dis Ther 2023):
true infection was ~141x reported in ages 6-14 and ~3452x in ages >=15, i.e. the
reporting fraction collapses in adults. reported_a  ∝  infections_a × d_a.
Then compare the reported age distribution to documented US reported-case shares."""
import sys, os, json
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO = r"C:\Users\MarkRoutson\epydemix-dashboard"
OUT = os.path.join(REPO, "docs", "validation", "figs")
BANDS = ["0-1", "1-6", "7-10", "11-19", "20-49", "50-64", "65+"]
ADULT = [4, 5, 6]

# real-coverage modeled INFECTIONS distribution (from val_agedist4)
inf = np.array(json.load(open(os.path.join(OUT, "val3d_realcoverage.json")))["infections_real"]["dist"], float)

# age-specific reporting fraction d_a (reporting = 1 / under-ascertainment):
#   infant <1 : severe, hospitalised (IHR 0.30) -> near-complete detection   ~0.50 (assumption)
#   1-6       : symptomatic young children, frequently tested                 ~0.05 (assumption)
#   7-10      : Italy 6-14 band  -> 1/141                                       0.0071
#   11-19     : mix of 6-14 (1/141) and >=15 (1/3452)                          ~0.002
#   20-49/50-64/65+ : Italy >=15 band -> 1/3452                                0.00029
D = np.array([0.50, 0.05, 1/141, 0.002, 1/3452, 1/3452, 1/3452])

rep = inf * D
rep = rep / rep.sum()

def shares(d): return float(d[0]), float(sum(d[i] for i in ADULT))
inf_i, inf_a = shares(inf)
rep_i, rep_a = shares(rep)
print(f"{'':34s}{'infant<1yr':>12s}{'adult 20+':>12s}")
print(f"{'modeled infections (real cov.)':34s}{inf_i:11.1%}{inf_a:11.1%}")
print(f"{'reported (severity+ascertainment)':34s}{rep_i:11.1%}{rep_a:11.1%}")
print(f"{'documented reported cases':34s}{'5.6-12.5%':>12s}{'~19%':>12s}")
print("\nbands:        ", BANDS)
print("infections :  ", [round(float(x),3) for x in inf])
print("reported dist:", [round(float(x),3) for x in rep])

res = {"bands": BANDS, "reporting_fraction": {b: float(d) for b, d in zip(BANDS, D)},
       "infections": {"infant": inf_i, "adult": inf_a, "dist": [float(x) for x in inf]},
       "reported": {"infant": rep_i, "adult": rep_a, "dist": [float(x) for x in rep]},
       "documented": {"infant_national": 0.125, "infant_oregon": 0.056, "adult": 0.19},
       "source": "age gradient from Italy seroprevalence (Infect Dis Ther 2023): 6-14 ~141x, >=15 ~3452x under-ascertained"}
json.dump(res, open(os.path.join(OUT, "val3e_ascertainment.json"), "w"), indent=2)

fig, ax = plt.subplots(figsize=(8.6, 4.2))
x = np.arange(len(BANDS)); w = 0.38
ax.bar(x - w/2, inf, w, color="#94a3b8", label="Modeled infections (real coverage)")
ax.bar(x + w/2, rep, w, color="#0e7a78", label="Reported (severity + ascertainment)")
ax.set_xticks(x); ax.set_xticklabels(BANDS)
ax.set_ylabel("Share of cases"); ax.set_xlabel("Age band")
ax.set_title("Full reporting model: infections → reported cases by age", weight="bold", fontsize=11)
ax.annotate(f"infant {inf_i:.0%}→{rep_i:.0%} (doc ~12.5%)\nadult {inf_a:.0%}→{rep_a:.0%} (doc ~19%)",
            (2.5, max(rep)*0.8), fontsize=9, color="#182028")
ax.legend(fontsize=9, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val3e_ascertainment.png"), dpi=140)
print("\nwrote val3e_ascertainment.json and .png")
