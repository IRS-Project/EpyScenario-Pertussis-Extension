# -*- coding: utf-8 -*-
"""§6 steps 1 & 3: age-structured demographic model -> realistic infection age
distribution, and a reporting gradient calibrated to local (Lane County) data."""
import sys, os, json
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO = r"C:\Users\MarkRoutson\epydemix-dashboard"
sys.path.insert(0, REPO)
from constants import DEFAULT_AGE_GROUPS
from helpers import load_population
from engine.demography import build_age_params, infection_age_distribution
from engine.hospitalization import _DEFAULT_IHR
from data.observed_datasets import aggregate_to_bands, LANE_COUNTY_RAW

OUT = os.path.join(REPO, "docs", "validation", "figs")
BANDS = DEFAULT_AGE_GROUPS
ADULT = [4, 5, 6]
IHR = np.array([_DEFAULT_IHR[a] for a in BANDS])

pop = load_population("United_States__Oregon")
M = sum(pop.contact_matrices[l] for l in pop.contact_matrices)
inf = np.array(infection_age_distribution(build_age_params(M))["dist"], float)

# single-season (no births) infections, for contrast
try:
    ss = np.array(json.load(open(os.path.join(OUT, "val3d_realcoverage.json")))["infections_real"]["dist"], float)
except Exception:
    ss = None

# severity-weighted "reported" from the age-structured infections
sev = inf * IHR; sev = sev/sev.sum()

# Lane County observed reported distribution + implied reporting gradient
b = aggregate_to_bands(LANE_COUNTY_RAW)
obs = np.array([b[a]["naive"]+b[a]["partial"]+b[a].get("Unknown", 0) for a in BANDS], float)
obs_dist = obs/obs.sum()
with np.errstate(divide="ignore", invalid="ignore"):
    d_raw = np.where(obs_dist > 0, obs_dist/inf, np.nan)
d_norm = d_raw/np.nanmax(d_raw)

def shares(x): return float(x[0]), float(sum(x[i] for i in ADULT))
print("bands:        ", BANDS)
print("infections (age-structured):", np.round(inf, 4), f"-> infant {shares(inf)[0]:.1%}, adult {shares(inf)[1]:.1%}")
if ss is not None:
    print("infections (single-season): ", np.round(ss, 4), f"-> infant {shares(ss)[0]:.1%}, adult {shares(ss)[1]:.1%}")
print("severity-weighted reported: ", np.round(sev, 4), f"-> infant {shares(sev)[0]:.1%}, adult {shares(sev)[1]:.1%}")
print("Lane County observed:       ", np.round(obs_dist, 4))
print("implied reporting gradient: ", np.round(d_norm, 3))

res = {"bands": BANDS,
       "infections_age_structured": {"dist": [float(x) for x in inf], "infant": shares(inf)[0], "adult": shares(inf)[1]},
       "infections_single_season": ({"infant": shares(ss)[0], "adult": shares(ss)[1]} if ss is not None else None),
       "severity_weighted": {"dist": [float(x) for x in sev], "infant": shares(sev)[0], "adult": shares(sev)[1]},
       "observed_lane_county": {"dist": [float(x) for x in obs_dist]},
       "reporting_gradient_norm": [None if np.isnan(v) else round(float(v), 3) for v in d_norm],
       "documented": {"infant_national": 0.125, "infant_oregon": 0.056, "adult": 0.19}}
json.dump(res, open(os.path.join(OUT, "val6c_agestructured.json"), "w"), indent=2)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2))
x = np.arange(len(BANDS)); w = 0.4
if ss is not None:
    a1.bar(x - w/2, ss, w, color="#c2410c", label="Single-season (no births)")
a1.bar(x + w/2, inf, w, color="#0e7a78", label="Age-structured + demography")
a1.set_xticks(x); a1.set_xticklabels(BANDS, fontsize=8); a1.set_ylabel("Share of infections")
a1.set_title("Infection age distribution: births raise the infant share", weight="bold", fontsize=10.5)
a1.legend(fontsize=8, frameon=False); a1.grid(axis="y", alpha=.3)
for s in ("top", "right"): a1.spines[s].set_visible(False)
a1.annotate(f"infant {shares(ss)[0]:.1%}→{shares(inf)[0]:.1%}" if ss is not None else "",
            (0, inf[0]), textcoords="offset points", xytext=(2, 14), fontsize=8)

a2.bar(x, d_norm, color="#185fa5")
a2.set_xticks(x); a2.set_xticklabels(BANDS, fontsize=8); a2.set_ylabel("Relative reporting fraction")
a2.set_title("Reporting gradient implied by Lane County", weight="bold", fontsize=10.5)
a2.grid(axis="y", alpha=.3)
for s in ("top", "right"): a2.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val6c_agestructured.png"), dpi=140)
print("\nwrote val6c_agestructured.json and .png")
