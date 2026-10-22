# -*- coding: utf-8 -*-
"""§3d, non-circular coverage: replace the 'vaccinated-share-of-cases' profile with
REAL population coverage by age (CDC NIS-Child DTaP, NIS-Teen Tdap, BRFSS adult
Tdap), pulled from data.cdc.gov this session. Re-run infections + severity-weighted
age distribution and compare to the circular profile and documented reported cases.

Oregon coverage (pertussis-protective doses), with source + band mapping:
  0-1   DTaP primary series building over the first year (NIS >=4 by 19mo 69.8%,
        partial-year protection)                                   -> 35%
  1-6   DTaP >=4 doses, full childhood series (NIS OR 24mo 78.7 / 35mo 84.2)  -> 80%
  7-10  DTaP series persists before the adolescent booster                   -> 82%
  11-19 Tdap adolescent booster (NIS-Teen OR 88.0)                           -> 88%
  20-49 adult Tdap (BRFSS aetd-68ew, OR 18-64 ~53.9)                         -> 54%
  50-64 adult Tdap, waning                                                   -> 45%
  65+   adult Tdap, lowest                                                   -> 25%
"""
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
from helpers import load_population, daily_doses_by_age
from engine.run import run_scenario
from engine import run as runmod
from schemas import MODEL_PARAM_SCHEMAS, INITIAL_CONDITION_DEFAULTS
from engine.calibration import modeled_case_summary
from engine.hospitalization import _DEFAULT_IHR
from data.coverage_store import get_coverage_profile, DEFAULT_PROFILE_NAME

OUT = os.path.join(REPO, "docs", "validation", "figs")
mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS["SEIRS (Pertussis)"]}
ic0 = dict(INITIAL_CONDITION_DEFAULTS["SEIRS (Pertussis)"])
runmod.N_SIM = 12
ADULT = [i for i, a in enumerate(DEFAULT_AGE_GROUPS) if a in ("20-49", "50-64", "65+")]
IHR = np.array([_DEFAULT_IHR[a] for a in DEFAULT_AGE_GROUPS])

REAL = {"0-1": 35.0, "1-6": 80.0, "7-10": 82.0, "11-19": 88.0,
        "20-49": 54.0, "50-64": 45.0, "65+": 25.0}  # CDC NIS-Child/Teen + BRFSS (Oregon)
CIRC = get_coverage_profile(DEFAULT_PROFILE_NAME)    # vaccinated-share-of-cases (circular)

GEO = "United_States__Oregon__Lane_County"
pop = load_population(GEO)

def run(profile):
    doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=300, age_groups=DEFAULT_AGE_GROUPS, dt=0.5)
    cfg = {"model": "SEIRS (Pertussis)", "geography": GEO, "population": pop, "sim_length": 300,
           "time_step": 0.5, "initial_conditions": dict(ic0), "model_params": dict(mp0),
           "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
           "vaccination_campaigns": [], "daily_doses_by_age": doses, "observed_mode": "Off",
           "use_age_immunity": True, "age_immunity_pct": {a: float(profile[a]) for a in DEFAULT_AGE_GROUPS}}
    _c, trans, *_ = run_scenario(cfg)
    inf = np.array(modeled_case_summary(trans)["age_dist"], float)
    sev = inf * IHR; sev = sev / sev.sum()
    return inf, sev

def shares(d): return float(d[0]), float(sum(d[i] for i in ADULT))

inf_c, sev_c = run(CIRC)
inf_r, sev_r = run(REAL)
print("profile (circular, share-of-cases):", {a: round(CIRC[a],1) for a in DEFAULT_AGE_GROUPS})
print("profile (REAL, NIS/BRFSS):         ", REAL)
print(f"\n{'':34s}{'infant<1yr':>12s}{'adult 20+':>12s}")
print(f"{'infections — circular profile':34s}{shares(inf_c)[0]:11.1%}{shares(inf_c)[1]:11.1%}")
print(f"{'infections — REAL coverage':34s}{shares(inf_r)[0]:11.1%}{shares(inf_r)[1]:11.1%}")
print(f"{'reported (sev-wt) — REAL coverage':34s}{shares(sev_r)[0]:11.1%}{shares(sev_r)[1]:11.1%}")
print(f"{'documented reported cases':34s}{'5.6-12.5%':>12s}{'~19%':>12s}")

res = {"bands": DEFAULT_AGE_GROUPS, "profile_circular": {a: CIRC[a] for a in DEFAULT_AGE_GROUPS},
       "profile_real": REAL,
       "infections_circular": {"infant": shares(inf_c)[0], "adult": shares(inf_c)[1]},
       "infections_real": {"infant": shares(inf_r)[0], "adult": shares(inf_r)[1], "dist": [float(x) for x in inf_r]},
       "reported_real": {"infant": shares(sev_r)[0], "adult": shares(sev_r)[1], "dist": [float(x) for x in sev_r]},
       "documented": {"infant_national": 0.125, "infant_oregon": 0.056, "adult": 0.19},
       "sources": "CDC NIS-Child DTaP (fhky-rtsk), NIS-Teen Tdap (ee48-w5t6), BRFSS adult Tdap (aetd-68ew)"}
json.dump(res, open(os.path.join(OUT, "val3d_realcoverage.json"), "w"), indent=2)

# figure: coverage profiles + resulting infant share
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
x = np.arange(len(DEFAULT_AGE_GROUPS)); w = 0.38
a1.bar(x - w/2, [CIRC[a] for a in DEFAULT_AGE_GROUPS], w, color="#c2410c", label="Circular (share-of-cases)")
a1.bar(x + w/2, [REAL[a] for a in DEFAULT_AGE_GROUPS], w, color="#0e7a78", label="Real (NIS/BRFSS)")
a1.set_xticks(x); a1.set_xticklabels(DEFAULT_AGE_GROUPS, fontsize=8); a1.set_ylabel("% immune / covered")
a1.set_title("Coverage-by-age input: circular vs real", weight="bold", fontsize=11)
a1.legend(fontsize=8, frameon=False); a1.grid(axis="y", alpha=.3)
for s in ("top", "right"): a1.spines[s].set_visible(False)

labels = ["infections\n(real cov.)", "reported\n(sev-wt, real)", "documented"]
inf_sh = [shares(inf_r)[0], shares(sev_r)[0], 0.125]
ad_sh = [shares(inf_r)[1], shares(sev_r)[1], 0.19]
xx = np.arange(3)
a2.bar(xx - w/2, inf_sh, w, color="#c2410c", label="infant <1yr")
a2.bar(xx + w/2, ad_sh, w, color="#0e7a78", label="adult 20+")
a2.set_xticks(xx); a2.set_xticklabels(labels, fontsize=8); a2.set_ylabel("share of cases")
a2.set_title("Age shares under real coverage", weight="bold", fontsize=11)
a2.legend(fontsize=8, frameon=False); a2.grid(axis="y", alpha=.3)
for s in ("top", "right"): a2.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val3d_realcoverage.png"), dpi=140)
print("\nwrote val3d_realcoverage.json and .png")
