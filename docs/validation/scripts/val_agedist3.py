# -*- coding: utf-8 -*-
"""§3d refined: reconcile MODELED INFECTIONS with the documented REPORTED-CASE
age distribution via an age-specific detection/severity filter (the app's own
IHRs). reported_a  ∝  infections_a × detection_a.
  - infections        : model output (true infections)
  - severity-weighted : infections × IHR_a  (detection scales with severity)
Compare infant<1yr and adult 20+ shares to documented reported-case shares."""
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
prof = get_coverage_profile(DEFAULT_PROFILE_NAME)
GEO = "United_States__Oregon__Lane_County"
pop = load_population(GEO)

doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=300, age_groups=DEFAULT_AGE_GROUPS, dt=0.5)
cfg = {"model": "SEIRS (Pertussis)", "geography": GEO, "population": pop, "sim_length": 300,
       "time_step": 0.5, "initial_conditions": dict(ic0), "model_params": dict(mp0),
       "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
       "vaccination_campaigns": [], "daily_doses_by_age": doses, "observed_mode": "Off",
       "use_age_immunity": True,
       "age_immunity_pct": {a: float(prof[a]) for a in DEFAULT_AGE_GROUPS}}
_c, trans, *_ = run_scenario(cfg)
inf = np.array(modeled_case_summary(trans)["age_dist"], float)   # true infections
sev = inf * IHR; sev = sev / sev.sum()                           # severity-weighted (≈ detected)

def shares(d): return float(d[0]), float(sum(d[i] for i in ADULT))
inf_i, inf_a = shares(inf)
sev_i, sev_a = shares(sev)
DOC = {"infant_national": 0.125, "infant_oregon": 0.056, "adult": 0.19}

print(f"{'':26s}{'infant<1yr':>12s}{'adult 20+':>12s}")
print(f"{'modeled infections':26s}{inf_i:11.1%}{inf_a:11.1%}")
print(f"{'severity-weighted (reported)':26s}{sev_i:11.1%}{sev_a:11.1%}")
print(f"{'documented reported cases':26s}{'5.6-12.5%':>12s}{'~19%':>12s}")
print("\nbands:", DEFAULT_AGE_GROUPS)
print("infections :", [round(float(x),3) for x in inf])
print("sev-weighted:", [round(float(x),3) for x in sev])

res = {"bands": DEFAULT_AGE_GROUPS, "ihr": {a: _DEFAULT_IHR[a] for a in DEFAULT_AGE_GROUPS},
       "infections": {"dist": [float(x) for x in inf], "infant": inf_i, "adult": inf_a},
       "severity_weighted": {"dist": [float(x) for x in sev], "infant": sev_i, "adult": sev_a},
       "documented": DOC}
json.dump(res, open(os.path.join(OUT, "val3c_severity.json"), "w"), indent=2)

# figure: infections vs severity-weighted across bands + infant/adult callouts
fig, ax = plt.subplots(figsize=(8.6, 4.2))
x = np.arange(len(DEFAULT_AGE_GROUPS)); w = 0.38
ax.bar(x - w/2, inf, w, color="#94a3b8", label="Modeled infections")
ax.bar(x + w/2, sev, w, color="#0e7a78", label="Severity-weighted (≈ reported)")
ax.set_xticks(x); ax.set_xticklabels(DEFAULT_AGE_GROUPS)
ax.set_ylabel("Share of cases"); ax.set_xlabel("Age band")
ax.set_title("Reconciling modeled infections with reported cases via age-specific detection",
             weight="bold", fontsize=11)
ax.annotate(f"infant: {inf_i:.0%} → {sev_i:.0%}\n(documented ~12.5%)", (0, max(inf[0], sev[0])),
            textcoords="offset points", xytext=(4, 24), fontsize=8.5, color="#182028")
ax.legend(fontsize=9, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val3c_severity.png"), dpi=140)
print("\nwrote val3c_severity.json and .png")
