# -*- coding: utf-8 -*-
"""Does the integrated immunization-by-age table close the §3d gap?
Rerun the age-distribution check with use_age_immunity = True and the
'Lane County (observed)' coverage profile, vs defaults (uniform immunity)."""
import sys, os, json
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
REPO = r"C:\Users\MarkRoutson\epydemix-dashboard"
sys.path.insert(0, REPO)
from constants import DEFAULT_AGE_GROUPS
from helpers import load_population, daily_doses_by_age
from engine.run import run_scenario
from engine import run as runmod
from schemas import MODEL_PARAM_SCHEMAS, INITIAL_CONDITION_DEFAULTS
from engine.calibration import modeled_case_summary
from data.coverage_store import get_coverage_profile, DEFAULT_PROFILE_NAME

mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS["SEIRS (Pertussis)"]}
ic0 = dict(INITIAL_CONDITION_DEFAULTS["SEIRS (Pertussis)"])
runmod.N_SIM = 10
ADULT = [i for i, a in enumerate(DEFAULT_AGE_GROUPS) if a in ("20-49", "50-64", "65+")]
prof = get_coverage_profile(DEFAULT_PROFILE_NAME)
print("coverage profile:", DEFAULT_PROFILE_NAME, "->", {a: prof[a] for a in DEFAULT_AGE_GROUPS})

GEO = "United_States__Oregon__Lane_County"
pop = load_population(GEO)

def run(use_age):
    doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=300, age_groups=DEFAULT_AGE_GROUPS, dt=0.5)
    cfg = {"model": "SEIRS (Pertussis)", "geography": GEO, "population": pop, "sim_length": 300,
           "time_step": 0.5, "initial_conditions": dict(ic0), "model_params": dict(mp0),
           "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
           "vaccination_campaigns": [], "daily_doses_by_age": doses, "observed_mode": "Off",
           "use_age_immunity": use_age,
           "age_immunity_pct": {a: float(prof[a]) for a in DEFAULT_AGE_GROUPS} if use_age else None}
    _c, trans, *_ = run_scenario(cfg)
    s = modeled_case_summary(trans)
    ad = s["age_dist"]
    return float(ad[0]), float(sum(ad[i] for i in ADULT)), [round(float(x), 3) for x in ad]

inf_def, ad_def, dist_def = run(False)
inf_age, ad_age, dist_age = run(True)
print(f"\n{'':22s}{'infant<1yr':>12s}{'adult 20+':>12s}")
print(f"{'uniform (defaults)':22s}{inf_def:11.1%}{ad_def:11.1%}")
print(f"{'age-immunity table':22s}{inf_age:11.1%}{ad_age:11.1%}")
print(f"{'documented':22s}{'5.6-12.5%':>12s}{'~19%':>12s}")
print("\nband dist (defaults):   ", dict(zip(DEFAULT_AGE_GROUPS, dist_def)))
print("band dist (age-immun):  ", dict(zip(DEFAULT_AGE_GROUPS, dist_age)))

res = {"profile": {a: prof[a] for a in DEFAULT_AGE_GROUPS},
       "uniform": {"infant": inf_def, "adult": ad_def, "dist": dist_def},
       "age_immunity": {"infant": inf_age, "adult": ad_age, "dist": dist_age},
       "documented": {"infant_oregon": 0.056, "infant_national": 0.125, "adult": 0.19}}
json.dump(res, open(os.path.join(REPO, "docs", "validation", "figs", "val3b_ageimmunity.json"), "w"), indent=2)
print("\nwrote val3b_ageimmunity.json")
