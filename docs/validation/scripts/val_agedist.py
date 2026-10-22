# -*- coding: utf-8 -*-
"""VALIDATION Part 3d - age-distribution face validity.

Compare the model's MODELED case age-distribution (infant <1yr share, adult 20+
share) against documented US pertussis epidemiology encoded in the EPICC
Omaha-2008 validation model / MMWR 60(1) 2011 and the EPICC Oregon model:
    infant share of cases : 0.125 (national) / 0.056 (Oregon)
    adult (20+) share      : 0.19
Run under model defaults (uniform immunity) to test whether the structure alone
reproduces the documented pattern, or whether the age-structured-immunity inputs
are needed (an identifiability/limitation statement either way).
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

OUT = os.path.join(REPO, "docs", "validation", "figs")
os.makedirs(OUT, exist_ok=True)
MODEL = "SEIRS (Pertussis)"
mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS[MODEL]}
ic0 = dict(INITIAL_CONDITION_DEFAULTS[MODEL])
runmod.N_SIM = 10

# documented targets (EPICC pertussis models / MMWR 60(1) 2011)
DOC = {"infant_national": 0.125, "infant_oregon": 0.056, "adult": 0.19}
print("age bands:", DEFAULT_AGE_GROUPS)
ADULT_IDX = [i for i, a in enumerate(DEFAULT_AGE_GROUPS) if a in ("20-49", "50-64", "65+")]

rows = {}
for name, geo in [("Oregon (Lane County)", "United_States__Oregon__Lane_County"),
                  ("Alabama", "United_States__Alabama")]:
    pop = load_population(geo)
    doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=300,
                               age_groups=DEFAULT_AGE_GROUPS, dt=0.5)
    cfg = {"model": MODEL, "geography": geo, "population": pop, "sim_length": 300,
           "time_step": 0.5, "initial_conditions": dict(ic0), "model_params": dict(mp0),
           "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
           "vaccination_campaigns": [], "daily_doses_by_age": doses,
           "observed_mode": "Off", "use_age_immunity": False}
    _c, trans, *_ = run_scenario(cfg)
    summ = modeled_case_summary(trans)
    ad = summ["age_dist"]
    infant = float(ad[0]); adult = float(sum(ad[i] for i in ADULT_IDX))
    rows[name] = {"age_dist": [float(x) for x in ad], "infant_share": infant,
                  "adult_share": adult, "overall_partial_share": float(summ["overall_partial_share"])}
    print(f"{name:22s} infant<1yr share {infant:.1%} | adult 20+ share {adult:.1%} | "
          f"partial(vaccinated)-track share {summ['overall_partial_share']:.1%}")

print("\nDocumented: infant 12.5% (national) / 5.6% (Oregon); adult 19%")
res = {"documented": DOC, "modeled": rows, "age_bands": DEFAULT_AGE_GROUPS}
with open(os.path.join(OUT, "val3_agedist.json"), "w") as f:
    json.dump(res, f, indent=2)

# figure: modeled infant & adult shares vs documented
fig, ax = plt.subplots(figsize=(8.2, 4.0))
labels = list(rows.keys())
x = np.arange(len(labels)); w = 0.35
inf = [rows[k]["infant_share"] for k in labels]
adu = [rows[k]["adult_share"] for k in labels]
ax.bar(x - w/2, inf, w, color="#c2410c", label="Modeled infant <1yr")
ax.bar(x + w/2, adu, w, color="#0e7a78", label="Modeled adult 20+")
ax.axhline(DOC["infant_national"], color="#c2410c", ls="--", lw=1.2, label="Documented infant (natl 12.5%)")
ax.axhline(DOC["adult"], color="#0e7a78", ls="--", lw=1.2, label="Documented adult (19%)")
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_ylabel("Share of modeled cases"); ax.set_ylim(0, max(max(inf+adu)*1.3, 0.35))
ax.set_title("Modeled case age-distribution vs documented pertussis epidemiology",
             weight="bold", fontsize=11)
ax.legend(fontsize=8, frameon=False, ncol=2);
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val3_agedist.png"), dpi=140)
print("wrote val3_agedist.json and .png")
print("VAL3_AGEDIST_DONE")
