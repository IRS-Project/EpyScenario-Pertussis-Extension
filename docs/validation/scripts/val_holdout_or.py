# -*- coding: utf-8 -*-
"""Improved hold-out: Oregon OUTBREAK-to-OUTBREAK (fit 2024 -> score 2025), so the
fit year actually has an epidemic signal (unlike 2023's 12 cases). Appends to
val2_holdout.json under 'Oregon (2024->2025)'."""
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
from data.weekly_sources import load_weekly_pertussis, weekly_series_to_array
from engine.calibration import weekly_incidence_from_trans

OUT = os.path.join(REPO, "docs", "validation", "figs")
MODEL = "SEIRS (Pertussis)"; SIM = 371; DT = 0.5
R0_GRID = np.arange(6.0, 18.01, 0.5)
mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS[MODEL]}
ic0 = dict(INITIAL_CONDITION_DEFAULTS[MODEL])
runmod.N_SIM = 6
GEO = "United_States__Oregon"; NAME = "Oregon (2024->2025)"
YEAR_FIT, YEAR_TEST = 2024, 2025


def shape_dist(obs, mod):
    L = min(len(obs), len(mod))
    a = np.asarray(obs[:L], float); b = np.asarray(mod[:L], float)
    if a.max() <= 0 or b.max() <= 0: return None, None
    an, bn = a / a.max(), b / b.max()
    rmse = float(np.sqrt(np.mean((an - bn) ** 2)))
    off = int(abs(int(np.argmax(a)) - int(np.argmax(b))))
    return rmse + 0.5 * off / L, {"rmse": rmse, "peak_week_offset": off}


pop = load_population(GEO)
def modeled_weekly(R0):
    mp = dict(mp0); mp["R0"] = float(R0)
    doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=SIM, age_groups=DEFAULT_AGE_GROUPS, dt=DT)
    cfg = {"model": MODEL, "geography": "x", "population": pop, "sim_length": SIM, "time_step": DT,
           "initial_conditions": dict(ic0), "model_params": mp, "contact_interventions": [],
           "vaccination_settings": {"target_compartments": ["S"]}, "vaccination_campaigns": [],
           "daily_doses_by_age": doses, "observed_mode": "Off", "use_age_immunity": False}
    _c, trans, *_ = run_scenario(cfg)
    return np.asarray(weekly_incidence_from_trans(trans, MODEL), float)

obs_fit = np.asarray(weekly_series_to_array(load_weekly_pertussis("Oregon", years=[YEAR_FIT])[0]), float)
obs_test = np.asarray(weekly_series_to_array(load_weekly_pertussis("Oregon", years=[YEAR_TEST])[0]), float)
print(f"{YEAR_FIT}: {int(obs_fit.sum())} cases | {YEAR_TEST}: {int(obs_test.sum())} cases", flush=True)

grid = {round(float(R0), 2): modeled_weekly(R0) for R0 in R0_GRID}
def best(obs):
    b, bd = None, np.inf
    for R0, m in grid.items():
        d, _ = shape_dist(obs, m)
        if d is not None and d < bd: bd, b = d, R0
    return b, bd

R0_fit, _ = best(obs_fit)
R0_test, d_floor = best(obs_test)
d_hold, det_h = shape_dist(obs_test, grid[R0_fit])
d_pers, det_p = shape_dist(obs_test, obs_fit)
row = {"fit_year": YEAR_FIT, "test_year": YEAR_TEST,
       "R0_fit_on_%d" % YEAR_FIT: R0_fit, "R0_insample_%d" % YEAR_TEST: R0_test,
       "holdout_distance": d_hold, "insample_floor": d_floor, "persistence_baseline": d_pers,
       "holdout_detail": det_h, "persistence_detail": det_p, "beats_persistence": bool(d_hold < d_pers)}
print(f"R0 fit {YEAR_FIT}={R0_fit} | in-sample {YEAR_TEST}={R0_test}")
print(f"hold-out {d_hold:.3f} (peak off {det_h['peak_week_offset']}) | floor {d_floor:.3f} | persistence {d_pers:.3f} | beats={row['beats_persistence']}")

jf = os.path.join(OUT, "val2_holdout.json")
data = json.load(open(jf)) if os.path.exists(jf) else {}
data[NAME] = row
json.dump(data, open(jf, "w"), indent=2)

mod = grid[R0_fit]; L = min(len(obs_test), len(mod))
norm = lambda x: np.asarray(x[:L], float) / max(np.asarray(x[:L]).max(), 1e-9)
fig, ax = plt.subplots(figsize=(8.2, 4.0))
ax.plot(norm(obs_test), color="#182028", lw=2.4, label=f"Observed {YEAR_TEST} (held out)")
ax.plot(norm(mod), color="#0e7a78", lw=2.2, label=f"Model (R₀ fit on {YEAR_FIT})")
ax.plot(norm(obs_fit), color="#c2410c", lw=1.6, ls="--", label=f"Persistence ({YEAR_FIT})")
ax.set_title("Oregon — outbreak-to-outbreak hold-out 2024→2025 (peak-normalised)", weight="bold", fontsize=11)
ax.set_xlabel("Week index"); ax.set_ylabel("Weekly cases (normalised)")
ax.grid(alpha=.3); ax.legend(fontsize=9, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val2_holdout_oregon_outbreak.png"), dpi=140)
print("VAL2_OR_DONE")
