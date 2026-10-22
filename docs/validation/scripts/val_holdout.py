# -*- coding: utf-8 -*-
"""VALIDATION Part 2 - out-of-sample temporal hold-out.

For each state: fit R0 (1-D grid, transparent) to one NNDSS season's weekly
curve, then SCORE the model on the *untouched* next season. Compare against:
  - in-sample floor  : fit the test year, score the test year (best case)
  - persistence base : use the prior year's observed curve to predict the test year

Metric: peak-normalised weekly-incidence shape distance (RMSE + peak-week offset),
scale-free (absolute magnitude is not expected to generalise; shape/timing is).
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
from data.weekly_sources import load_weekly_pertussis, weekly_series_to_array
from engine.calibration import weekly_incidence_from_trans

OUT = os.path.join(REPO, "docs", "validation", "figs")
os.makedirs(OUT, exist_ok=True)
MODEL = "SEIRS (Pertussis)"
SIM = 371; DT = 0.5
R0_GRID = np.arange(6.0, 18.01, 0.5)
mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS[MODEL]}
ic0 = dict(INITIAL_CONDITION_DEFAULTS[MODEL])
runmod.N_SIM = 6

STATES = [("Oregon", "United_States__Oregon"), ("Alabama", "United_States__Alabama")]
YEAR_FIT, YEAR_TEST = 2023, 2024


def shape_dist(obs, mod):
    """Peak-normalised shape distance between two weekly series."""
    L = min(len(obs), len(mod))
    a = np.asarray(obs[:L], float); b = np.asarray(mod[:L], float)
    if a.max() <= 0 or b.max() <= 0:
        return None, None
    an = a / a.max(); bn = b / b.max()
    rmse = float(np.sqrt(np.mean((an - bn) ** 2)))
    peak_off = int(abs(int(np.argmax(a)) - int(np.argmax(b))))
    return rmse + 0.5 * peak_off / L, {"rmse": rmse, "peak_week_offset": peak_off}


def modeled_weekly(geo_pop, R0):
    mp = dict(mp0); mp["R0"] = float(R0)
    doses = daily_doses_by_age(campaigns=[], Nk=geo_pop.Nk, sim_length=SIM,
                               age_groups=DEFAULT_AGE_GROUPS, dt=DT)
    cfg = {"model": MODEL, "geography": "x", "population": geo_pop, "sim_length": SIM,
           "time_step": DT, "initial_conditions": dict(ic0), "model_params": mp,
           "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
           "vaccination_campaigns": [], "daily_doses_by_age": doses,
           "observed_mode": "Off", "use_age_immunity": False}
    _comp, trans, *_ = run_scenario(cfg)
    return np.asarray(weekly_incidence_from_trans(trans, MODEL), float)


report = {}
for name, geo in STATES:
    print(f"\n=== {name} ===", flush=True)
    try:
        df_fit, src_fit = load_weekly_pertussis(name, years=[YEAR_FIT])
        df_test, src_test = load_weekly_pertussis(name, years=[YEAR_TEST])
        obs_fit = np.asarray(weekly_series_to_array(df_fit), float)
        obs_test = np.asarray(weekly_series_to_array(df_test), float)
    except Exception as e:
        print("  data load FAILED:", e, flush=True)
        report[name] = {"error": str(e)}
        continue
    print(f"  {YEAR_FIT}: {len(obs_fit)} wks / {int(obs_fit.sum())} cases ({src_fit}) | "
          f"{YEAR_TEST}: {len(obs_test)} wks / {int(obs_test.sum())} cases ({src_test})", flush=True)

    pop = load_population(geo)
    # simulate the R0 grid once, store modeled weekly curves
    grid = {}
    for R0 in R0_GRID:
        grid[round(float(R0), 2)] = modeled_weekly(pop, R0)

    def best_R0(obs):
        best, bd = None, np.inf
        for R0, mod in grid.items():
            d, _ = shape_dist(obs, mod)
            if d is not None and d < bd:
                bd, best = d, R0
        return best, bd

    R0_fit, d_fit_insample = best_R0(obs_fit)     # fit on 2023
    R0_test, d_test_floor = best_R0(obs_test)     # in-sample floor on 2024
    d_holdout, det_h = shape_dist(obs_test, grid[R0_fit])   # 2023-fit model vs 2024
    d_persist, det_p = shape_dist(obs_test, obs_fit)        # persistence baseline

    report[name] = {
        "fit_year": YEAR_FIT, "test_year": YEAR_TEST,
        "R0_fit_on_%d" % YEAR_FIT: R0_fit, "R0_insample_%d" % YEAR_TEST: R0_test,
        "holdout_distance": d_holdout, "insample_floor": d_test_floor,
        "persistence_baseline": d_persist,
        "holdout_detail": det_h, "persistence_detail": det_p,
        "beats_persistence": bool(d_holdout < d_persist),
    }
    print(f"  R0 fit on {YEAR_FIT} = {R0_fit} | in-sample-best R0 for {YEAR_TEST} = {R0_test}")
    print(f"  hold-out dist (fit {YEAR_FIT} -> score {YEAR_TEST}) = {d_holdout:.3f}  (peak-wk off {det_h['peak_week_offset']})")
    print(f"  in-sample floor ({YEAR_TEST})                      = {d_test_floor:.3f}")
    print(f"  persistence baseline ({YEAR_FIT}->{YEAR_TEST})          = {d_persist:.3f}")
    print(f"  model beats persistence: {d_holdout < d_persist}")

    # figure: observed test-year vs model(R0_fit) vs persistence, peak-normalised
    mod = grid[R0_fit]; L = min(len(obs_test), len(mod))
    norm = lambda x: np.asarray(x[:L], float) / (max(np.asarray(x[:L]).max(), 1e-9))
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    ax.plot(norm(obs_test), color="#182028", lw=2.4, label=f"Observed {YEAR_TEST} (held out)")
    ax.plot(norm(mod), color="#0e7a78", lw=2.2, label=f"Model (R₀ fit on {YEAR_FIT})")
    ax.plot(norm(obs_fit), color="#c2410c", lw=1.6, ls="--", label=f"Persistence ({YEAR_FIT} curve)")
    ax.set_title(f"{name} — out-of-sample hold-out (peak-normalised weekly)", weight="bold", fontsize=12)
    ax.set_xlabel("Week index"); ax.set_ylabel("Weekly cases (normalised)")
    ax.grid(alpha=.3); ax.legend(fontsize=9, frameon=False)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"val2_holdout_{name.lower()}.png"), dpi=140); plt.close(fig)

with open(os.path.join(OUT, "val2_holdout.json"), "w") as f:
    json.dump(report, f, indent=2)
print("\nwrote val2_holdout.json")
print("VAL2_DONE")
