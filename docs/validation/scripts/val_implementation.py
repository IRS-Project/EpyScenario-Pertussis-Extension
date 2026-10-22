# -*- coding: utf-8 -*-
"""VALIDATION Part 1 - implementation correctness.

  (a) Population conservation: sum of the 8 compartments == N every day.
  (b) Non-negativity: no compartment goes negative.
  (c) R0 recovery: in the clean naive-SEIR limit (no prior immunity, no
      seasonality, empty partial track) the realized early growth rate r must
      imply the R0 that was fed in, via R0 = (1 + r/gamma)(1 + r/mu).
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

OUT = os.path.join(REPO, "docs", "validation", "figs")
os.makedirs(OUT, exist_ok=True)
MODEL = "SEIRS (Pertussis)"
GEO = "United_States__Oregon__Lane_County"
pop = load_population(GEO)
N = float(pop.Nk.sum())
COMPS = ["S", "E", "I", "R", "Sp", "Ep", "Ip", "Rp"]

mp0 = {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS[MODEL]}


def make_cfg(mp, ic, sim=250, dt=0.5):
    doses = daily_doses_by_age(campaigns=[], Nk=pop.Nk, sim_length=sim,
                               age_groups=DEFAULT_AGE_GROUPS, dt=dt)
    return {"model": MODEL, "geography": GEO, "population": pop, "sim_length": sim,
            "time_step": dt, "initial_conditions": ic, "model_params": mp,
            "contact_interventions": [], "vaccination_settings": {"target_compartments": ["S"]},
            "vaccination_campaigns": [], "daily_doses_by_age": doses,
            "observed_mode": "Off", "use_age_immunity": False}


results = {"conservation": {}, "nonnegativity": {}, "r0_recovery": []}
runmod.N_SIM = 20  # smooth stochastic mean

# ---- (a)+(b) conservation & non-negativity on a default run ----
ic_def = dict(INITIAL_CONDITION_DEFAULTS[MODEL])
df = run_scenario(make_cfg(dict(mp0), ic_def))[0]
tot = sum(df[f"{c}_total"].to_numpy() for c in COMPS)
rel_err = np.abs(tot - N) / N
min_comp = min(float(df[f"{c}_total"].to_numpy().min()) for c in COMPS)
results["conservation"] = {"N": N, "max_rel_error": float(rel_err.max()),
                           "mean_rel_error": float(rel_err.mean())}
results["nonnegativity"] = {"min_compartment_value": min_comp,
                            "ok": bool(min_comp >= -1e-6 * N)}
print(f"[a] conservation: max rel error {rel_err.max():.2e} (N={N:,.0f})")
print(f"[b] non-negativity: min compartment value {min_comp:.3g} -> ok={results['nonnegativity']['ok']}")

# ---- (c) R0 recovery in clean naive-SEIR limit ----
# no prior immunity, no partial seeding, seasonality off; tiny seed.
ic_clean = {"infected_pct": 0.05, "immune_pct": 0.0,
            "partial_infection_pct": 0.0, "partial_immune_pct": 0.0}
incub = mp0["incubation_period"]; infper = mp0["infectious_period"]
gamma = 1.0 / incub   # E->I rate
mu = 1.0 / infper     # I->R rate
runmod.N_SIM = 30
for R0_in in [6.0, 8.0, 12.0]:
    mp = dict(mp0); mp["R0"] = R0_in; mp["seasonality_amplitude"] = "None"
    d = run_scenario(make_cfg(mp, dict(ic_clean), sim=120, dt=0.5))[0]
    t = d["t"].to_numpy()
    inf = (d["I_total"] + d["Ip_total"] + d["E_total"] + d["Ep_total"]).to_numpy()
    cum_removed = d["R_total"].to_numpy() + d["Rp_total"].to_numpy()
    # early exponential window: susceptible depletion < 2% of N
    depleted = (cum_removed + inf) / N
    mask = (depleted < 0.02) & (inf > 0)
    # need a few points; take the earliest contiguous window
    idx = np.where(mask)[0]
    if len(idx) < 5:
        idx = np.arange(5, 40)
    w = idx[(idx >= 5)][:40]  # skip the very first noisy days
    r = np.polyfit(t[w], np.log(inf[w]), 1)[0]
    R0_implied = (1 + r / gamma) * (1 + r / mu)
    results["r0_recovery"].append({"R0_input": R0_in, "growth_rate_per_day": float(r),
                                   "R0_implied": float(R0_implied),
                                   "rel_error": float(abs(R0_implied - R0_in) / R0_in)})
    print(f"[c] R0 in={R0_in:5.1f}  r={r:.4f}/day  R0 implied={R0_implied:5.2f}  "
          f"rel err={abs(R0_implied-R0_in)/R0_in:.1%}")

with open(os.path.join(OUT, "val1_implementation.json"), "w") as f:
    json.dump(results, f, indent=2)

# ---- figures ----
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(df["t"], rel_err, color="#0e7a78", lw=1.6)
ax[0].set_title("Population conservation", weight="bold", fontsize=12)
ax[0].set_xlabel("Day"); ax[0].set_ylabel("|Σcompartments − N| / N")
ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0)); ax[0].grid(alpha=.3)
for s in ("top", "right"): ax[0].spines[s].set_visible(False)

rr = results["r0_recovery"]
xi = [r["R0_input"] for r in rr]; yi = [r["R0_implied"] for r in rr]
lim = [5, 13]
ax[1].plot(lim, lim, ls="--", color="#999", label="identity (perfect recovery)")
ax[1].scatter(xi, yi, color="#c2410c", s=70, zorder=3, label="engine")
for r in rr:
    ax[1].annotate(f"{r['rel_error']:.0%}", (r["R0_input"], r["R0_implied"]),
                   textcoords="offset points", xytext=(8, -4), fontsize=9)
ax[1].set_title("R₀ recovery (clean naive-SEIR limit)", weight="bold", fontsize=12)
ax[1].set_xlabel("R₀ set in model"); ax[1].set_ylabel("R₀ implied by realized growth")
ax[1].set_xlim(lim); ax[1].set_ylim(lim); ax[1].grid(alpha=.3); ax[1].legend(fontsize=9, frameon=False)
for s in ("top", "right"): ax[1].spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val1_implementation.png"), dpi=140)
print("\nwrote", os.path.join(OUT, "val1_implementation.json"), "and .png")
print("VAL1_DONE")
