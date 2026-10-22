# -*- coding: utf-8 -*-
"""VALIDATION §6 (demography): the Wearing-Rohani structural test.

A standalone DETERMINISTIC version of the app's 8-compartment partial-immunity
pertussis model, with vital dynamics added (births into S, deaths from all
compartments) — the ingredient the single-season app model lacks. Integrated over
~160 years; the intrinsic INTER-EPIDEMIC period is read from the annual-incidence
spectrum and compared to the literature (2-4.6 yr; 3.3-4.1 yr post-vaccination).

Parameters are the app's SEIRS (Pertussis) defaults. Birth rate = US crude birth
rate 11.6/1000/yr (NCHS e6fc-ccez, 2018). Mild seasonal forcing sustains the
oscillation; annual aggregation isolates the multi-year cycle.
"""
import sys, os, json
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\MarkRoutson\epydemix-dashboard\docs\validation\figs"
os.makedirs(OUT, exist_ok=True)

# ---- parameters (app SEIRS Pertussis defaults) ----
R0 = 8.0
gamma = 1/15.0      # I->R  (naive infectious period 15 d)
gamma_p = 1/12.0    # Ip->Rp (partial infectious 12 d)
eps = 1/9.0         # E->I, Ep->Ip (incubation 9 d)
sigma = 0.5         # rel. infectiousness partial
delta = 0.3         # rel. susceptibility partial
w1 = 1/(4*365.0)    # R  -> Sp
w2 = 1/(12*365.0)   # Rp -> S
w3 = 1/(7*365.0)    # Sp -> S
beta0 = R0 * gamma  # naive-track R0 = beta/gamma
b = 0.0116/365.0    # births/day (US CBR 11.6/1000/yr, NCHS)
d = b               # deaths/day = births -> N constant
SEAS = 0.08         # mild seasonal forcing amplitude

def rhs(t, y, p, seas=0.0):
    S,E,I,R,Sp,Ep,Ip,Rp = y
    beta = beta0*(1 + seas*np.cos(2*np.pi*t/365.0))
    lam = beta*(I + sigma*Ip)
    dS  = (1-p)*b - lam*S + w2*Rp + w3*Sp - d*S
    dE  = lam*S - eps*E - d*E
    dI  = eps*E - gamma*I - d*I
    dR  = gamma*I - w1*R - d*R
    dSp = p*b + w1*R - delta*lam*Sp - w3*Sp - d*Sp
    dEp = delta*lam*Sp - eps*Ep - d*Ep
    dIp = eps*Ep - gamma_p*Ip - d*Ip
    dRp = gamma_p*Ip - w2*Rp - d*Rp
    return [dS,dE,dI,dR,dSp,dEp,dIp,dRp]


def equilibrium(p):
    """Integrate unforced to the endemic steady state."""
    y0 = np.array([0.1,0.001,0.001,0.3,0.2,0.0005,0.0005,0.2]); y0 /= y0.sum()
    sol = solve_ivp(lambda t,y: rhs(t,y,p,0.0), (0, 600*365), y0,
                    method="LSODA", rtol=1e-9, atol=1e-12, t_eval=[600*365])
    return sol.y[:,-1]


def jacobian(yeq, p):
    f0 = np.array(rhs(0, yeq, p, 0.0)); n = len(yeq); J = np.zeros((n,n))
    for j in range(n):
        h = 1e-8*max(1e-8, abs(yeq[j]))
        yp = yeq.copy(); yp[j] += h
        J[:,j] = (np.array(rhs(0, yp, p, 0.0)) - f0)/h
    return J


def intrinsic_period(p):
    """Dominant damped-oscillation period (yr) from the Jacobian eigenvalues."""
    yeq = equilibrium(p)
    ev = np.linalg.eigvals(jacobian(yeq, p))
    osc = ev[np.abs(ev.imag) > 1e-9]            # complex (oscillatory) modes
    osc = osc[osc.imag > 0]
    dom = osc[np.argmax(osc.real)]              # slowest-decaying oscillation
    T_years = (2*np.pi/abs(dom.imag))/365.0
    damp_years = (-1.0/dom.real)/365.0 if dom.real < 0 else float("inf")
    return T_years, damp_years, yeq


def damped_trace(p):
    """Unforced transient from a perturbed IC, to show the damped cycle."""
    yeq = equilibrium(p); y0 = yeq.copy()
    y0[2] *= 3.0                                 # perturb I
    t_eval = np.arange(0, 60*365, 2.0)
    sol = solve_ivp(lambda t,y: rhs(t,y,p,0.0), (0, 60*365), y0,
                    method="LSODA", rtol=1e-8, atol=1e-11, t_eval=t_eval)
    return sol.t/365.0, eps*(sol.y[1]+sol.y[5])


res = {}
fig, axes = plt.subplots(2,1, figsize=(9,6), sharex=True)
for ax,(p,label,lit) in zip(axes, [(0.0,"Pre-vaccination","2-3 yr (pre-vax)"),
                                    (0.85,"With vaccination (85%)","3.3-4.1 yr (post-vax)")]):
    T_years, damp_years, yeq = intrinsic_period(p)
    tt, inc = damped_trace(p)
    res[label] = {"vaccinated_fraction": p, "intrinsic_period_years": round(float(T_years),2),
                  "damping_time_years": round(float(damp_years),1), "literature": lit}
    ax.plot(tt, inc, color="#0e7a78", lw=1.3)
    ax.set_title(f"{label}: intrinsic period ≈ {T_years:.1f} yr "
                 f"(damping ~{damp_years:.0f} yr) · literature {lit}", fontsize=10, weight="bold")
    ax.set_ylabel("Incidence"); ax.grid(alpha=.3)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    print(f"{label:28s} intrinsic period {T_years:.2f} yr | damping {damp_years:.1f} yr | lit {lit}")
axes[-1].set_xlabel("Year (after a perturbation)")
fig.suptitle("Demography added → intrinsic inter-epidemic cycle (Wearing-Rohani signature)",
             fontsize=12, weight="bold")
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(os.path.join(OUT,"val6_demography.png"), dpi=140)

# ---- sensitivity: intrinsic period vs natural-immunity duration (1/w2) ----
imm_years = [4, 8, 12, 20, 30, 40]
periods = []
_w2_default = w2
for iy in imm_years:
    globals()['w2'] = 1.0/(iy*365.0)
    Tp, _, _ = intrinsic_period(0.0)
    periods.append(round(float(Tp), 2))
globals()['w2'] = _w2_default
res["immunity_sweep"] = {"natural_immunity_years": imm_years, "intrinsic_period_years": periods}
print("\nperiod vs natural-immunity duration (pre-vax):")
for iy, pp in zip(imm_years, periods):
    print(f"  {iy:3d} yr immunity -> {pp:.2f} yr cycle")

fig2, ax = plt.subplots(figsize=(7.2, 4.0))
ax.axhspan(2.0, 4.6, color="#0e7a78", alpha=.10, label="empirical range 2–4.6 yr")
ax.plot(imm_years, periods, "o-", color="#0e7a78", lw=2)
ax.axvline(12, color="#c2410c", ls="--", lw=1.2, label="model default ω₂ = 12 yr")
ax.axvline(30, color="#475569", ls=":", lw=1.2, label="Wearing & Rohani ≥30 yr")
ax.set_xlabel("Natural-immunity duration 1/ω₂ (years)")
ax.set_ylabel("Intrinsic inter-epidemic period (years)")
ax.set_title("Cycle period lengthens with immunity duration", weight="bold", fontsize=11)
ax.legend(fontsize=8, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(alpha=.3)
fig2.tight_layout(); fig2.savefig(os.path.join(OUT, "val6b_immunity_sweep.png"), dpi=140)

res["birth_rate_per_1000"] = 11.6
res["source"] = "US CBR NCHS e6fc-ccez (2018); literature periods 2-4.6 yr (empirical), 3.3-4.1 yr post-vaccination"
json.dump(res, open(os.path.join(OUT,"val6_demography.json"),"w"), indent=2)
print("\nwrote val6_demography.json, val6_demography.png, val6b_immunity_sweep.png")
