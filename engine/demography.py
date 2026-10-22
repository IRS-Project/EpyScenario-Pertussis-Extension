# engine/demography.py
#
# Demographic (vital-dynamics) companion to the single-season epydemix pertussis
# stub. epydemix's transition framework is compartment-to-compartment only, so it
# cannot express births (a source with no origin compartment) or aging (movement
# between demographic groups). This module therefore provides a self-contained,
# age-aggregated deterministic implementation of the SAME 8-compartment partial-
# immunity structure WITH vital dynamics added:
#
#     births -> S (naive)           deaths <- every compartment (rate d)
#
# It is the engine path for multi-year / inter-epidemic analysis — e.g. the
# Wearing & Rohani "epidemiological signatures" (intrinsic cycle period, fade-outs)
# — which the closed-population single-season model cannot produce. Parameters
# default to the app's SEIRS (Pertussis) schema so the two stay in sync.
#
# Scope: age-AGGREGATED (homogeneous mixing). The inter-epidemic cycle is governed
# by susceptible turnover, which this captures; age-structured demography (births
# into the youngest band + aging across the 7 bands) is a further extension.

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp

from schemas import MODEL_PARAM_SCHEMAS

_PERT = "SEIRS (Pertussis)"
COMPARTMENTS = ["S", "E", "I", "R", "Sp", "Ep", "Ip", "Rp"]

# US crude birth rate, NCHS "Births and General Fertility Rates" (e6fc-ccez), 2018.
US_CRUDE_BIRTH_RATE_PER_1000 = 11.6


def _schema_defaults() -> dict:
    return {p["key"]: p["default"] for p in MODEL_PARAM_SCHEMAS[_PERT]}


def build_params(
    R0: float | None = None,
    incubation_period: float | None = None,
    infectious_period: float | None = None,
    infectious_period_partial: float | None = None,
    rel_infectiousness_partial: float | None = None,
    rel_susceptibility_partial: float | None = None,
    waning_full_to_partial_years: float | None = None,
    waning_partial_to_susceptible_years: float | None = None,
    waning_vaccine_to_susceptible_years: float | None = None,
    birth_rate_per_1000: float = US_CRUDE_BIRTH_RATE_PER_1000,
    vaccinated_fraction: float = 0.0,
    seasonality_amp: float = 0.0,
) -> dict:
    """Assemble the rate dictionary (per DAY) for the demographic model.

    Any epidemiological argument left None falls back to the app's SEIRS
    (Pertussis) schema default. `birth_rate_per_1000` is an annual crude rate;
    deaths are set equal to births so total population is conserved.
    `vaccinated_fraction` p routes that fraction of births into Sp (partially
    immune) instead of S — p=0 is the pre-vaccination regime.
    """
    d = _schema_defaults()
    g = lambda v, k: float(d[k] if v is None else v)
    R0 = g(R0, "R0")
    inc = g(incubation_period, "incubation_period")
    infd = g(infectious_period, "infectious_period")
    infp = g(infectious_period_partial, "infectious_period_partial")
    sigma = g(rel_infectiousness_partial, "rel_infectiousness_partial")
    delta = g(rel_susceptibility_partial, "rel_susceptibility_partial")
    w1y = g(waning_full_to_partial_years, "waning_full_to_partial_years")
    w2y = g(waning_partial_to_susceptible_years, "waning_partial_to_susceptible_years")
    w3y = g(waning_vaccine_to_susceptible_years, "waning_vaccine_to_susceptible_years")

    gamma = 1.0 / infd          # I  -> R
    b = (birth_rate_per_1000 / 1000.0) / 365.0
    return {
        "beta0": R0 * gamma,    # naive-track R0 = beta / gamma
        "eps": 1.0 / inc,       # E -> I, Ep -> Ip
        "gamma": gamma,         # I  -> R
        "gamma_p": 1.0 / infp,  # Ip -> Rp
        "sigma": sigma, "delta": delta,
        "w1": 1.0 / (w1y * 365.0),  # R  -> Sp
        "w2": 1.0 / (w2y * 365.0),  # Rp -> S
        "w3": 1.0 / (w3y * 365.0),  # Sp -> S
        "b": b, "d": b,             # births = deaths -> N constant
        "p": float(vaccinated_fraction),
        "seas": float(seasonality_amp),
    }


def rhs(t: float, y, p: dict):
    """8-compartment partial-immunity SEIRS with vital dynamics (N normalised to 1)."""
    S, E, I, R, Sp, Ep, Ip, Rp = y
    beta = p["beta0"] * (1.0 + p["seas"] * np.cos(2 * np.pi * t / 365.0))
    lam = beta * (I + p["sigma"] * Ip)
    b, d, pp = p["b"], p["d"], p["p"]
    dS = (1 - pp) * b - lam * S + p["w2"] * Rp + p["w3"] * Sp - d * S
    dE = lam * S - p["eps"] * E - d * E
    dI = p["eps"] * E - p["gamma"] * I - d * I
    dR = p["gamma"] * I - p["w1"] * R - d * R
    dSp = pp * b + p["w1"] * R - p["delta"] * lam * Sp - p["w3"] * Sp - d * Sp
    dEp = p["delta"] * lam * Sp - p["eps"] * Ep - d * Ep
    dIp = p["eps"] * Ep - p["gamma_p"] * Ip - d * Ip
    dRp = p["gamma_p"] * Ip - p["w2"] * Rp - d * Rp
    return [dS, dE, dI, dR, dSp, dEp, dIp, dRp]


def equilibrium(p: dict, settle_years: int = 600) -> np.ndarray:
    """Integrate (unforced) to the endemic steady state."""
    q = dict(p); q["seas"] = 0.0
    y0 = np.array([0.1, 1e-3, 1e-3, 0.3, 0.2, 5e-4, 5e-4, 0.2]); y0 /= y0.sum()
    sol = solve_ivp(lambda t, y: rhs(t, y, q), (0, settle_years * 365), y0,
                    method="LSODA", rtol=1e-9, atol=1e-12, t_eval=[settle_years * 365])
    return sol.y[:, -1]


def _jacobian(yeq: np.ndarray, p: dict) -> np.ndarray:
    q = dict(p); q["seas"] = 0.0
    f0 = np.array(rhs(0.0, yeq, q)); n = len(yeq); J = np.zeros((n, n))
    for j in range(n):
        h = 1e-8 * max(1e-8, abs(yeq[j]))
        yp = yeq.copy(); yp[j] += h
        J[:, j] = (np.array(rhs(0.0, yp, q)) - f0) / h
    return J


def intrinsic_period(p: dict) -> dict:
    """Wearing-&-Rohani 'natural period': the dominant damped oscillation at the
    endemic equilibrium, from the Jacobian's leading complex eigenvalue.

    Returns {'period_years', 'damping_years', 'equilibrium'}.
    """
    yeq = equilibrium(p)
    ev = np.linalg.eigvals(_jacobian(yeq, p))
    osc = ev[(np.abs(ev.imag) > 1e-9) & (ev.imag > 0)]
    if len(osc) == 0:
        return {"period_years": float("nan"), "damping_years": float("inf"), "equilibrium": yeq}
    dom = osc[np.argmax(osc.real)]
    return {
        "period_years": (2 * np.pi / abs(dom.imag)) / 365.0,
        "damping_years": (-1.0 / dom.real) / 365.0 if dom.real < 0 else float("inf"),
        "equilibrium": yeq,
    }


# ---------------------------------------------------------------------------
# Age-structured variant (§6 step 1): 7 bands x 8 compartments, births into the
# youngest band, aging through bands, deaths from the oldest (type-I survival),
# age-mixing via a contact matrix. Continuous infant births make the infant
# INFECTION share realistic (the single-season model cannot, having no influx).
# ---------------------------------------------------------------------------
AGE_WIDTHS_YEARS = [1.0, 5.0, 4.0, 9.0, 30.0, 15.0, 13.0]  # last = life exp − 65


def build_age_params(contact_matrix, R0: float | None = None,
                     width_years=AGE_WIDTHS_YEARS, vaccinated_fraction: float = 0.0,
                     **overrides) -> dict:
    """Age-structured parameters. `contact_matrix` is the combined (summed-layer)
    7×7 matrix; β is calibrated to R0 via its spectral radius, as the app does."""
    base = build_params(R0=R0, vaccinated_fraction=vaccinated_fraction, **overrides)
    M = np.asarray(contact_matrix, float)
    nb = M.shape[0]
    w = np.asarray(width_years, float)
    a = np.array([1.0/(w[i]*365.0) for i in range(nb)]); a[-1] = 0.0
    dvec = np.zeros(nb); dvec[-1] = 1.0/(w[-1]*365.0)
    B = 1.0/(w.sum()*365.0)
    base.update({"M": M, "nb": nb, "a": a, "dvec": dvec, "B": B, "Nstat": w/w.sum()})
    # β from R0 and the contact-matrix spectral radius (matches engine.run.compute_beta)
    R0_eff = R0 if R0 is not None else _schema_defaults()["R0"]
    base["beta0"] = R0_eff * base["gamma"] / np.linalg.eigvals(M).real.max()
    return base


def age_rhs(t: float, y, p: dict):
    nb = p["nb"]; Y = y.reshape(8, nb)
    S, E, I, R, Sp, Ep, Ip, Rp = Y
    Nband = np.where(Y.sum(axis=0) > 1e-12, Y.sum(axis=0), 1e-12)
    lam = p["beta0"] * (p["M"] @ ((I + p["sigma"]*Ip)/Nband))
    a, dvec, B = p["a"], p["dvec"], p["B"]
    def age(X):
        inn = np.zeros(nb); inn[1:] = a[:-1]*X[:-1]
        return inn - a*X
    births = np.zeros(nb); births[0] = B
    dS  = births - lam*S + p["w2"]*Rp + p["w3"]*Sp + age(S) - dvec*S
    dE  = lam*S - p["eps"]*E + age(E) - dvec*E
    dI  = p["eps"]*E - p["gamma"]*I + age(I) - dvec*I
    dR  = p["gamma"]*I - p["w1"]*R + age(R) - dvec*R
    dSp = p["w1"]*R - p["delta"]*lam*Sp - p["w3"]*Sp + age(Sp) - dvec*Sp
    dEp = p["delta"]*lam*Sp - p["eps"]*Ep + age(Ep) - dvec*Ep
    dIp = p["eps"]*Ep - p["gamma_p"]*Ip + age(Ip) - dvec*Ip
    dRp = p["gamma_p"]*Ip - p["w2"]*Rp + age(Rp) - dvec*Rp
    return np.concatenate([dS, dE, dI, dR, dSp, dEp, dIp, dRp])


def age_equilibrium(p: dict, settle_years: int = 400) -> np.ndarray:
    nb = p["nb"]; Nstat = p["Nstat"]
    Y0 = np.zeros((8, nb)); Y0[0] = Nstat*0.5; Y0[3] = Nstat*0.5
    Y0[2] = Nstat*1e-4; Y0[0] -= Y0[2]
    sol = solve_ivp(lambda t, y: age_rhs(t, y, p), (0, settle_years*365), Y0.flatten(),
                    method="LSODA", rtol=1e-8, atol=1e-11, t_eval=[settle_years*365])
    return sol.y[:, -1].reshape(8, nb)


def infection_age_distribution(p: dict) -> dict:
    """Endemic incidence (ε·(E+Ep)) by age band, normalised."""
    Yeq = age_equilibrium(p)
    E, Ep = Yeq[1], Yeq[5]
    inc = p["eps"]*(E + Ep); dist = inc/inc.sum()
    adult = [i for i in range(p["nb"]) if i >= 4]
    return {"dist": dist, "infant": float(dist[0]), "adult": float(sum(dist[adult])),
            "pop_fraction": Yeq.sum(axis=0)}


def simulate(p: dict, years: int = 120, perturb: bool = True, dt_days: float = 2.0) -> dict:
    """Integrate the forced/unforced model and return trajectories + incidence.

    Returns {'t_years', 'incidence', 'compartments' (dict band->array)}.
    """
    y0 = equilibrium(p).copy()
    if perturb:
        y0[2] *= 3.0  # kick I to excite the damped oscillation
    t_eval = np.arange(0, years * 365, dt_days)
    sol = solve_ivp(lambda t, y: rhs(t, y, p), (0, years * 365), y0,
                    method="LSODA", rtol=1e-8, atol=1e-11, t_eval=t_eval)
    comp = {c: sol.y[i] for i, c in enumerate(COMPARTMENTS)}
    incidence = p["eps"] * (comp["E"] + comp["Ep"])
    return {"t_years": sol.t / 365.0, "incidence": incidence, "compartments": comp}


# ---------------------------------------------------------------------------
# Stochastic version (§6 step 2): tau-leaping on the age-aggregated model with
# vital dynamics + rare external imports. Used to measure FADE-OUT frequency
# (fraction of weeks with zero infectives) vs population size — the second
# Wearing & Rohani signature / critical community size.
# ---------------------------------------------------------------------------
def simulate_stochastic(p: dict, population: int = 500_000, years: int = 60,
                        dt: float = 1.0, import_rate_per_year: float = 3.0,
                        burn_in_years: int = 10, seas: float = 0.12, seed: int = 0) -> dict:
    """Tau-leaping SEIRS+demography with seasonal forcing + rare imports. Returns
    weekly infective totals and the fade-out fraction (weeks whose TROUGH hits zero
    I+Ip after burn-in) — the fade-out / critical-community-size signature."""
    rng = np.random.default_rng(seed)
    yeq = equilibrium(p)
    x = np.maximum(np.round(yeq * population).astype(np.int64), 0)  # S E I R Sp Ep Ip Rp
    imp = import_rate_per_year / 365.0
    nsteps = int(years * 365 / dt)
    weekly_min, wk = [], 0.0
    acc = None
    for step in range(nsteps):
        N = max(int(x.sum()), 1)
        S, E, I, R, Sp, Ep, Ip, Rp = x
        beta_t = p["beta0"] * (1.0 + seas * np.cos(2 * np.pi * (step * dt) / 365.0))
        lam = beta_t * (I + p["sigma"] * Ip) / N
        def po(rate):
            return rng.poisson(max(rate, 0.0) * dt)
        # events (cap outflows at the compartment)
        se = min(S, po(lam * S) + po(imp))            # infection + imports -> E
        spep = min(Sp, po(p["delta"] * lam * Sp))
        ei = min(E, po(p["eps"] * E)); epip = min(Ep, po(p["eps"] * Ep))
        ir = min(I, po(p["gamma"] * I)); iprp = min(Ip, po(p["gamma_p"] * Ip))
        rsp = min(R, po(p["w1"] * R)); rps = min(Rp, po(p["w2"] * Rp)); sps = min(Sp - spep, po(p["w3"] * Sp))
        births = po(p["b"] * N)
        # deaths (independent per compartment, capped)
        dd = np.array([min(x[k], po(p["d"] * x[k])) for k in range(8)], np.int64)
        nx = x.copy()
        nx[0] += births - se + rps + sps - dd[0]
        nx[1] += se - ei - dd[1]
        nx[2] += ei - ir - dd[2]
        nx[3] += ir - rsp - dd[3]
        nx[4] += rsp - spep - sps - dd[4]
        nx[5] += spep - epip - dd[5]
        nx[6] += epip - iprp - dd[6]
        nx[7] += iprp - rps - dd[7]
        x = np.maximum(nx, 0)
        cur = int(x[2] + x[6])
        acc = cur if acc is None else min(acc, cur)   # weekly TROUGH
        wk += dt
        if wk >= 7.0:
            weekly_min.append(acc); acc = None; wk = 0.0
    weekly = np.array(weekly_min)
    burn_wk = int(burn_in_years * 52)
    post = weekly[burn_wk:] if len(weekly) > burn_wk else weekly
    fadeout = float(np.mean(post == 0)) if len(post) else float("nan")
    return {"weekly_trough_infectives": weekly, "fadeout_fraction": fadeout, "population": population}
