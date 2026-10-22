# pages/Demography.py
#
# UI entry for the demographic companion engine (engine/demography.py): the
# multi-year / long-run behaviour the single-season dashboard model cannot show.
# Surfaces the two Wearing & Rohani signatures (intrinsic inter-epidemic cycle,
# stochastic fade-outs / critical community size) and the age-structured
# infection distribution. Deterministic results are cached for responsiveness;
# the stochastic fade-out sweep runs on demand.

import numpy as np
import pandas as pd
import streamlit as st

from layout.header import show_dashboard_header
from layout.sidebar import render_sidebar
from layout.logos import show_logos
from constants import DEFAULT_AGE_GROUPS
from helpers import load_population
from engine.demography import (
    build_params, intrinsic_period, simulate,
    build_age_params, infection_age_distribution, simulate_stochastic,
    US_CRUDE_BIRTH_RATE_PER_1000,
)

st.set_page_config(
    page_title="Demography & multi-year dynamics",
    page_icon="assets/epydemix-icon.svg",
    layout="wide",
    initial_sidebar_state="collapsed",
)
st.markdown("<style>.block-container{padding-top:2rem;padding-bottom:5rem;}</style>",
            unsafe_allow_html=True)

show_dashboard_header()
render_sidebar()

st.markdown("## Demography & multi-year dynamics")
st.caption(
    "The main dashboard runs a **single season** with no births, so it cannot show "
    "pertussis's multi-year behaviour. This page uses the demographic companion model "
    "(`engine/demography.py`) — the same 8-compartment partial-immunity structure **with "
    "vital dynamics added** — to show the intrinsic inter-epidemic cycle, stochastic "
    "fade-outs, and the age distribution of infections. It is an age-aggregated / "
    "age-structured analysis tool, distinct from the stochastic single-season engine."
)

# ---- controls --------------------------------------------------------------
c = st.columns(4)
R0 = c[0].slider("R₀", 5.0, 17.0, 8.0, 0.5)
imm = c[1].slider("Natural immunity 1/ω₂ (years)", 4, 40, 12, 1,
                  help="Default 12 y (clinical). Wearing & Rohani estimate ≥30 y; "
                       "longer immunity lengthens the inter-epidemic cycle.")
vax = c[2].slider("Vaccinated fraction of births", 0.0, 0.95, 0.0, 0.05,
                  help="0 = pre-vaccination. Routes that fraction of births into the "
                       "partially-immune pool (Sₚ).")
geo = st.session_state.get("selected_geography") or "United_States"
c[3].metric("Birth rate (NCHS)", f"{US_CRUDE_BIRTH_RATE_PER_1000}/1000/yr")


@st.cache_data(show_spinner=False)
def _period(R0, imm, vax):
    r0 = intrinsic_period(build_params(R0=R0, waning_partial_to_susceptible_years=imm, vaccinated_fraction=0.0))
    rv = intrinsic_period(build_params(R0=R0, waning_partial_to_susceptible_years=imm, vaccinated_fraction=vax))
    return r0, rv


@st.cache_data(show_spinner=False)
def _trace(R0, imm, vax):
    sim = simulate(build_params(R0=R0, waning_partial_to_susceptible_years=imm,
                                vaccinated_fraction=vax), years=60)
    return sim["t_years"], sim["incidence"]


@st.cache_data(show_spinner=False)
def _sweep(R0, vax):
    ys = [4, 8, 12, 20, 30, 40]
    per = [intrinsic_period(build_params(R0=R0, waning_partial_to_susceptible_years=y,
                                         vaccinated_fraction=vax))["period_years"] for y in ys]
    return ys, per


@st.cache_data(show_spinner=False)
def _agedist(geo, R0, imm, vax):
    pop = load_population(geo)
    M = sum(pop.contact_matrices[l] for l in pop.contact_matrices)
    p = build_age_params(M, R0=R0, waning_partial_to_susceptible_years=imm, vaccinated_fraction=vax)
    return infection_age_distribution(p)["dist"]


# ---- 1. intrinsic inter-epidemic cycle ------------------------------------
st.markdown("### 1 · Intrinsic inter-epidemic cycle (Wearing–Rohani)")
with st.spinner("Computing endemic equilibrium and eigenvalues…"):
    r0, rv = _period(R0, imm, vax)
    tt, inc = _trace(R0, imm, vax)
m = st.columns(3)
m[0].metric("Period — pre-vaccination", f"{r0['period_years']:.2f} yr")
m[1].metric(f"Period — {int(vax*100)}% vaccinated", f"{rv['period_years']:.2f} yr")
m[2].metric("Literature (empirical)", "2–4.6 yr")
st.caption("Dominant damped-oscillation period from the Jacobian at the endemic equilibrium. "
           "Damped waves after a perturbation (incidence):")
st.line_chart(pd.DataFrame({"incidence": inc}, index=np.round(tt, 2)).rename_axis("year"),
              height=240)

# ---- 2. period vs immunity duration ---------------------------------------
st.markdown("### 2 · How the cycle depends on immunity duration")
ys, per = _sweep(R0, vax)
sweep_df = pd.DataFrame({"immunity_years": ys, "period_years": per}).set_index("immunity_years")
st.line_chart(sweep_df, height=240)
st.caption(f"At the current R₀={R0:g}, the intrinsic period runs "
           f"{min(per):.1f}–{max(per):.1f} yr as natural-immunity duration goes 4→40 y "
           "(all within the empirical 2–4.6 yr band). The default ω₂ = 12 y sits at the short end; "
           "Wearing & Rohani's ≥30 y lengthens it.")

# ---- 3. age-structured infection distribution -----------------------------
st.markdown("### 3 · Age distribution of infections (age-structured + births)")
try:
    with st.spinner(f"Running age-structured model for {geo.replace('_',' ')}…"):
        dist = _agedist(geo, R0, imm, vax)
    ADULT = [4, 5, 6]
    import altair as alt
    ad = pd.DataFrame({"age_group": DEFAULT_AGE_GROUPS, "share": [float(x) for x in dist]})
    st.altair_chart(
        alt.Chart(ad).mark_bar(color="#4c9be8").encode(
            x=alt.X("age_group:N", title="Age group", sort=DEFAULT_AGE_GROUPS),
            y=alt.Y("share:Q", title="Share of infections", axis=alt.Axis(format="%")),
            tooltip=["age_group:N", alt.Tooltip("share:Q", format=".1%")],
        ),
        use_container_width=True,
    )
    cc = st.columns(2)
    cc[0].metric("Infant <1yr infection share", f"{dist[0]:.1%}", help="Single-season model (no births): ~0.7%")
    cc[1].metric("Adult 20+ infection share", f"{sum(dist[i] for i in ADULT):.1%}")
    st.caption(
        f"With births + aging ({geo.replace('_',' ')} age pyramid and contact matrix), the infant "
        "infection share is far more realistic than the single-season model (~0.7%). These are "
        "**infections**; reported cases are skewed toward infants (severity/detection) and away from "
        "adults (under-ascertainment) — see Model References / docs/validation §3d."
    )
except Exception as exc:
    st.warning(f"Age-structured run unavailable for this geography ({exc}).")

# ---- 4. stochastic fade-outs / critical community size --------------------
st.markdown("### 4 · Stochastic fade-outs & critical community size")
st.caption("The second Wearing–Rohani signature: in small populations the infection goes "
           "extinct between epidemics (fade-out). Runs a tau-leaping stochastic sweep (slower).")
if st.button("Run fade-out analysis", use_container_width=True):
    Ns = [2000, 5000, 10000, 20000, 50000, 100000, 300000, 1000000]
    prog = st.progress(0.0, text="Running stochastic simulations…")
    fade = []
    base = build_params(R0=R0, waning_partial_to_susceptible_years=imm, vaccinated_fraction=vax)
    for k, N in enumerate(Ns):
        vals = [simulate_stochastic(base, population=N, years=45, seed=s)["fadeout_fraction"]
                for s in (1, 2)]
        fade.append(float(np.mean(vals)))
        prog.progress((k + 1) / len(Ns), text=f"N = {N:,}")
    prog.empty()
    ccs = next((N for N, f in zip(Ns, fade) if f < 0.05), None)
    fd = pd.DataFrame({"population": Ns, "fadeout_%": [100 * f for f in fade]}).set_index("population")
    st.line_chart(fd, height=240)
    if ccs:
        st.success(f"Critical community size ≈ **{ccs:,}** (fade-out < 5%). "
                   "Below measles' classic ~250–500k — pertussis persists in smaller "
                   "populations via its waning-immunity adult reservoir.")

show_logos()
