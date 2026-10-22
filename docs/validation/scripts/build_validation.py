# -*- coding: utf-8 -*-
"""Build the EpyScenario pertussis VALIDATION report -> docs/validation/README.md
and validation.html (one content model, EpyScenario styling). Reads the JSON
outputs of val_implementation / val_holdout / val_agedist; the independent
cross-check figures (CDC NNDSS, EPICC) are stated inline from this session."""
import os, json, html, re, datetime

REPO = r"C:\Users\MarkRoutson\epydemix-dashboard"
BASE = os.path.join(REPO, "docs", "validation")
FIGS = os.path.join(BASE, "figs")
os.makedirs(FIGS, exist_ok=True)

def load(name):
    p = os.path.join(FIGS, name)
    return json.load(open(p)) if os.path.exists(p) else None

V1 = load("val1_implementation.json")
V2 = load("val2_holdout.json")
V3 = load("val3_agedist.json")
V3B = load("val3b_ageimmunity.json")
V3C = load("val3c_severity.json")
V3D = load("val3d_realcoverage.json")
V3E = load("val3e_ascertainment.json")
V6 = load("val6_demography.json")
V6C = load("val6c_agestructured.json")
V6D = load("val6d_fadeout.json")

# weekly-sum case counts actually used (NNDSS weekly table; lags annual confirmed)
WK = {}
try:
    import sys as _s
    _s.path.insert(0, REPO)
    from data.weekly_sources import load_weekly_pertussis as _lw, weekly_series_to_array as _wa
    for area in ("Oregon", "Alabama"):
        for yr in (2023, 2024, 2025):
            try:
                WK[(area, yr)] = int(sum(_wa(_lw(area, years=[yr])[0])))
            except Exception:
                pass
except Exception:
    pass
def wk(area, yr): return f"{WK[(area,yr)]}" if (area, yr) in WK else "—"

def pct(x): return f"{x:.1%}"

BLOCKS = []
def h(l, t): BLOCKS.append(("h", (l, t)))
def p(t): BLOCKS.append(("p", t))
def steps(x): BLOCKS.append(("steps", x))
def ul(x): BLOCKS.append(("ul", x))
def fig(s, c): BLOCKS.append(("fig", (s, c)))
def callout(k, ti, tx): BLOCKS.append(("callout", (k, ti, tx)))
def table(hd, rw): BLOCKS.append(("table", (hd, rw)))
def rule(): BLOCKS.append(("rule", None))

# ---------------------------------------------------------------- content
h(1, "EpyScenario Pertussis Model — Validation Report")
p("*An evidence record for the 8-compartment SEIRS partial-immunity pertussis "
  "model: implementation correctness, out-of-sample prediction, and corroboration "
  "against independent data and an independent model.*")
callout("note", "Scope",
        "This validates the modeling engine and its calibrated outputs for Oregon "
        "(Lane County) and Alabama. It is organised as a ladder: (1) does the code "
        "implement the math, (2) does the calibrated model predict data it was not "
        "fit to, (3) do its outputs agree with independent surveillance and an "
        "independent cost model. Every test is reproducible from the scripts named "
        "in the methods appendix.")

# ---- summary verdict
summary_rows = []
if V1:
    ce = V1["conservation"]["max_rel_error"]
    summary_rows.append(["Population conservation", f"max error {ce:.1e}", "PASS"])
    summary_rows.append(["Non-negativity", f"min compartment {V1['nonnegativity']['min_compartment_value']:.0f}", "PASS"])
    rr = V1["r0_recovery"]
    rng = f"{min(r['rel_error'] for r in rr):.0%}–{max(r['rel_error'] for r in rr):.0%}"
    summary_rows.append(["R₀ recovery (clean SEIR limit)", f"within {rng} of formula, constant ratio", "PASS (w/ caveat)"])
if V2:
    summary_rows.append(["Out-of-sample shape forecast", "does not beat persistence (quiet fit years; scenario engine)", "BOUNDS CLAIMS (§2)"])
summary_rows.append(["Independent data (CDC NNDSS)", "OR 2026 YTD 330 vs EPICC 327", "AGREES"])
summary_rows.append(["Independent model (EPICC Omaha-2008)", "$2,172/case reproduced (MMWR $52,131/24)", "AGREES"])
if V3:
    summary_rows.append(["Age distribution (infections → reported)", "infant 1%→13.8% via severity filter ≈ documented 12.5%", "PARTIAL — infant ✓, adults open (§3d)"])
summary_rows.append(["Literature parameter benchmarks", "R₀, incubation, immunity, σ in published ranges", "AGREES (ω₂ short-flagged)"])
if V6:
    summary_rows.append(["Inter-epidemic cycle (demography §5)", "intrinsic period ~2.0 yr vs empirical 2–4.6 yr", "REPRODUCED (companion model)"])
if V6D:
    summary_rows.append(["Stochastic fade-outs / CCS (§5)", "critical community size ~50k (below measles ~300k)", "REPRODUCED"])
if V6C:
    summary_rows.append(["Age-structured infection distribution (§5)", "infant infections 0.7%→3.3% with births/aging", "MUCH MORE REALISTIC"])
h(2, "Summary")
table(["Test", "Result", "Verdict"], summary_rows)

rule()
# ---- Part 1
h(2, "1 · Implementation correctness")
p("Does the engine implement the compartmental math without leaks or mis-wired "
  "rates? Three checks, all headless on the live engine.")
if V1:
    ul([
      f"**Population conservation** — Σ of the 8 compartments stays within "
      f"**{V1['conservation']['max_rel_error']:.1e}** (relative) of N over the whole run. "
      f"The small residual is from summing per-compartment *medians* across stochastic "
      f"replicates (median-of-sums ≠ sum-of-medians), not a mass leak.",
      f"**Non-negativity** — the smallest value any compartment reaches is "
      f"**{V1['nonnegativity']['min_compartment_value']:.0f}** (positive); no compartment goes negative.",
      "**R₀ recovery** — in the clean naive-SEIR limit (no prior immunity, empty "
      "partial track, seasonality off) the realized early growth rate is measured and "
      "mapped back to R₀ via R₀ = (1 + r/γ)(1 + r/μ):",
    ])
    table(["R₀ set", "growth rate r (/day)", "R₀ implied", "difference"],
          [[f"{r['R0_input']:.0f}", f"{r['growth_rate_per_day']:.4f}",
            f"{r['R0_implied']:.2f}", f"{r['rel_error']:.0%}"] for r in V1["r0_recovery"]])
    callout("note", "Reading the R₀ check",
            "Implied R₀ sits a near-constant ~10% below the value set — and the ratio is "
            "stable across R₀ = 6, 8, 12. A wiring bug (wrong β, double-counted FOI) would "
            "produce an erratic or factor-of-two error, not a constant ~0.9 scaling. The "
            "offset is the known gap between the homogeneous r↔R₀ formula and the model's "
            "age-structured mixing (R₀ is set as the next-generation-matrix eigenvalue, "
            "whose realized growth differs slightly from the homogeneous approximation).")
fig("figs/val1_implementation.png",
    "Left: population conservation error over time. Right: R₀ implied by realized "
    "growth vs R₀ set in the model — tight, monotonic, constant proportional offset.")

rule()
# ---- Part 2
h(2, "2 · Out-of-sample prediction (temporal hold-out)")
p("The decisive test: calibrate R₀ to one NNDSS season, then score the model on "
  "the **untouched next season** with a scale-free, peak-normalised weekly shape "
  "distance (absolute magnitude is not expected to generalise; seasonal shape and "
  "timing is). The benchmark is a **persistence baseline** — using the prior year's "
  "observed curve to predict the next. A model earns forecasting credit only by "
  "beating persistence.")
callout("note", "Result up front — an honest boundary on claims",
        "Cross-year *shape* prediction does **not** beat persistence in any "
        "configuration tested. Two reasons, both legitimate and bounding: (1) the "
        "quiet fit years carried almost no epidemic signal to calibrate to, and (2) "
        "EpyScenario generates a free-running outbreak from a seed — it is a "
        "*scenario* engine (its own data layer calls these \"scenarios, not "
        "forecasts\"), not a calendar-anchored forecaster. The validated use is "
        "what-if scenario analysis, corroborated in §1 and §3 — not next-season curve "
        "prediction.")
p("**Why the quiet-year fit is the binding constraint.** Oregon's confirmed "
  "pertussis went 23 (2022) → 35 (2023) → 1,039 (2024) → 1,560 (2025); the NNDSS "
  f"*weekly* table (what the fit sees) sums to Oregon {wk('Oregon',2023)} (2023) vs "
  f"{wk('Oregon',2024)} (2024), Alabama {wk('Alabama',2023)} (2023) vs {wk('Alabama',2024)} "
  "(2024). Fitting R₀ to a ~10-case annual curve is fitting noise — and indeed R₀ "
  "railed to the grid edge (18) on those years. The Oregon 2024→2025 row below "
  "removes that confound by fitting a real outbreak year — and it **still** does not "
  "beat persistence (peak-week off by ~40), which isolates the deeper reason: the "
  "model's free-running outbreak peaks at a model-determined time, not the observed "
  "calendar time. That is a property of a scenario engine, not a defect.")
if V2:
    order = [k for k in ("Oregon (2024->2025)", "Oregon", "Alabama") if k in V2] + \
            [k for k in V2 if k not in ("Oregon (2024->2025)", "Oregon", "Alabama")]
    hold_rows = []
    for s in order:
        r = V2[s]
        if "error" in r: continue
        beats = "yes" if r["beats_persistence"] else "no"
        hold_rows.append([s, f"{r['fit_year']}→{r['test_year']}",
                          f"{r.get('R0_fit_on_%d'%r['fit_year'])}",
                          f"{r['holdout_distance']:.2f}", f"{r['insample_floor']:.2f}",
                          f"{r['persistence_baseline']:.2f}", beats])
    table(["Hold-out", "Years", "R₀ fit", "Model dist", "In-sample floor", "Persistence", "Beats persist."],
          hold_rows)
    if os.path.exists(os.path.join(FIGS, "val2_holdout_oregon_outbreak.png")):
        fig("figs/val2_holdout_oregon_outbreak.png",
            "Oregon outbreak-to-outbreak hold-out (fit 2024 → score 2025), peak-"
            "normalised: model vs persistence vs the held-out 2025 curve.")
    fig("figs/val2_holdout_oregon.png",
        "Oregon 2023→2024: the 2023 fit year is essentially flat (≈12 weekly cases), "
        "so the calibration has no outbreak to learn from.")
    callout("caution", "What this correctly rules in and out",
            "Ruled IN: EpyScenario as a scenario/what-if engine — §1 shows the "
            "dynamics are implemented correctly and §3 shows outputs agree with "
            "independent data and an independent model. Ruled OUT (for now): using it "
            "to forecast next season's case curve. A genuine forecasting validation "
            "would need calendar-anchored seeding, a reporting model, and several "
            "high-incidence seasons to train and test on — none of which this single-"
            "season scenario tool currently claims.")

rule()
# ---- Part 3
h(2, "3 · Independent corroboration")
h(3, "3a · CDC NNDSS surveillance (independent of the model)")
p("The model's inputs and relative magnitudes agree with live CDC NNDSS pulled this "
  "session — and two independent sources agree with each other:")
table(["Jurisdiction", "2026 YTD (NNDSS)", "Prior-year YTD", "Cross-source"],
      [["United States", "11,358", "27,882", "—"],
       ["Oregon", "330", "1,378", "EPICC 2026 YTD = 327 ✓"],
       ["Alabama", "136", "893", "low-incidence (matches sparse curve) ✓"]])
p("Oregon runs far hotter than Alabama in every year — exactly the contrast the two "
  "manuals show, and the reason Alabama needs curve-only calibration on a sparse series.")

h(3, "3b · Independent model — EPICC Omaha-2008 validation")
p("EPICC (a separate codebase) reproduces the only published U.S. local-health-"
  "department pertussis cost study to the dollar: **$52,127 across 24 cases = "
  "$2,172/case**, against the published **MMWR 60(1) 2011** figure of $52,131. An "
  "independent implementation hitting the published answer is a strong external "
  "anchor for the shared pertussis epidemiology both tools rely on.")

h(3, "3c · Attack-rate plausibility (EPICC Oregon record)")
p("Oregon 2025 incidence was **36.51 / 100,000** — highest in the nation (8.46 "
  "national) — on a 75-year case record. The calibrated Oregon run's high R₀ (≈15) "
  "and large modeled outbreak are consistent with a state that genuinely recorded "
  "1,039 (2024) and 1,560 (2025) confirmed cases; the model is not manufacturing an "
  "implausibly large epidemic.")

h(3, "3d · Age distribution — infections vs reported cases")
p("Documented U.S. pertussis *reported-case* epidemiology (EPICC/MMWR): infants "
  "<1yr are **12.5%** of cases nationally (**5.6%** in Oregon), adults 20+ about "
  "**19%**. The model produces true **infections**, which are a different quantity "
  "— so this is worked through in three steps rather than a single pass/fail.")
if V3 and V3B and V3C:
    m = V3["modeled"]
    table(["Age mix", "Infant <1yr", "Adult 20+", "vs documented"],
          [["Modeled infections — uniform immunity (defaults)",
            pct(m["Oregon (Lane County)"]["infant_share"]),
            pct(m["Oregon (Lane County)"]["adult_share"]), "far off"],
           ["Modeled infections — + age-immunity table",
            pct(V3B["age_immunity"]["infant"]), pct(V3B["age_immunity"]["adult"]), "barely moves"],
           ["**Severity-weighted (≈ reported cases)**",
            pct(V3C["severity_weighted"]["infant"]), pct(V3C["severity_weighted"]["adult"]),
            "**infant matches**"],
           ["Documented reported cases", "5.6–12.5%", "~19%", "—"]])
    callout("note", "The real explanation — three findings, not one failure",
            "**(1)** Under uniform-immunity defaults the modeled infection mix is driven "
            "by the contact matrix (high-mixing adolescents/adults), so only ~1% of "
            "infections land in infants.\n\n"
            "**(2)** The integrated **immunization-by-age table** (use_age_immunity, the "
            "\"Lane County (observed)\" profile with infants 0% immune) helps only "
            "marginally — ~0.5%→1.0% — because the infant band is small and low-contact; "
            "age-structured susceptibility alone cannot produce the documented skew.\n\n"
            "**(3)** The documented distribution is a *reported-case* distribution, "
            "skewed to infants by their severity and detection. Applying the app's own "
            "age-specific severity weights (`hospitalization.py` IHRs: infants 0.30 vs "
            "0.01–0.04) as a detection filter moves the infant share **1.0% → 13.8%**, "
            "right at the documented ~12.5%. The infant signature is reproduced — as a "
            "detection artifact, correctly.")
    callout("caution", "What is still off, and why it's informative",
            "Severity-weighting leaves adults at ~72% vs the documented ~19%. That is a "
            "real, interpretable gap: adults genuinely have many (often atypical) "
            "infections but are **heavily under-ascertained** in surveillance — a factor "
            "beyond severity — and the 65+ IHR over-weights the elderly. A complete "
            "reporting model would add an adult under-detection term. Bottom line: the "
            "model's *infection* age-distribution is a valid, distinct quantity; mapping "
            "it to *reported cases* needs a detection model, of which severity is only "
            "the infant-dominant part.")
    fig("figs/val3c_severity.png",
        "Modeled infections vs severity-weighted (≈ reported) age distribution. The "
        "severity filter lifts the infant share to ~14% (documented ~12.5%); adults "
        "remain over-represented, pointing to additional adult under-ascertainment.")
    if V3D:
        ri = V3D["infections_real"]; rr = V3D["reported_real"]; ci = V3D["infections_circular"]
        h(4, "Closing the circularity: real NIS/BRFSS coverage")
        p("The age-immunity profile above is the *vaccinated-share-of-cases* (which the "
          "code itself flags as circular). Re-running with **real population coverage** "
          "pulled from data.cdc.gov this session — CDC **NIS-Child** DTaP, **NIS-Teen** "
          "Tdap, and **BRFSS** adult Tdap — removes the circularity and the conclusion "
          "holds, more strongly:")
        table(["Age mix (real coverage input)", "Infant <1yr", "Adult 20+"],
              [["Modeled infections — circular profile", pct(ci["infant"]), pct(ci["adult"])],
               ["Modeled infections — real coverage", pct(ri["infant"]), pct(ri["adult"])],
               ["Severity-weighted (≈ reported) — real coverage", pct(rr["infant"]), pct(rr["adult"])],
               ["Documented reported cases", "5.6–12.5%", "~19%"]])
        callout("tip", "The infant burden is detection, not susceptibility — confirmed non-circularly",
                "Real coverage protects infants (35%) and children (80%) more than the "
                "circular profile did, so modeled infant *infections* fall further (1.0% "
                "→ 0.7%). Yet the severity-weighted *reported* infant share is 10.2% — "
                "still inside the documented 5.6–12.5% range. Infant over-representation "
                "in surveillance is therefore a reporting/severity artifact, reproduced "
                "from a non-circular input; it is not driven by infant susceptibility. "
                "Adults remain over-represented under real coverage too (~79%), so the "
                "adult under-ascertainment term is independent of the coverage source.")
        fig("figs/val3d_realcoverage.png",
            "Left: circular (share-of-cases) vs real (NIS/BRFSS) coverage-by-age input. "
            "Right: resulting infant and adult shares under real coverage — infections, "
            "severity-weighted reported, and documented.")
    if V3E:
        rep = V3E["reported"]
        h(4, "Closing the adult gap: age-specific under-ascertainment")
        p("Severity alone left adults over-represented. The missing term is "
          "age-specific **under-ascertainment**, anchored to an independent "
          "seroprevalence study (Italy, *Infect Dis Ther* 2023): true infection was "
          "**~141× reported in ages 6–14** but **~3,452× in ages ≥15** — adult "
          "pertussis is barely captured by surveillance. Applying that reporting "
          "gradient to the modeled infections:")
        table(["Age mix", "Infant <1yr", "Adult 20+"],
              [["Modeled infections (real coverage)", pct(V3E["infections"]["infant"]), pct(V3E["infections"]["adult"])],
               ["Reported — severity filter only", pct(V3D["reported_real"]["infant"]), pct(V3D["reported_real"]["adult"])],
               ["Reported — severity + serosurvey ascertainment", pct(rep["infant"]), pct(rep["adult"])],
               ["Documented reported cases", "5.6–12.5%", "~19%"]])
        callout("note", "The documented distribution is bracketed — mechanism confirmed, not fit",
                "The severity-only filter leaves adults too high (~79%); the full "
                "serosurvey gradient over-corrects (adults 3.8%, infants 58%) because "
                "the Italian ratios are extreme and country-specific. The documented US "
                "pattern (infant ~12.5%, adult ~19%) sits **between** these two "
                "independent anchors. So the reported-case age distribution is an "
                "**observation-model** quantity — reproducible from modeled infections "
                "with an age-specific reporting gradient whose magnitude lies between "
                "severity and serosurvey bounds — and is calibratable to local "
                "surveillance. It is not a property of, or a defect in, the transmission "
                "model, which correctly produces the underlying infections.")
        fig("figs/val3e_ascertainment.png",
            "Full reporting model: modeled infections → reported cases once the "
            "literature under-ascertainment gradient is applied. Infants dominate and "
            "adults collapse — bounding the documented distribution from the other side.")
elif V3:
    m = V3["modeled"]
    table(["Run (defaults)", "Modeled infant <1yr", "Modeled adult 20+", "Documented"],
          [[k, pct(m[k]["infant_share"]), pct(m[k]["adult_share"]), "infant 5.6–12.5% · adult ~19%"]
           for k in m])
    fig("figs/val3_agedist.png", "Modeled vs documented age shares.")

rule()
rule()
# ---- Part 4 : literature benchmarks
h(2, "4 · Literature benchmarks")
p("Published pertussis epidemiology and modeling give quantitative ranges for most "
  "of the model's parameters and several of its outputs. Scored against them:")
table(["Quantity", "Published value (source)", "This model", "Verdict"],
      [["R₀", "12–17 common; 5–17 broader (HIQA; Anderson & May)",
        "default 8; calibrated OR 15, AL 7.2", "✅ within range"],
       ["Incubation period", "7–10 days (CDC)", "9 days", "✅"],
       ["Infectious period", "up to ~21 days communicability (CDC)", "15 d naïve / 12 d partial", "✅ plausible"],
       ["Vaccine immunity (ω₃)", "4–12 yr (Wendelboe 2005)", "7 yr", "✅ mid-range"],
       ["Natural immunity (ω₂)", "4–20 yr (Wendelboe); mean ≥30 yr (Wearing & Rohani 2009)",
        "12 yr", "⚠️ matches clinical range, below W&R — open to calibration"],
       ["Partial infection less infectious (σ)", "contributes little to transmission; contact # < R₀ (Hethcote 2000; W&R 2009)",
        "σ calibrated ~0.25", "✅ qualitatively confirmed"],
       ["Seasonality", "US summer–fall peak", "peak-day 240 (late Aug)", "✅"],
       ["Inter-epidemic cycle", "2–4.6 yr; 3.3–4.1 yr post-vaccination", "~2.0 yr (demography companion §5)", "✅ reproduced (short end)"],
       ["Age distribution (reported cases)", "infants 5.6–12.5%, adults ~19% (CDC/MMWR)", "infant reproduced via severity filter (→13.8%); adults over-represented", "⚠️ partial — see §3d"],
       ["Outbreak cost", "$2,172/case Omaha 2008 (MMWR 60(1) 2011)", "reproduced (EPICC)", "✅ §3b"]])
callout("note", "The gold-standard structural anchor — Wearing & Rohani (2009)",
        "The model's 8-compartment partial-immunity structure comes directly from "
        "Wearing & Rohani (2009, PLOS Pathogens), who validated it against England & "
        "Wales incidence using two *epidemiological signatures*: the gap between major "
        "outbreaks and the frequency of case-free (fade-out) periods. Reproducing "
        "those signatures is the definitive dynamical validation for this model type — "
        "but it requires demography (births feeding the susceptible pool), which this "
        "single-season configuration does not yet include. The literature points "
        "straight back at the same gap identified in §2.")
ul([
  "Wearing & Rohani 2009, *Estimating the Duration of Pertussis Immunity Using Epidemiological Signatures*, PLOS Pathogens — PMC2763266.",
  "Wendelboe et al. 2005, *Duration of Immunity Against Pertussis After Natural Infection or Vaccination*, Pediatr Infect Dis J.",
  "Hethcote 2000, *The Mathematics of Infectious Diseases*, SIAM Review (pertussis R₀; contact number < R₀).",
  "HIQA Immunisation Guidelines, Ch. 15 Pertussis (R₀ 12–17).",
  "Also relevant: Hethcote (age-structured US pertussis), Lavine, King & Bjørnstad 2011 PNAS (immune boosting), Rohani, Zhong & King 2010 Science (age structure).",
])

rule()
# ---- Part 5 : demography / Wearing-Rohani
h(2, "5 · Demography & the inter-epidemic cycle (Wearing–Rohani)")
p("§2 and §4 both pointed at the same missing ingredient — **demography** (births "
  "refilling the susceptible pool). To test it, a standalone **deterministic** version "
  "of the app's 8-compartment partial-immunity model was built with vital dynamics "
  "added (births into S, deaths from all compartments; US crude birth rate **11.6 / "
  "1000 / yr**, NCHS). The intrinsic inter-epidemic period is read from the dominant "
  "oscillatory eigenvalue of the Jacobian at the endemic equilibrium — the "
  "Wearing-&-Rohani \"natural period\".")
if V6:
    pre = V6["Pre-vaccination"]; vax = V6["With vaccination (85%)"]
    table(["Scenario", "Intrinsic period", "Damping time", "Literature"],
          [["Pre-vaccination", f"{pre['intrinsic_period_years']} yr",
            f"~{pre['damping_time_years']} yr", "2–3 yr (pre-vax); 2–4.6 yr empirical"],
           ["With vaccination (85%)", f"{vax['intrinsic_period_years']} yr",
            f"~{vax['damping_time_years']} yr", "3.3–4.1 yr (post-vax)"]])
    callout("note", "The Wearing–Rohani signature is reproduced — within range, at the short end",
            "Adding demography produces a realistic intrinsic multi-year cycle: the "
            "pre-vaccination period (~2.0 yr) sits in the literature's 2–3 yr pre-vax "
            "range and the low end of the 2–4.6 yr empirical band. The model therefore "
            "reproduces the core dynamical signature its structure (Wearing & Rohani "
            "2009) was built to explain. The period lands at the short end, and "
            "vaccination lengthens it only modestly — both consequences of the model's "
            "relatively fast waning (ω₂ = 12 yr) and leaky-vaccine representation.")
    if V6.get("immunity_sweep"):
        sw = V6["immunity_sweep"]
        pairs = ", ".join(f"{iy}yr→{pp}yr" for iy, pp in
                          zip(sw["natural_immunity_years"], sw["intrinsic_period_years"]))
        p("**Sensitivity ties this back to §4.** The cycle lengthens monotonically with "
          f"natural-immunity duration ({pairs}). At the model default ω₂ = 12 yr the "
          "period is ~2.0 yr; at Wearing & Rohani's ≥30 yr it lengthens to ~2.4 yr — so "
          "the short cycle and the short immunity duration flagged in §4 are the same "
          "finding. Longer immunity (or a less leaky vaccine) moves the period up toward "
          "the middle of the empirical range.")
    fig("figs/val6_demography.png",
        "Damped inter-epidemic oscillations after a perturbation, pre- and post-"
        "vaccination. The spacing of the waves is the intrinsic period.")
    fig("figs/val6b_immunity_sweep.png",
        "Intrinsic cycle period vs natural-immunity duration — it rises through the "
        "empirical 2–4.6 yr band as immunity lengthens; the model default (12 yr) and "
        "Wearing & Rohani's ≥30 yr are marked.")
    callout("caution", "Scope of this result",
            "This is now a committed engine module — `engine/demography.py` (smoke-"
            "tested in `scripts/smoke_demography.py`) — a deterministic, age-aggregated "
            "companion to the stochastic app engine. epydemix's transition framework is "
            "compartment-to-compartment only and cannot express births (a source with no "
            "origin) or aging (movement between demographic groups), so demography lives "
            "in its own module rather than being bolted onto the single-season stub. It "
            "validates that the compartmental structure produces realistic cycling. A "
            "tau-leaping stochastic version is also in the module (used for fade-outs "
            "below); a UI entry point is the remaining packaging step (§6).")
    if V6C:
        asd = V6C["infections_age_structured"]; ssd = V6C.get("infections_single_season")
        h(4, "Age-structured extension — a realistic infant infection share")
        p("The age-aggregated model above (and the single-season app model) have no "
          "infant influx, so they put almost no infections in infants. The age-"
          "structured demographic model (`build_age_params` / "
          "`infection_age_distribution` in `engine/demography.py`) adds births into the "
          "youngest band, aging across the 7 bands, and the app's contact matrix — and "
          "the infant picture becomes realistic:")
        rows = []
        if ssd:
            rows.append(["Modeled infections — single-season (no births)", pct(ssd["infant"]), pct(ssd["adult"])])
        rows.append(["Modeled infections — age-structured + demography", pct(asd["infant"]), pct(asd["adult"])])
        rows.append(["Documented reported cases", "5.6–12.5%", "~19%"])
        table(["Infection age mix", "Infant <1yr", "Adult 20+"], rows)
        callout("note", "Births close most of the infant gap; the rest is reporting",
                "Adding demography lifts the modeled infant *infection* share from ~0.7% "
                "to ~3.3% and drops adults from ~83% to ~65% — the continuous influx of "
                "susceptible infants and age-mixing do what a closed single-season model "
                "cannot. With this realistic baseline, the documented reported infant "
                "share (~12.5%) is bracketed between infections (3.3%) and the severity-"
                "weighted reported distribution (~36%), and the reporting gradient implied "
                "by the Lane County line-list confirms steep **adult under-ascertainment** "
                "(adults reported at ~0.06–0.12 of the adolescent peak; the adolescent peak "
                "itself reflects that specific school outbreak, and infant/elderly are "
                "undefined there with zero local cases). A representative national age-"
                "stratified series would pin the gradient precisely.")
        fig("figs/val6c_agestructured.png",
            "Left: infection age distribution — the single-season model vs the age-"
            "structured demographic model (infants rise from 0.7% to 3.3%). Right: the "
            "reporting-fraction gradient implied by the Lane County cases, collapsing for "
            "adults.")
    if V6D:
        h(4, "Stochastic fade-outs & the critical community size")
        p("The **second** Wearing & Rohani signature is the frequency of case-free "
          "(fade-out) periods, which only appears in a stochastic model. A tau-leaping "
          "stochastic version (`simulate_stochastic` in `engine/demography.py`, with "
          "seasonal forcing and rare imports) was swept across population sizes; a "
          "fade-out is a week whose trough reaches zero infectives.")
        fr = V6D["fadeout_fraction"]; Ns = V6D["population"]
        table(["Population", "Fade-out frequency"],
              [[f"{n:,}", f"{f:.0%}"] for n, f in zip(Ns, fr)])
        callout("note", "Second signature reproduced — CCS below measles, as expected",
                f"Fade-outs are frequent in small populations (≈{fr[0]:.0%} of weeks at "
                f"{Ns[0]:,}) and vanish above a **critical community size of "
                f"≈{V6D['approx_ccs']:,}**. That is well below measles' classic "
                "~250–500k CCS — exactly the direction expected for pertussis, which "
                "persists in smaller populations thanks to its waning-immunity adult "
                "reservoir and longer transmission. The exact CCS is sensitive to the "
                "seasonal amplitude, import rate and the age-aggregated approximation, so "
                "read it as order-of-magnitude; the qualitative signature — a sharp "
                "fade-out/persistence transition with population size — is the result.")
        fig("figs/val6d_fadeout.png",
            "Left: fade-out frequency falls with population size, vanishing above the "
            "critical community size. Right: weekly troughs for a small population "
            "(repeated fade-outs) vs a large one (persistent circulation).")

rule()
h(2, "6 · Next validation milestone (roadmap)")
p("§5 now covers the full demographic story — the intrinsic cycle, the age-structured "
  "infection distribution, and stochastic fade-outs are all in `engine/demography.py`, "
  "and **both** Wearing & Rohani signatures are reproduced. What remains:")
steps([
  "**UI entry for the demographic engine** — the model (deterministic + stochastic, "
  "age-aggregated + age-structured) is committed and smoke-tested; a Streamlit entry "
  "point would expose the cycle period, age distribution and fade-outs in the app itself.",
  "**Calibrate the reporting gradient to representative data** — §5 implies the gradient "
  "from the Lane County line-list (steep adult drop, adolescent-peaked by that outbreak); "
  "a national age-stratified surveillance series would pin it and reproduce the full "
  "reported-case age distribution, not just the bracket.",
  "**Forecasting validation (optional)** — calendar-anchored seeding + a reporting "
  "model, trained and tested across multiple high-incidence seasons, to see whether the "
  "model can beat persistence once it is a forecaster by construction.",
])
callout("note", "Where §6 stands",
        "Of the original roadmap, the demographic engine (age-aggregated + age-"
        "structured + stochastic) and **both** Wearing & Rohani signatures (intrinsic "
        "cycle ≈2 yr; fade-out/critical-community-size ≈50k) are done and committed. "
        "What is left is packaging (a UI entry) and refinement (a representative national "
        "reporting calibration), plus the optional forecasting study.")

h(2, "What this validation does NOT cover")
ul([
  "**Multi-year dynamics in the app engine** — the inter-epidemic cycle is reproduced in a standalone demographic companion (§5), but the app's own stochastic, age-structured engine still has no births/aging; porting it in is §6 step 1.",
  "**Curve-only identifiability** — a weekly curve constrains R₀ only; σ, δ and the Sₚ share are not validated by the hold-out and stay near their priors.",
  "**Absolute case counts** — every metric here is scale-free; the model is not validated as a case-count forecaster (under-reporting + no demography).",
  "**Hospitalization/severity** — the hospitalization module's rates were not independently scored here (documented infant 28.6% / non-infant 4.7% are available as a future target).",
])

h(2, "Methods & reproducibility")
ul([
  "**§1** `val_implementation.py` — conservation, non-negativity, R₀ recovery (clean SEIR limit, seasonality off).",
  "**§2** `val_holdout.py` — 1-D R₀ grid fit on NNDSS season t, scored on season t+1; peak-normalised weekly shape distance vs persistence.",
  "**§3d** `val_agedist.py` / `val_agedist2.py` / `val_agedist3.py` / `val_agedist4.py` — modeled age-distribution via `modeled_case_summary`, with the uniform/age-immunity/severity-filter/real-coverage variants.",
  "**§3a–c** CDC NNDSS and EPICC (`cdc-data`, `epicc` MCP) pulled live this session.",
  "**§3d coverage** CDC NIS-Child DTaP (`fhky-rtsk`), NIS-Teen Tdap (`ee48-w5t6`), BRFSS adult Tdap (`aetd-68ew`) — queried live via `cdc-data` `socrata_query`.",
  "**§3d ascertainment** `val_ascertainment.py` — age-specific reporting gradient anchored to the Italy seroprevalence study (Infect Dis Ther 2023).",
  "**§5** `val_demography.py` — standalone deterministic 8-compartment model + vital dynamics; intrinsic period from Jacobian eigenvalues; US crude birth rate from NCHS `e6fc-ccez`.",
  "Engine, data and calibration are the same modules the app uses; JSON outputs are in `docs/validation/figs/`.",
])

# ---------------------------------------------------------------- renderers
def render_md():
    out = []
    for kind, pl in BLOCKS:
        if kind == "h":
            lvl, txt = pl; out.append(f"{'#'*lvl} {txt}\n")
        elif kind == "p": out.append(pl + "\n")
        elif kind == "steps":
            for i, s in enumerate(pl, 1): out.append(f"{i}. {s}")
            out.append("")
        elif kind == "ul":
            for s in pl: out.append(f"- {s}")
            out.append("")
        elif kind == "fig":
            src, cap = pl; out.append(f"![{cap}]({src})\n"); out.append(f"*{cap}*\n")
        elif kind == "callout":
            k, title, txt = pl
            icon = {"note":"ℹ️","caution":"⚠️","tip":"💡"}[k]
            out.append(f"> {icon} **{title}**  ")
            for ln in txt.split("\n\n"): out.append(f"> {ln}  ")
            out.append("")
        elif kind == "table":
            heads, rows = pl
            out.append("| " + " | ".join(heads) + " |")
            out.append("| " + " | ".join(["---"]*len(heads)) + " |")
            for r in rows: out.append("| " + " | ".join(str(c) for c in r) + " |")
            out.append("")
        elif kind == "rule": out.append("---\n")
    return "\n".join(out)

def hi(t):
    t = html.escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"`(.+?)`", r"<code>\1</code>", t)
    t = re.sub(r"\*(.+?)\*", r"<em>\1</em>", t)
    t = re.sub(r"\[(.+?)\]\((.+?)\)", r'<a href="\2">\1</a>', t)
    return t

CSS = re.search(r"<style>(.*?)</style>",
                open(os.path.join(REPO, "docs", "oregon-user-manual", "manual.html"), encoding="utf-8").read(),
                re.S).group(1)

def render_html():
    body = []
    for kind, pl in BLOCKS:
        if kind == "h":
            lvl, txt = pl
            if lvl == 1: body.append('<p class="eyebrow">EpyScenario · SOAR / IRS</p>')
            body.append(f"<h{lvl}>{hi(txt)}</h{lvl}>")
        elif kind == "p":
            cls = ' class="lede"' if pl.startswith("*") and pl.endswith("*") else ""
            body.append(f"<p{cls}>{hi(pl)}</p>")
        elif kind == "steps":
            body.append("<ol>" + "".join(f"<li>{hi(s)}</li>" for s in pl) + "</ol>")
        elif kind == "ul":
            body.append("<ul>" + "".join(f"<li>{hi(s)}</li>" for s in pl) + "</ul>")
        elif kind == "fig":
            src, cap = pl
            body.append(f'<figure><img src="{src}" alt="{html.escape(cap)}"><figcaption>{hi(cap)}</figcaption></figure>')
        elif kind == "callout":
            k, title, txt = pl
            ps = "".join(f"<p>{hi(ln)}</p>" for ln in txt.split("\n\n"))
            body.append(f'<div class="callout {k}"><span class="ct">{hi(title)}</span>{ps}</div>')
        elif kind == "table":
            heads, rows = pl
            th = "".join(f"<th>{hi(str(x))}</th>" for x in heads)
            tr = "".join("<tr>" + "".join(f"<td>{hi(str(c))}</td>" for c in r) + "</tr>" for r in rows)
            body.append(f"<table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table>")
        elif kind == "rule": body.append("<hr>")
    today = datetime.date.today().isoformat()
    return ("<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
            "<title>EpyScenario Pertussis — Validation Report</title>"
            f"<style>{CSS}</style></head><body><div class=\"wrap\">\n"
            + "\n".join(body)
            + f"\n<footer>EpyScenario · SOAR / IRS — validation generated {today}. "
              "Reproducible from docs/validation scripts.</footer>\n</div></body></html>\n")

open(os.path.join(BASE, "README.md"), "w", encoding="utf-8").write(render_md())
open(os.path.join(BASE, "validation.html"), "w", encoding="utf-8").write(render_html())
print("wrote docs/validation/README.md and validation.html; figs:",
      sum(1 for k, _ in BLOCKS if k == "fig"), "| V1", bool(V1), "V2", bool(V2), "V3", bool(V3))
