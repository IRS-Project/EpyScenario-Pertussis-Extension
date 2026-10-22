# EpyScenario — Pertussis Extension (IRS)

A bespoke **pertussis (whooping cough) modeling tool** for the IRS/SOAR project,
built as a fork of [EpyScenario](https://scenario.epydemix.org/)
([`ngozzi/epydemix-dashboard`](https://github.com/ngozzi/epydemix-dashboard)),
which is itself built on [Epydemix](https://epydemix.org).

This fork adds an age-structured, two-track **SEIRS partial-immunity** pertussis
model with real-data calibration, a vaccination planner, demographic analysis, and
a full documentation set — on top of the upstream scenario platform.

> **Private IRS tool.** Licensed **GPL v3** (inherited from upstream); the LICENSE
> and upstream attribution are preserved. See [Credits](#credits--upstream).

---

## What this fork adds

- **Pertussis model** — 8-compartment SEIRS with a naive track (S→E→I→R) and a
  partial-immunity track (Sₚ→Eₚ→Iₚ→Rₚ). Force of infection `λ = β·(I + σ·Iₚ)`
  drives S→E, `δλ` drives Sₚ→Eₚ; waning (ω₁/ω₂/ω₃) and vaccination S→Sₚ. 7 age
  bands. Closed population (no births/deaths — susceptibility is replenished by
  waning). See [`docs/Pertussis_Model_Diagram.pdf`](docs/Pertussis_Model_Diagram.pdf).
- **Calibration** — ABC-SMC posterior and posterior-predictive fits to weekly or
  cumulative CDC NNDSS pertussis curves (weekly series derived from the YTD
  cumulative column so annual totals match), plus an optional hospitalization overlay.
- **Demography / vital dynamics** — an age-structured demographic engine
  (Wearing–Rohani style), stochastic fade-outs and critical community size.
- **Vaccination planner** — age-stratified DTaP/Tdap scheduling with persistent
  campaign and coverage-by-age stores.
- **In-app provenance** — a Model References & Provenance page with
  data→parameter calculators and tooltips explaining how each default is derived.
- **Documentation** — Oregon (Lane County) and Alabama user manuals, a validation
  report, UAT evidence, and release notes under [`docs/`](docs/).

## Installation

Requires Python 3.11+.

1. **Clone the repository** (private — you need access to the IRS-Project org):
    ```bash
    git clone https://github.com/IRS-Project/EpyScenario-Pertussis-Extension.git
    cd EpyScenario-Pertussis-Extension
    ```
2. **Create and activate a virtual environment:**
    ```bash
    python -m venv venv
    # Windows (PowerShell):
    venv\Scripts\Activate.ps1
    # macOS / Linux:
    source venv/bin/activate
    ```
3. **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

## Running the dashboard

Use the **supervised launcher** (auto-restarts the Streamlit server and applies the
runtime shims this app needs):

```bash
# Windows (PowerShell)
./run_dashboard.ps1

# macOS / Linux
./run_dashboard.sh
```

Then open **http://localhost:8501** in your browser. (You can also run it directly
with `streamlit run Dashboard.py`, but the launcher is recommended.)

### Using the pertussis model

1. In the setup panel, choose the **SEIRS (Pertussis)** disease preset and a
   **United States** geography (state or county).
2. Set model parameters (sensible pertussis defaults are pre-filled; each has a
   tooltip and a provenance reference).
3. Optionally calibrate to observed NNDSS data, plan a vaccination campaign, or
   open the **Demography** and **Model References & Provenance** pages for deeper
   analysis.

> **Note:** scenarios are exploratory, not forecasts. See the validation report in
> [`docs/`](docs/) for what the model has and hasn't been checked against.

## Repository layout

| Path | Contents |
|------|----------|
| `Dashboard.py`, `pages/` | Streamlit entry point and app pages |
| `engine/` | Simulation, calibration, demography, hospitalization models |
| `ui/` | Setup, parameters, visualisation, vaccination panels |
| `data/` | Observed data, weekly NNDSS sources, campaign/coverage stores |
| `schemas.py`, `constants.py` | Model presets, parameters, defaults |
| `docs/` | Manuals, validation report, model diagram, provenance, release notes |

## Credits & upstream

- **EpyScenario** — the scenario platform this is forked from:
  [`ngozzi/epydemix-dashboard`](https://github.com/ngozzi/epydemix-dashboard) ·
  live demo [scenario.epydemix.org](https://scenario.epydemix.org/)
- **Epydemix** — the underlying epidemic modeling library:
  [epydemix.org](https://epydemix.org)
- Pertussis model structure adapted from **Wearing & Rohani (2009)**.

## License

**GNU General Public License v3.0** — see [`LICENSE`](LICENSE). As a GPL v3 work,
any distributed derivative must remain under GPL v3 with source available.
