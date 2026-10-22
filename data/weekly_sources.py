# data/weekly_sources.py
#
# Public weekly pertussis incidence time series, for use as a *shape* target in
# calibration (growth rate, peak timing, seasonality — the things a cumulative
# cross-tab cannot constrain).
#
# Source: CDC National Notifiable Diseases Surveillance System (NNDSS) "Weekly
# Data" table, Socrata dataset x9gk-5huc on data.cdc.gov. Fields used:
#   label   = disease name ("Pertussis")
#   states  = reporting area (state name; UPPERCASE in 2022-2024, Title Case in
#             2025-2026 — matched case-insensitively so both eras merge cleanly)
#   year    = MMWR year
#   week    = MMWR week
#   m1      = current-week case count (absent  => 0 that week; m3 is cumulative)
#
# Works for every U.S. state, DC, and the national total. Results are cached to
# data/observed/weekly_cache/<area>.csv so the app works offline after a first
# fetch and does not hit the network on every Streamlit rerun.

from __future__ import annotations
import json
import os
import re
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

CDC_HOST = "https://data.cdc.gov"
CDC_DATASET = "x9gk-5huc"
CDC_LABEL = "Pertussis"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "observed", "weekly_cache")

# National total is stored under the reporting area "TOTAL" in this dataset.
NATIONAL_AREA = "TOTAL"
NATIONAL_LABEL = "United States (national total)"

US_STATES = [
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
    "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia",
    "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
    "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
    "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
    "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
    "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island",
    "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont",
    "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming",
]

# Selectable reporting areas: national total first, then the states + DC.
WEEKLY_AREAS = [NATIONAL_LABEL] + US_STATES


def _area_query_value(area: str) -> str:
    """Map a friendly area label to the dataset's `states` value (matched
    case-insensitively, so era-specific casing does not matter)."""
    if area == NATIONAL_LABEL:
        return NATIONAL_AREA
    return area


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def _soql(params: dict, timeout: int = 30) -> list:
    query = urllib.parse.urlencode(params)
    url = f"{CDC_HOST}/resource/{CDC_DATASET}.json?{query}"
    req = urllib.request.Request(url, headers={"User-Agent": "epyscenario/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_cdc_weekly_pertussis(area: str, years=None, timeout: int = 30) -> pd.DataFrame:
    """Fetch the weekly pertussis series for a reporting area from CDC NNDSS.

    Returns a DataFrame [year, week, cases] sorted chronologically.

    IMPORTANT: the dataset's `m1` column (current-week new cases) systematically
    UNDER-counts — CDC adds many cases straight to the year-to-date cumulative
    without attributing them to a specific week (reporting lag). So this derives
    the weekly series from the YTD cumulative column `m3` (week-over-week
    differences), which matches the authoritative annual totals (e.g. Oregon 2025
    = 1,560, not the ~390 you get summing `m1`). Falls back to `m1` only if `m3`
    is entirely absent for a year. Raises on a network/HTTP error."""
    q_area = _area_query_value(area).replace("'", "''")
    where = f"label='{CDC_LABEL}' AND upper(states)=upper('{q_area}')"
    if years:
        yrs = ",".join(f"'{int(y)}'" for y in years)
        where += f" AND year in({yrs})"
    rows = _soql({
        "$select": "year,week,m1,m3",   # m1 = current week; m3 = cumulative YTD
        "$where": where,
        "$order": "year,week",
        "$limit": "5000",
    }, timeout=timeout)

    def _num(raw):
        try:
            return float(raw) if raw not in (None, "", ".", "-") else float("nan")
        except (TypeError, ValueError):
            return float("nan")

    recs = []
    for r in rows:
        try:
            y, w = int(r["year"]), int(r["week"])
        except (KeyError, TypeError, ValueError):
            continue
        recs.append((y, w, _num(r.get("m1")), _num(r.get("m3"))))

    raw = pd.DataFrame(recs, columns=["year", "week", "m1", "m3"])
    if raw.empty:
        return pd.DataFrame(columns=["year", "week", "cases"])
    # Eras do not overlap, but guard against any duplicate (year, week).
    raw = raw.groupby(["year", "week"], as_index=False).agg(m1=("m1", "max"), m3=("m3", "max"))

    parts = []
    for y, g in raw.groupby("year", sort=True):
        g = g.sort_values("week")
        cum = g["m3"].ffill()                     # YTD cumulative, forward-filled
        if cum.notna().any():
            cum = cum.fillna(0.0)
            weekly = cum.diff()
            weekly.iloc[0] = cum.iloc[0]          # first week = first cumulative
            # NOT clipped: sum of diffs == final YTD exactly (a downward data
            # revision shows as a small negative week; calibration clips those).
        else:                                     # no YTD column for this year -> use m1
            weekly = g["m1"].fillna(0.0).clip(lower=0.0)
        parts.append(pd.DataFrame({"year": int(y), "week": g["week"].astype(int),
                                   "cases": weekly.round().astype(int)}))
    return pd.concat(parts, ignore_index=True).sort_values(["year", "week"]).reset_index(drop=True)


def load_weekly_pertussis(area: str, years=None, refresh: bool = False):
    """Cached wrapper around fetch_cdc_weekly_pertussis.

    Returns (DataFrame[year, week, cases], source_str). Falls back to the cached
    CSV if the network fetch fails; raises only if there is no cache either."""
    slug = _slug(area)
    if years:
        slug += "_" + "_".join(str(int(y)) for y in sorted(years))
    path = os.path.join(CACHE_DIR, f"{slug}.csv")

    if not refresh and os.path.exists(path):
        return pd.read_csv(path), "cache"

    try:
        df = fetch_cdc_weekly_pertussis(area, years)
        os.makedirs(CACHE_DIR, exist_ok=True)
        df.to_csv(path, index=False)
        return df, "CDC NNDSS (network)"
    except Exception as exc:  # network/HTTP/parse failure
        if os.path.exists(path):
            return pd.read_csv(path), f"cache (refresh failed: {exc})"
        raise


def weekly_series_to_array(df: pd.DataFrame):
    """Flatten a [year, week, cases] frame to a chronological cases array."""
    if df is None or df.empty:
        return []
    d = df.sort_values(["year", "week"])
    return [float(c) for c in d["cases"].to_numpy()]


# Seasonality-amplitude categories -> trough/peak ratio (matches engine SEASONALITY_OPTIONS
# and the schema options for seasonality_amplitude).
_AMPLITUDE_CATS = {"Strong": 0.5, "Moderate": 0.65, "Medium": 0.75, "Low": 0.9, "None": 1.0}


def estimate_seasonality(area: str, years=None, smooth: int = 5) -> dict | None:
    """Estimate the seasonal peak timing and amplitude from multi-year weekly
    NNDSS pertussis data.

    Each year's weekly counts are normalised by that year's mean (so an epidemic
    year does not dominate the *shape*), then averaged by MMWR week across years
    to give a typical seasonal cycle. The peak week of the (smoothed) average
    cycle -> seasonality_peak_day (day-of-year, since START_DATE is Jan 1); the
    trough/peak ratio -> the nearest seasonality-amplitude category.

    Returns a dict with the estimate and the per-year / average frames for
    plotting, or None if no data. Pins seasonality empirically rather than
    letting R0 absorb seasonal misfit during calibration.
    """
    df, src = load_weekly_pertussis(area, years=years)
    if df is None or df.empty:
        return None

    per_rows = []
    for y, g in df.groupby("year"):
        g = g.sort_values("week")
        m = float(g["cases"].mean())
        norm = (g["cases"] / m) if m > 0 else g["cases"] * 0.0
        for wk, c in zip(g["week"].astype(int), norm):
            per_rows.append({"year": int(y), "week": int(wk), "norm": float(c)})
    per_year = pd.DataFrame(per_rows)
    if per_year.empty:
        return None

    avg = per_year.groupby("week", as_index=False)["norm"].mean().sort_values("week")
    if smooth and smooth > 1:
        avg["smooth"] = avg["norm"].rolling(smooth, center=True, min_periods=1).mean()
    else:
        avg["smooth"] = avg["norm"]

    peak_week = int(avg.loc[avg["smooth"].idxmax(), "week"])
    peak_day = int(min(365, max(1, round((peak_week - 0.5) * 7))))
    # Robust trough/peak (10th/90th percentile of the smoothed cycle) so a single
    # near-zero week in sparse data does not force the amplitude to "Strong".
    sm = avg["smooth"].to_numpy(dtype=float)
    peak = float(np.nanpercentile(sm, 90))
    trough = float(np.nanpercentile(sm, 10))
    ratio = (trough / peak) if peak > 0 else 1.0
    label = min(_AMPLITUDE_CATS, key=lambda k: abs(_AMPLITUDE_CATS[k] - ratio))
    # Spread of each year's own peak week — a proxy for how consistent (trustworthy)
    # the seasonal timing is across years.
    peak_weeks_by_year = (per_year.sort_values("norm")
                          .groupby("year")["week"].last())
    peak_week_spread = float(peak_weeks_by_year.std()) if len(peak_weeks_by_year) > 1 else 0.0

    return {
        "area": area,
        "years": sorted(per_year["year"].unique().tolist()),
        "n_years": int(per_year["year"].nunique()),
        "peak_week": peak_week,
        "peak_day": peak_day,
        "amplitude_ratio": ratio,
        "amplitude_label": label,
        "peak_week_spread": peak_week_spread,
        "source": src,
        "per_year": per_year,
        "avg": avg,
    }
