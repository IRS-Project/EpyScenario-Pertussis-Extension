# scripts/smoke_demography.py
#
# Smoke test for engine/demography.py (the vital-dynamics companion model).
# Checks the intrinsic inter-epidemic period lands in the literature range and
# lengthens with natural-immunity duration. Run:  python scripts/smoke_demography.py
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from engine.demography import (build_params, intrinsic_period, simulate,
                               build_age_params, infection_age_distribution,
                               simulate_stochastic)


def main() -> int:
    ok = True

    pre = intrinsic_period(build_params(vaccinated_fraction=0.0))["period_years"]
    vax = intrinsic_period(build_params(vaccinated_fraction=0.85))["period_years"]
    print(f"intrinsic period  pre-vax={pre:.2f} yr  vax85%={vax:.2f} yr  (empirical 2-4.6 yr)")
    ok &= 1.5 <= pre <= 5.0           # within the plausible pertussis cycle band

    periods = [intrinsic_period(build_params(waning_partial_to_susceptible_years=iy))["period_years"]
               for iy in (4, 12, 30)]
    print(f"period vs immunity (4/12/30 yr): {[round(x, 2) for x in periods]}")
    ok &= periods[0] < periods[1] < periods[2]   # monotone increasing

    sim = simulate(build_params(), years=30)
    print(f"simulate: {len(sim['incidence'])} points, mean incidence {sim['incidence'].mean():.5f}")
    ok &= sim["incidence"].min() >= 0 and sim["incidence"].mean() > 0

    # age-structured: synthetic contact matrix (assortative), check a valid dist
    M = 2.0 * np.ones((7, 7)) + 8.0 * np.eye(7)
    ad = infection_age_distribution(build_age_params(M))
    print(f"age-structured infection dist sums to {ad['dist'].sum():.3f}; infant={ad['infant']:.1%}")
    ok &= abs(ad["dist"].sum() - 1.0) < 1e-6 and ad["infant"] > 0

    # stochastic fade-outs: a small population fades more than a large one
    sm = simulate_stochastic(build_params(), population=3000, years=40, seed=1)["fadeout_fraction"]
    lg = simulate_stochastic(build_params(), population=500000, years=40, seed=1)["fadeout_fraction"]
    print(f"fade-out: N=3k {sm:.1%}, N=500k {lg:.1%}")
    ok &= sm > lg

    print("SMOKE_OK" if ok else "SMOKE_FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
