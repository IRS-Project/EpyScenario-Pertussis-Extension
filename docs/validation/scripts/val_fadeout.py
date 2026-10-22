# -*- coding: utf-8 -*-
"""§6 step 2: stochastic fade-outs (Wearing-Rohani second signature). Sweep
population size; measure the fraction of weeks whose trough hits zero infectives.
Fade-outs vanish above a critical community size (CCS)."""
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
from engine.demography import build_params, simulate_stochastic

OUT = os.path.join(REPO, "docs", "validation", "figs")
p = build_params()
Ns = [2000, 5000, 10000, 20000, 50000, 100000, 300000, 1000000]
SEEDS = [1, 2, 3]

fade = []
for N in Ns:
    vals = [simulate_stochastic(p, population=N, years=50, seed=s)["fadeout_fraction"] for s in SEEDS]
    fade.append(float(np.mean(vals)))
    print(f"N={N:>9,}  fade-out {fade[-1]:.1%}")

# CCS ~ where fade-out drops below 5%
ccs = next((N for N, f in zip(Ns, fade) if f < 0.05), None)
print("approx critical community size (<5% fade-out):", f"{ccs:,}" if ccs else "n/a")

# illustrative weekly-trough series: small (fades) vs large (persists)
small = simulate_stochastic(p, population=10000, years=40, seed=7)["weekly_trough_infectives"]
large = simulate_stochastic(p, population=500000, years=40, seed=7)["weekly_trough_infectives"]

res = {"population": Ns, "fadeout_fraction": [round(f, 3) for f in fade],
       "approx_ccs": ccs, "seeds": SEEDS,
       "note": "fade-out = fraction of weeks whose trough reaches zero I+Ip after burn-in; "
               "measles CCS ~250-500k, pertussis generally lower (waning-immunity reservoir)"}
json.dump(res, open(os.path.join(OUT, "val6d_fadeout.json"), "w"), indent=2)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.2))
a1.semilogx(Ns, np.array(fade) * 100, "o-", color="#0e7a78", lw=2)
if ccs:
    a1.axvline(ccs, color="#c2410c", ls="--", lw=1.2, label=f"CCS ≈ {ccs:,} (<5%)")
a1.axvline(300000, color="#475569", ls=":", lw=1.2, label="measles CCS ~300k")
a1.set_xlabel("Population size"); a1.set_ylabel("Fade-out frequency (% of weeks)")
a1.set_title("Fade-outs vanish above a critical community size", weight="bold", fontsize=10.5)
a1.legend(fontsize=8, frameon=False); a1.grid(alpha=.3, which="both")
for s in ("top", "right"): a1.spines[s].set_visible(False)

a2.plot(np.arange(len(small)) / 52, small, color="#c2410c", lw=0.8, label="N = 10,000 (fades out)")
a2.plot(np.arange(len(large)) / 52, large, color="#0e7a78", lw=0.8, label="N = 500,000 (persists)")
a2.set_xlabel("Year"); a2.set_ylabel("Weekly trough infectives")
a2.set_title("Weekly troughs: small vs large population", weight="bold", fontsize=10.5)
a2.legend(fontsize=8, frameon=False); a2.grid(alpha=.3)
for s in ("top", "right"): a2.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "val6d_fadeout.png"), dpi=140)
print("\nwrote val6d_fadeout.json and .png")
