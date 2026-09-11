# Descriptive analytics on the merged state table. No ML: rankings, gaps,
# correlations, group comparisons, and the cross-checks the article relies on.
#
# Reads : data/processed/state_merged.csv
#         data/processed/mortality_by_state_year.csv
#         data/raw/cdc/opioid_by_race_2018_2024.csv
# Writes: nothing (prints everything; the article build recomputes what it needs)

from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
m = pd.read_csv(ROOT / "data" / "processed" / "state_merged.csv", dtype={"fips": str})
long = pd.read_csv(ROOT / "data" / "processed" / "mortality_by_state_year.csv", dtype={"fips": str})
race = pd.read_csv(ROOT / "data" / "raw" / "cdc" / "opioid_by_race_2018_2024.csv")

pd.set_option("display.width", 220)
LATEST = int(m["year_latest"].iloc[0])
FIVE = ["SC", "WA", "TN", "NV", "AL"]

# ---------- drug-specificity adjusted burden (2022-2024) ----------
w3 = long[long.year >= LATEST - 2].groupby("fips").agg(
    d3=("deaths", "sum"), adj3=("opioid_adj", "sum"), p3=("population", "sum"),
    od3=("overdose_total", "sum"), spec3=("overdose_specified", "sum"))
w3["rate_3yr"] = w3.d3 / w3.p3 * 1e5
w3["rate_3yr_adj"] = w3.adj3 / w3.p3 * 1e5
w3["unspec_share"] = (w3.od3 - w3.spec3) / w3.od3 * 100
m = m.drop(columns=["rate_3yr"]).merge(
    w3[["rate_3yr", "rate_3yr_adj", "unspec_share"]], on="fips")

m["r_burden"] = m.rate_3yr.rank(ascending=False).astype(int)
m["r_burden_adj"] = m.rate_3yr_adj.rank(ascending=False).astype(int)
m["r_cap_moud"] = m.moud_per_100k.rank(ascending=False).astype(int)
m["r_cap_otp"] = m.otp_per_100k.rank(ascending=False).astype(int)
m["mismatch"] = m.r_cap_moud - m.r_burden_adj
m["deaths_per_moud_facility"] = m.deaths_latest / m.moud_count

print("=== 1. Drug-specificity correction ===")
print(f"national unspecified share 2022-24: "
      f"{(w3.od3.sum()-w3.spec3.sum())/w3.od3.sum()*100:.1f}%  |  "
      f"rate {w3.d3.sum()/w3.p3.sum()*1e5:.1f} -> {w3.adj3.sum()/w3.p3.sum()*1e5:.1f}")
m["rank_move"] = m.r_burden - m.r_burden_adj
mv = m.reindex(m.rank_move.abs().sort_values(ascending=False).index).head(6)
print(mv[["abbrev", "unspec_share", "rate_3yr", "rate_3yr_adj",
          "r_burden", "r_burden_adj", "rank_move"]].round(1).to_string(index=False))
print("worst-reporting states:",
      ", ".join(f"{r.abbrev} {r.unspec_share:.0f}%"
                for r in m.nlargest(5, "unspec_share").itertuples()))

print("\n=== 2. Burden (adjusted, 2022-24) ===")
cols = ["abbrev", "rate_3yr_adj", "rate_3yr", "deaths_latest", "crude_lo", "crude_hi"]
print(m.nlargest(10, "rate_3yr_adj")[cols].round(1).to_string(index=False))

print("\n=== 3. Capacity ===")
print("official OTPs:", int(m.otp_count.sum()), "| self-reported flag:", int(m.otp_selfreport.sum()),
      "| MOUD facilities:", int(m.moud_count.sum()))
print("states with no certified OTP:", ", ".join(m.loc[m.otp_count == 0, "state"]) or "none")
print(m.nsmallest(8, "moud_per_100k")[["abbrev", "moud_per_100k", "otp_count",
                                       "otp_per_100k"]].round(2).to_string(index=False))

print("\n=== 4. Does capacity track burden? ===")
def sp(x, y, label):
    r, p = stats.spearmanr(m[x], m[y])
    print(f"   {label:44s} r={r:+.3f}  p={p:.4f}")
sp("rate_3yr_adj", "moud_per_100k", "adjusted burden vs MOUD/100k")
sp("rate_3yr_adj", "otp_per_100k", "adjusted burden vs OTP/100k (official)")
sp("rate_3yr", "otp_per_100k", "reported burden vs OTP/100k (official)")

print("\n=== 5. The gap: which states have least, relative to burden ===")
g = ["abbrev", "rate_3yr_adj", "r_burden_adj", "moud_per_100k", "r_cap_moud",
     "mismatch", "deaths_per_moud_facility", "otp_count"]
print(m.nlargest(8, "mismatch")[g].round(2).to_string(index=False))
print("\nhighest deaths per MOUD facility:")
print(m.nlargest(8, "deaths_per_moud_facility")[
    ["abbrev", "deaths_latest", "moud_count", "deaths_per_moud_facility",
     "r_burden_adj"]].round(1).to_string(index=False))

# regression residual as the statistician's version of the rank gap
x = m.rate_3yr_adj.values
slope, icpt, rr, pp, se = stats.linregress(x, m.moud_per_100k.values)
m["resid"] = m.moud_per_100k - (icpt + slope * x)
print("\nlargest negative residuals (least capacity for their burden):")
print(m.nsmallest(8, "resid")[["abbrev", "rate_3yr_adj", "moud_per_100k",
                               "resid", "mismatch"]].round(2).to_string(index=False))
print("rank-gap vs residual agreement (Spearman):",
      f"{stats.spearmanr(m.mismatch, -m.resid).statistic:+.3f}")

print("\n=== 6. Medicaid expansion ===")
for col, lbl in [("moud_per_100k", "MOUD facilities /100k"),
                 ("otp_per_100k", "certified OTPs /100k"),
                 ("deaths_per_moud_facility", "deaths per MOUD facility"),
                 ("rate_3yr_adj", "adjusted death rate")]:
    a, b = m.loc[m.expanded_medicaid, col], m.loc[~m.expanded_medicaid, col]
    U, p = stats.mannwhitneyu(a, b)
    print(f"   {lbl:26s} expanded {a.median():7.2f} | non-expanded {b.median():7.2f} | p={p:.4f}")
print("   non-expansion states among the five:",
      ", ".join(sorted(m.loc[m.abbrev.isin(FIVE) & ~m.expanded_medicaid, "abbrev"])))

print("\n=== 7. Socioeconomics (the null) ===")
for x_, y_, lbl in [("mismatch", "median_hh_income", "gap vs median income"),
                    ("mismatch", "pct_uninsured", "gap vs % uninsured"),
                    ("deaths_per_moud_facility", "median_hh_income", "deaths/facility vs income"),
                    ("deaths_per_moud_facility", "pct_uninsured", "deaths/facility vs uninsured")]:
    r, p = stats.spearmanr(m[x_], m[y_])
    print(f"   {lbl:34s} r={r:+.3f}  p={p:.3f}")

print("\n=== 8. Uncertainty: do ranks survive their confidence intervals? ===")
top = m.nlargest(12, "crude_rate")[["abbrev", "crude_rate", "crude_lo", "crude_hi"]]
print(top.round(1).to_string(index=False))
wv = m.loc[m.abbrev == "WV"].iloc[0]
n_overlap = int(((m.crude_lo <= wv.crude_hi) & (m.crude_hi >= wv.crude_lo)).sum() - 1)
print(f"states whose 2024 CI overlaps West Virginia's: {n_overlap}")

print("\n=== 9. Race (national, 2018-2024) ===")
piv = race.pivot(index="race", columns="year", values="rate")
piv["peak_year"] = piv.idxmax(axis=1)
piv["vs_2018_pct"] = (piv[2024] / piv[2018] - 1) * 100
print(piv[[2018, 2021, 2023, 2024, "peak_year", "vs_2018_pct"]].round(1).to_string())
tot = race.groupby("year").deaths.sum()
blk = race[race.race.str.startswith("Black")].set_index("year").deaths / tot * 100
print("Black share of opioid deaths: "
      f"{blk[2018]:.1f}% (2018) -> {blk[2023]:.1f}% (2023) -> {blk[2024]:.1f}% (2024)")

print("\n=== 10. Trend ===")
nat = long.groupby("year").agg(d=("deaths", "sum"), p=("population", "sum"))
nat["rate"] = nat.d / nat.p * 1e5
print(nat.round(1).to_string())
