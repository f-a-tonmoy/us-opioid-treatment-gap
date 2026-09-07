# Descriptive analytics on the merged state table. No ML - rankings, gaps,
# correlations, trends, and the planned self-correction cross-checks.
#
# Reads : data/processed/state_merged.csv, data/processed/mortality_by_state_year.csv
# Writes: nothing (prints everything; the article build consumes these numbers)

from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
m = pd.read_csv(ROOT / "data" / "processed" / "state_merged.csv", dtype={"fips": str})
long = pd.read_csv(ROOT / "data" / "processed" / "mortality_by_state_year.csv", dtype={"fips": str})

pd.set_option("display.width", 200)
LATEST = int(m["year_latest"].iloc[0])

# ---------- 1. rankings ----------
# burden ranks: 1 = highest mortality; capacity ranks: 1 = most facilities per capita
m["r_burden_aa"] = m["aa_rate"].rank(ascending=False).astype(int)
m["r_burden_crude"] = m["crude_rate"].rank(ascending=False).astype(int)
m["r_burden_3yr"] = m["rate_3yr"].rank(ascending=False).astype(int)
m["r_cap_moud"] = m["moud_per_100k"].rank(ascending=False).astype(int)
m["r_cap_otp"] = m["otp_per_100k"].rank(ascending=False).astype(int)

print(f"=== 1. Burden ranking (age-adjusted opioid deaths /100k, {LATEST}) ===")
cols = ["abbrev", "aa_rate", "crude_rate", "rate_3yr", "deaths_latest"]
print("TOP 10:");    print(m.sort_values("aa_rate", ascending=False)[cols].head(10).round(1).to_string(index=False))
print("BOTTOM 5:");  print(m.sort_values("aa_rate")[cols].head(5).round(1).to_string(index=False))

print("\n=== 2. Capacity ranking (MOUD facilities /100k) ===")
cols = ["abbrev", "moud_per_100k", "otp_per_100k", "moud_count", "otp_count"]
print("TOP 10:");    print(m.sort_values("moud_per_100k", ascending=False)[cols].head(10).round(2).to_string(index=False))
print("BOTTOM 10:"); print(m.sort_values("moud_per_100k")[cols].head(10).round(2).to_string(index=False))

# ---------- 3. the gap ----------
# mismatch = capacity rank - burden rank; positive = burden outranks capacity (underserved)
m["mismatch"] = m["r_cap_moud"] - m["r_burden_3yr"]
print("\n=== 3. Burden-vs-capacity gap ===")
print("(burden = 3yr pooled rate rank; capacity = MOUD/100k rank; +mismatch = underserved)")
cols = ["abbrev", "rate_3yr", "r_burden_3yr", "moud_per_100k", "r_cap_moud", "mismatch", "deaths_per_moud_facility"]
print("MOST UNDERSERVED (high burden, low capacity):")
print(m.sort_values("mismatch", ascending=False)[cols].head(10).round(2).to_string(index=False))
print("MOST OVER-PROVISIONED (relative to burden):")
print(m.sort_values("mismatch")[cols].head(8).round(2).to_string(index=False))
print("\nDeaths per MOUD facility (burden per unit of capacity), top 10:")
print(m.sort_values("deaths_per_moud_facility", ascending=False)[["abbrev", "deaths_latest", "moud_count", "deaths_per_moud_facility", "r_burden_3yr"]].head(10).round(1).to_string(index=False))

# ---------- 4. does capacity track burden at all? + socioeconomics ----------
print("\n=== 4. Correlations (n=51) ===")
def corr(x, y, label):
    rs, ps = stats.spearmanr(m[x], m[y])
    rp, pp = stats.pearsonr(m[x], m[y])
    print(f"  {label:55s} Spearman r={rs:+.3f} (p={ps:.4f}) | Pearson r={rp:+.3f} (p={pp:.4f})")
corr("rate_3yr", "moud_per_100k", "burden (3yr rate) vs capacity (MOUD/100k)")
corr("rate_3yr", "otp_per_100k", "burden (3yr rate) vs capacity (OTP/100k)")
corr("mismatch", "median_hh_income", "mismatch vs median household income")
corr("mismatch", "pct_uninsured", "mismatch vs % uninsured")
corr("deaths_per_moud_facility", "median_hh_income", "deaths-per-facility vs median household income")
corr("deaths_per_moud_facility", "pct_uninsured", "deaths-per-facility vs % uninsured")
corr("rate_3yr", "median_hh_income", "burden vs median household income")
corr("rate_3yr", "pct_uninsured", "burden vs % uninsured")

# ---------- 5. time trend ----------
print("\n=== 5. Trend ===")
nat = long.groupby("year").agg(deaths=("deaths", "sum"), pop=("population", "sum"))
nat["rate"] = nat["deaths"] / nat["pop"] * 1e5
print(nat.round(1).to_string())
piv = long.pivot_table(index="fips", columns="year", values="crude_rate")
peak_col = piv[[2021, 2022, 2023]].max(axis=1)
m2 = m.set_index("fips")
m2["chg_peak_to_latest"] = (piv[LATEST] - peak_col) / peak_col * 100
print("\nBiggest decline, peak(2021-23) -> 2024 (%):")
print(m2.sort_values("chg_peak_to_latest")[["abbrev", "chg_peak_to_latest", "aa_rate"]].head(8).round(1).to_string(index=False))
print("Smallest decline / still rising:")
print(m2.sort_values("chg_peak_to_latest", ascending=False)[["abbrev", "chg_peak_to_latest", "aa_rate"]].head(8).round(1).to_string(index=False))

# ---------- 6. self-correction cross-checks ----------
print("\n=== 6. Self-correction cross-checks ===")
m["shift_crude_aa"] = m["r_burden_crude"] - m["r_burden_aa"]
m["shift_year_3yr"] = m["r_burden_aa"] - m["r_burden_3yr"]
m["shift_otp_moud"] = m["r_cap_otp"] - m["r_cap_moud"]
print("a) crude vs age-adjusted burden rank shifts >3:")
print(m[m["shift_crude_aa"].abs() > 3][["abbrev", "r_burden_crude", "r_burden_aa", "shift_crude_aa"]].to_string(index=False))
print("b) latest-year vs 3yr-pooled burden rank shifts >5:")
print(m[m["shift_year_3yr"].abs() > 5][["abbrev", "r_burden_aa", "r_burden_3yr", "shift_year_3yr", "aa_rate", "rate_3yr"]].round(1).to_string(index=False))
print("c) OTP vs MOUD capacity rank shifts >15:")
print(m[m["shift_otp_moud"].abs() > 15][["abbrev", "r_cap_otp", "r_cap_moud", "otp_count", "moud_count"]].to_string(index=False))
print("d) DC vs state median:", f"DC aa_rate={m.loc[m.abbrev=='DC','aa_rate'].iloc[0]}, median={m['aa_rate'].median():.1f}")
small = m[m["acs_pop"] < 1.2e6][["abbrev", "deaths_latest", "aa_rate", "rate_3yr", "moud_count"]]
print("e) small states (pop <1.2M) - volatility watch:")
print(small.round(1).to_string(index=False))
