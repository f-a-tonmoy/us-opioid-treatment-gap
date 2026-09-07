# Clean and merge the three sources into one state-level table.
#
# Inputs:
#   data/raw/cdc/*.txt|*.csv          CDC WONDER export (user-provided, manual)
#   data/raw/samhsa/facilities_sa_raw.json
#   data/raw/census/acs2023_detailed.json, acs2023_profile.json
# Outputs (written only when the CDC file is present):
#   data/processed/state_merged.csv            one row per state (latest year)
#   data/processed/mortality_by_state_year.csv long state x year table (trends)
#
# If the CDC file is missing, the script still builds and prints the
# capacity + ACS half (join-key audit included) so that part stays verified.
#
# Join key: state FIPS code. The STATES table below is the source of truth;
# every source must match it exactly (50 states + DC = 51 rows) or we stop.

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"

# fips -> (name, USPS abbrev): 50 states + DC
STATES = {
    "01": ("Alabama", "AL"), "02": ("Alaska", "AK"), "04": ("Arizona", "AZ"),
    "05": ("Arkansas", "AR"), "06": ("California", "CA"), "08": ("Colorado", "CO"),
    "09": ("Connecticut", "CT"), "10": ("Delaware", "DE"), "11": ("District of Columbia", "DC"),
    "12": ("Florida", "FL"), "13": ("Georgia", "GA"), "15": ("Hawaii", "HI"),
    "16": ("Idaho", "ID"), "17": ("Illinois", "IL"), "18": ("Indiana", "IN"),
    "19": ("Iowa", "IA"), "20": ("Kansas", "KS"), "21": ("Kentucky", "KY"),
    "22": ("Louisiana", "LA"), "23": ("Maine", "ME"), "24": ("Maryland", "MD"),
    "25": ("Massachusetts", "MA"), "26": ("Michigan", "MI"), "27": ("Minnesota", "MN"),
    "28": ("Mississippi", "MS"), "29": ("Missouri", "MO"), "30": ("Montana", "MT"),
    "31": ("Nebraska", "NE"), "32": ("Nevada", "NV"), "33": ("New Hampshire", "NH"),
    "34": ("New Jersey", "NJ"), "35": ("New Mexico", "NM"), "36": ("New York", "NY"),
    "37": ("North Carolina", "NC"), "38": ("North Dakota", "ND"), "39": ("Ohio", "OH"),
    "40": ("Oklahoma", "OK"), "41": ("Oregon", "OR"), "42": ("Pennsylvania", "PA"),
    "44": ("Rhode Island", "RI"), "45": ("South Carolina", "SC"), "46": ("South Dakota", "SD"),
    "47": ("Tennessee", "TN"), "48": ("Texas", "TX"), "49": ("Utah", "UT"),
    "50": ("Vermont", "VT"), "51": ("Virginia", "VA"), "53": ("Washington", "WA"),
    "54": ("West Virginia", "WV"), "55": ("Wisconsin", "WI"), "56": ("Wyoming", "WY"),
}
ABBREV_TO_FIPS = {ab: f for f, (_, ab) in STATES.items()}

# Facility filter definitions (confirmed 2026-09-05 against the enumerated
# service taxonomy; see scripts/01_fetch_samhsa.py output):
#   OTP  = License/Certification includes the SAMHSA OTP certification string
#          (2,037 nationally; matches SAMHSA's published ~2,100 OTPs)
#   MOUD = uses buprenorphine or methadone in treatment, OR is an OTP
#          (naltrexone-only and "accepts MAT prescribed elsewhere" excluded)
OTP_STRING = "SAMHSA certification for opioid treatment program (OTP)"
MOUD_OM_STRINGS = {"Buprenorphine used in Treatment", "Methadone used in Treatment"}
TERRITORIES = {"PR", "GU", "MP", "VI", "AS"}


def audit_join(source, fips_set):
    """Every source must cover exactly the 51 STATES fips codes."""
    missing = set(STATES) - fips_set
    extra = fips_set - set(STATES)
    if missing or extra:
        sys.exit(f"JOIN AUDIT FAILED for {source}: missing={sorted(missing)} extra={sorted(extra)}")
    print(f"  join audit [{source}]: 51/51 states matched, no extras")


def load_samhsa():
    """Per-state facility counts under the confirmed OTP/MOUD definitions."""
    rows = json.loads((RAW / "samhsa" / "facilities_sa_raw.json").read_text(encoding="utf-8"))

    # address-level dedupe (multiple program listings at one site)
    seen, uniq = set(), []
    for r in rows:
        k = (r["name1"], r.get("street1"), r.get("city"), r.get("state"), r.get("zip"))
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    dropped_t = [r for r in uniq if r.get("state") in TERRITORIES]
    usable = [r for r in uniq if r.get("state") not in TERRITORIES]
    print(f"SAMHSA: {len(rows)} rows -> {len(uniq)} after dedupe -> {len(usable)} after dropping {len(dropped_t)} territory rows")

    def svc(row, code):
        out = set()
        for s in row.get("services", []):
            if s.get("f2") == code:
                out.update(x.strip() for x in s.get("f3", "").split("; ") if x.strip())
        return out

    recs = []
    for r in usable:
        is_otp = OTP_STRING in svc(r, "LCA")
        recs.append({
            "fips": ABBREV_TO_FIPS[r["state"]],
            "otp": is_otp,
            "moud": is_otp or bool(svc(r, "OM") & MOUD_OM_STRINGS),
        })
    df = (pd.DataFrame(recs).groupby("fips").agg(otp_count=("otp", "sum"), moud_count=("moud", "sum"))
          .astype(int).reset_index())
    audit_join("SAMHSA", set(df["fips"]))
    return df


def load_census():
    """ACS 2023 1-year: population, median HH income, insurance coverage."""
    det = json.loads((RAW / "census" / "acs2023_detailed.json").read_text(encoding="utf-8"))
    prof = json.loads((RAW / "census" / "acs2023_profile.json").read_text(encoding="utf-8"))
    d = pd.DataFrame(det[1:], columns=det[0]).rename(columns={
        "B01003_001E": "acs_pop", "B19013_001E": "median_hh_income", "state": "fips"})
    p = pd.DataFrame(prof[1:], columns=prof[0]).rename(columns={
        "DP03_0096PE": "pct_insured", "DP03_0099PE": "pct_uninsured", "state": "fips"})
    df = d.merge(p[["fips", "pct_insured", "pct_uninsured"]], on="fips", validate="1:1")
    df = df[df["fips"].isin(STATES)]  # drops Puerto Rico (72)
    for c in ("acs_pop", "median_hh_income", "pct_insured", "pct_uninsured"):
        df[c] = pd.to_numeric(df[c])
    print(f"Census: {len(df)} state rows kept (PR dropped)")
    audit_join("Census", set(df["fips"]))
    return df[["fips", "acs_pop", "median_hh_income", "pct_insured", "pct_uninsured"]]


def find_cdc_files():
    """All WONDER exports, alphabetical: the 2014-2020 file (database D77) sorts
    before the 2018-2024 file (D157), so keep='last' below prefers the newer
    database for the overlap years 2018-2020."""
    return sorted(list((RAW / "cdc").glob("*.txt")) + list((RAW / "cdc").glob("*.csv")))


def load_wonder(path):
    """Parse a CDC WONDER export: tab-delimited, trailing Notes block.
    Returns long df: fips, state, year, deaths, population, crude_rate, aa_rate."""
    df = pd.read_csv(path, sep="\t", dtype=str)

    # normalize column names -> canonical keys
    def canon(name):
        return " ".join(str(name).lower().replace("-", " ").split())
    cols = {canon(c): c for c in df.columns}
    need = {
        "state": "state", "state code": "fips", "deaths": "deaths", "population": "population",
        "crude rate": "crude_rate",
    }
    year_col = cols.get("year code") or cols.get("year")
    aa_col = next((cols[k] for k in cols if k.startswith("age adjusted rate") and "ci" not in k
                   and "standard error" not in k), None)
    missing = [k for k in need if k not in cols] + ([] if year_col else ["year"])
    if missing:
        sys.exit(f"WONDER parse: expected columns not found: {missing}.\nColumns present: {list(df.columns)}")

    out = pd.DataFrame({
        "state": df[cols["state"]],
        "fips": df[cols["state code"]].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(2),
        "year": df[year_col],
        "deaths": df[cols["deaths"]],
        "population": df[cols["population"]],
        "crude_rate": df[cols["crude rate"]],
        "aa_rate": df[aa_col] if aa_col else pd.NA,
    })
    # notes block rows have no State Code
    out = out[df[cols["state code"]].notna() & (df[cols["state code"]].astype(str).str.strip() != "")]

    # numeric coercion; 'Suppressed'/'Unreliable'/'Not Applicable' -> NaN, reported
    for c in ("deaths", "population", "crude_rate", "aa_rate"):
        raw = out[c].astype(str).str.replace(",", "", regex=False)
        out[c] = pd.to_numeric(raw, errors="coerce")
        bad = out[c].isna().sum()
        if bad:
            print(f"  WONDER: {bad} non-numeric values in '{c}' (suppressed/unreliable) -> NaN")
    out["year"] = out["year"].astype(float).astype(int)
    if aa_col is None:
        print("  WONDER WARNING: no Age Adjusted Rate column found — export was made without it")

    # grain check
    years = sorted(out["year"].unique())
    print(f"WONDER: {len(out)} state-year rows | years {years[0]}-{years[-1]} | states {out['fips'].nunique()}")
    expect = 51 * len(years)
    if len(out) != expect:
        print(f"  WONDER WARNING: expected {expect} rows (51 x {len(years)}), got {len(out)}")
    audit_join("CDC WONDER", set(out["fips"]))
    return out


def main():
    print("=== Merge step ===")
    samhsa = load_samhsa()
    census = load_census()
    capacity = samhsa.merge(census, on="fips", validate="1:1")
    capacity["otp_per_100k"] = capacity["otp_count"] / capacity["acs_pop"] * 1e5
    capacity["moud_per_100k"] = capacity["moud_count"] / capacity["acs_pop"] * 1e5

    cdc_paths = find_cdc_files()
    if not cdc_paths:
        print("\nNO CDC FILE in data/raw/cdc/ — partial mode (nothing written).")
        print("Capacity + ACS table for verification:\n")
        show = capacity.copy()
        show["state"] = show["fips"].map(lambda f: STATES[f][1])
        print(show[["state", "otp_count", "moud_count", "otp_per_100k", "moud_per_100k",
                    "median_hh_income", "pct_uninsured"]]
              .sort_values("moud_per_100k", ascending=False).round(2).to_string(index=False))
        return

    wonder_frames = []
    for p in cdc_paths:
        print(f"\nCDC file: {p.name}")
        wonder_frames.append(load_wonder(p))
    wonder = pd.concat(wonder_frames, ignore_index=True)
    # overlap validation: any state-year present in two exports must agree on deaths
    dup = wonder[wonder.duplicated(["fips", "year"], keep=False)]
    if len(dup):
        per_cell = dup.groupby(["fips", "year"])["deaths"].nunique()
        print(f"CDC overlap: {len(per_cell)} state-year cells in both files; "
              f"{int((per_cell > 1).sum())} disagree on deaths")
    wonder = wonder.drop_duplicates(["fips", "year"], keep="last").reset_index(drop=True)
    print(f"CDC combined: {len(wonder)} state-year rows, years "
          f"{int(wonder.year.min())}-{int(wonder.year.max())}")

    latest = int(wonder["year"].max())
    w_latest = wonder[wonder["year"] == latest].set_index("fips")

    # 3-year pooled rate (sum deaths / sum WONDER population, damps small-state noise)
    last3 = wonder[wonder["year"] >= latest - 2]
    pooled = (last3.groupby("fips").agg(d3=("deaths", "sum"), p3=("population", "sum")))
    pooled["rate_3yr"] = pooled["d3"] / pooled["p3"] * 1e5

    merged = capacity.set_index("fips").join(w_latest[["deaths", "population", "crude_rate", "aa_rate"]])
    merged = merged.join(pooled[["d3", "rate_3yr"]]).reset_index()
    merged.insert(1, "state", merged["fips"].map(lambda f: STATES[f][0]))
    merged.insert(2, "abbrev", merged["fips"].map(lambda f: STATES[f][1]))
    merged.insert(3, "year_latest", latest)
    merged = merged.rename(columns={"deaths": "deaths_latest", "population": "wonder_pop",
                                    "d3": "deaths_3yr"})
    merged["deaths_per_100k_acs"] = merged["deaths_latest"] / merged["acs_pop"] * 1e5
    merged["deaths_per_moud_facility"] = merged["deaths_latest"] / merged["moud_count"]

    OUT.mkdir(parents=True, exist_ok=True)
    merged.round(3).to_csv(OUT / "state_merged.csv", index=False)
    wonder.round(3).to_csv(OUT / "mortality_by_state_year.csv", index=False)
    print(f"\nWrote {OUT / 'state_merged.csv'} ({len(merged)} rows)")
    print(f"Wrote {OUT / 'mortality_by_state_year.csv'} ({len(wonder)} rows)")
    print("\nSpot check (top 5 by age-adjusted mortality rate):")
    print(merged.sort_values("aa_rate", ascending=False)
          [["abbrev", "deaths_latest", "aa_rate", "crude_rate", "rate_3yr",
            "otp_per_100k", "moud_per_100k", "deaths_per_moud_facility"]]
          .head().round(2).to_string(index=False))


if __name__ == "__main__":
    main()
