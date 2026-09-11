# Correct opioid death counts for incomplete drug reporting on death
# certificates, and write the result into the processed data.
#
# Roughly one overdose death certificate in twenty names no drug at all, and
# the share varies enormously by state because it depends on local
# medical-examiner and coroner practice. States that report poorly therefore
# look healthier than they are, which matters here because the whole article
# ranks states against each other.
#
# The fix is the standard proportional redistribution used in the literature
# (Ruhm 2018): assume the unspecified deaths in a state contain opioids at the
# same rate as that state's deaths where a drug WAS named.
#
#     adjusted_opioid = opioid x (total_overdose / specified_drug)
#
# It is a simplification of Ruhm's published method, which also models county
# characteristics; the direction and rough size of the correction are the
# same, and it is transparent enough to state in one line.
#
# Inputs : data/raw/cdc/drug_specificity_capture.csv (two CDC WONDER queries)
#          data/processed/mortality_by_state_year.csv (opioid deaths)
# Outputs: data/raw/cdc/drug_specificity_2018_2024.csv  (audit trail)
#          adds unspecified_share / opioid_adj / rate_adj columns to the
#          processed mortality table

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CAPTURE = ROOT / "data" / "raw" / "cdc" / "drug_specificity_capture.csv"
RAW_OUT = ROOT / "data" / "raw" / "cdc" / "drug_specificity_2018_2024.csv"
LONG = ROOT / "data" / "processed" / "mortality_by_state_year.csv"


def main():
    if not CAPTURE.exists():
        sys.exit(f"missing {CAPTURE}\nRun the two WONDER queries in data/README.md "
                 f"(section 1) and save them there as fips,year,overdose_total,"
                 f"overdose_specified.")
    ds = pd.read_csv(CAPTURE, dtype={"fips": str})

    # integrity: a state can never name a drug in more deaths than it recorded
    bad = ds[ds.overdose_specified > ds.overdose_total]
    if len(bad):
        sys.exit(f"specified exceeds total in {len(bad)} state-years")
    ds["unspecified"] = ds.overdose_total - ds.overdose_specified
    ds["unspecified_share"] = ds.unspecified / ds.overdose_total
    ds["adj_factor"] = ds.overdose_total / ds.overdose_specified
    RAW_OUT.parent.mkdir(parents=True, exist_ok=True)
    ds.round(6).to_csv(RAW_OUT, index=False)
    print(f"wrote {RAW_OUT} ({len(ds)} state-years, 2018-2024)")

    long = pd.read_csv(LONG, dtype={"fips": str})
    # drop any columns from an earlier run, so re-running this step re-derives
    # them instead of colliding with itself and producing _x/_y suffixes
    added = ["overdose_total", "overdose_specified", "unspecified_share",
             "adj_factor", "opioid_adj", "rate_adj"]
    long = long.drop(columns=[c for c in added if c in long.columns])
    m = long.merge(ds[["fips", "year", "overdose_total", "overdose_specified",
                       "unspecified_share", "adj_factor"]],
                   on=["fips", "year"], how="left")
    m["opioid_adj"] = (m.deaths * m.adj_factor).round(1)
    m["rate_adj"] = m.opioid_adj / m.population * 1e5
    m.round(4).to_csv(LONG, index=False)
    print(f"updated {LONG} with adjusted columns "
          f"({m.adj_factor.notna().sum()} of {len(m)} rows covered, 2018-2024)")

    # --- what the correction does ---
    latest = ds[ds.year == 2024]
    nat_share = latest.unspecified.sum() / latest.overdose_total.sum()
    print(f"\nnational unspecified share, 2024: {nat_share*100:.1f}%")
    print(f"national opioid deaths 2024: reported "
          f"{int(m.loc[m.year == 2024, 'deaths'].sum()):,} -> adjusted "
          f"{int(m.loc[m.year == 2024, 'opioid_adj'].sum()):,}")

    p = latest.set_index("fips").unspecified_share.mul(100)
    names = pd.read_csv(ROOT / "data" / "processed" / "state_merged.csv",
                        dtype={"fips": str}).set_index("fips")
    p = p.to_frame("pct").join(names[["abbrev", "state"]])
    print("\nworst drug reporting, 2024 (share of overdose deaths naming no drug):")
    print(p.nlargest(8, "pct")[["abbrev", "pct"]].round(1).to_string(index=False))
    print("\nbest reporting:")
    print(p.nsmallest(5, "pct")[["abbrev", "pct"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
