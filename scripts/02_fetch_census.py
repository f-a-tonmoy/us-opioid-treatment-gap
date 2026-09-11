# Fetch state-level ACS 2023 1-year estimates from the Census API:
#   - B01003_001E   total population            (detailed tables endpoint)
#   - B19013_001E   median household income     (detailed tables endpoint)
#   - DP03_0096PE   % with health insurance     (data profile endpoint)
#   - DP03_0099PE   % without health insurance  (data profile endpoint)
#
# Vintage 2023 matches the latest final CDC mortality year. Data queries need
# a free API key; only the metadata endpoints stay keyless.
#
# Each variable's label is fetched from the API's own variables.json first and
# printed, so the meaning is verified from the source rather than assumed.

import json
import os
import sys
from pathlib import Path

import requests

YEAR = 2023
DETAILED = f"https://api.census.gov/data/{YEAR}/acs/acs1"
PROFILE = f"https://api.census.gov/data/{YEAR}/acs/acs1/profile"
OUTDIR = Path(__file__).resolve().parents[1] / "data" / "raw" / "census"

VARS = {  # variable -> (endpoint, sanity word that must appear in its label/concept)
    "B01003_001E": (DETAILED, "population"),
    "B19013_001E": (DETAILED, "income"),
    "DP03_0096PE": (PROFILE, "insurance"),
    "DP03_0099PE": (PROFILE, "insurance"),
}


def api_key():
    """The Census API requires a key for data queries (free, instant:
    https://api.census.gov/data/key_signup.html). Metadata endpoints stay
    keyless. Reads CENSUS_API_KEY from the environment, else from the
    gitignored .env file at the project root."""
    key = os.environ.get("CENSUS_API_KEY", "").strip()
    if not key:
        envfile = Path(__file__).resolve().parents[1] / ".env"
        if envfile.exists():
            for line in envfile.read_text(encoding="utf-8").splitlines():
                if line.strip().startswith("CENSUS_API_KEY="):
                    key = line.split("=", 1)[1].strip().strip("'\"")
    if not key:
        sys.exit("No Census API key. Set CENSUS_API_KEY, or put CENSUS_API_KEY=<key> in .env")
    return key


def verify_labels():
    """Pull each variable's official definition and sanity-check it."""
    print("Variable definitions (from the Census API's variables.json):")
    for var, (endpoint, expect) in VARS.items():
        meta = requests.get(f"{endpoint}/variables/{var}.json", timeout=60).json()
        label, concept = meta.get("label", ""), meta.get("concept", "")
        print(f"  {var}: {label}  [concept: {concept}]")
        assert expect.lower() in (label + concept).lower(), (
            f"{var} label does not mention '{expect}' — wrong variable?"
        )


def fetch(endpoint, varlist, fname, key):
    """Fetch NAME + varlist for all states, save raw JSON, return rows."""
    url = f"{endpoint}?get=NAME,{','.join(varlist)}&for=state:*&key={key}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    if "json" not in r.headers.get("content-type", ""):
        sys.exit(f"Census API returned non-JSON (bad/inactive key?): {r.text[:200]}")
    data = r.json()
    (OUTDIR / fname).write_text(json.dumps(data), encoding="utf-8")
    print(f"Saved {fname}: {len(data) - 1} state rows")
    return data


def main():
    key = api_key()
    verify_labels()
    print()
    detailed = fetch(DETAILED, ["B01003_001E", "B19013_001E"], f"acs{YEAR}_detailed.json", key)
    profile = fetch(PROFILE, ["DP03_0096PE", "DP03_0099PE"], f"acs{YEAR}_profile.json", key)

    # Display merged-by-FIPS table (real merge happens in 03_clean_merge.py).
    det = {row[-1]: row for row in detailed[1:]}   # last column = state FIPS
    prof = {row[-1]: row for row in profile[1:]}
    assert set(det) == set(prof), "FIPS mismatch between the two Census pulls"

    print(f"\n{'FIPS':>4} {'State':<22} {'Population':>11} {'MedHHInc':>9} {'%Insured':>8} {'%Uninsur':>8}")
    for fips in sorted(det):
        name, pop, inc = det[fips][0], det[fips][1], det[fips][2]
        ins, unins = prof[fips][1], prof[fips][2]
        print(f"{fips:>4} {name:<22} {int(pop):>11,} {int(inc):>9,} {float(ins):>8.1f} {float(unins):>8.1f}")


if __name__ == "__main__":
    main()
