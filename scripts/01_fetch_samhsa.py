# Fetch all substance-abuse treatment facilities from the findtreatment.gov
# locator API (SAMHSA) and save the raw JSON.
#
# API notes (verified by probing, 2026-09-05):
#   - Endpoint: https://findtreatment.gov/locator/exportsAsJson/v2
#   - sType=SA        -> substance-abuse facilities (opioid treatment lives here)
#   - sAddr           -> "lat,lon" anchor; only used for the 'miles' sort field
#   - omitting limitType -> NO radius filter, returns the full nationwide set
#   - pageSize up to 2000 works; recordCount ~12,220
#   - Each row's 'services' is a list of {f1: category name, f2: category code,
#     f3: '; '-joined service strings}. Opioid-relevant categories:
#       TC  = Type of Care, OM = Opioid Medications used in Treatment,
#       OT  = Type of Opioid Treatment, PHR = Pharmacotherapies,
#       LCA = License/Certification/Accreditation
#
# This script only ACQUIRES and PROFILES. Filtering to OTP / buprenorphine /
# methadone facilities happens in the merge step, after the filter strings are
# confirmed against the enumeration this script prints.

import json
import time
from collections import Counter
from pathlib import Path

import requests

BASE = "https://findtreatment.gov/locator/exportsAsJson/v2"
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) data-journalism research script"}
OUT = Path(__file__).resolve().parents[1] / "data" / "raw" / "samhsa" / "facilities_sa_raw.json"
PAGE_SIZE = 2000

# Categories whose distinct values we enumerate for the filter decision.
OPIOID_CATEGORIES = ("TC", "OM", "OT", "PHR", "LCA")


def fetch_all():
    """Page through the full nationwide SA facility list."""
    rows, page, expected = [], 1, None
    while True:
        r = requests.get(
            BASE,
            params={"sType": "SA", "sAddr": "39.83,-98.58", "pageSize": PAGE_SIZE, "page": page},
            headers=UA,
            timeout=120,
        )
        r.raise_for_status()
        j = r.json()
        expected = j["recordCount"]
        rows.extend(j["rows"])
        print(f"  page {page}/{j['totalPages']}: +{len(j['rows'])} rows (total {len(rows)})")
        if page >= j["totalPages"]:
            break
        page += 1
        time.sleep(0.5)  # be polite to the API
    return rows, expected


def main():
    print("Fetching nationwide SA facility list from findtreatment.gov ...")
    rows, expected = fetch_all()

    # --- integrity checks ---
    assert len(rows) == expected, f"got {len(rows)} rows, API reported {expected}"
    # No stable facility ID in the payload; dedupe on address identity.
    keys = [(r["name1"], r.get("street1"), r.get("city"), r.get("state"), r.get("zip")) for r in rows]
    dupes = len(keys) - len(set(keys))
    print(f"\nTotal facilities: {len(rows)} (API recordCount {expected}, address-dupes: {dupes})")

    OUT.write_text(json.dumps(rows), encoding="utf-8")
    print(f"Raw JSON saved -> {OUT} ({OUT.stat().st_size / 1e6:.1f} MB)")

    # --- profile: states ---
    states = Counter(r.get("state") for r in rows)
    print(f"\nDistinct 'state' values ({len(states)}):")
    print("  " + ", ".join(f"{s}:{n}" for s, n in sorted(states.items(), key=lambda x: str(x[0]))))

    # --- profile: distinct service strings in opioid-relevant categories ---
    print("\nDistinct service strings by category (for filter definition):")
    for code in OPIOID_CATEGORIES:
        vals = Counter()
        for r in rows:
            for svc in r.get("services", []):
                if svc.get("f2") == code:
                    for item in svc.get("f3", "").split("; "):
                        if item.strip():
                            vals[item.strip()] += 1
        name = next(
            (s["f1"] for r in rows for s in r.get("services", []) if s.get("f2") == code),
            code,
        )
        print(f"\n[{code}] {name.strip()}:")
        for val, n in vals.most_common():
            print(f"  {n:6d}  {val}")


if __name__ == "__main__":
    main()
