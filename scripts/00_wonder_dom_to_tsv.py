# Convert the CDC WONDER results captured from the browser DOM into the
# canonical WONDER tab-separated export format that 03_clean_merge.py expects.
#
# Provenance: CDC WONDER has no public API for sub-national mortality queries,
# so the query was run through the regular wonder.cdc.gov web interface
# (database D157, "Multiple Cause of Death, 2018-2024, Single Race"), driven
# via browser automation on 2026-09-05 with the user's consent to the data-use
# terms. Query criteria (echoed by WONDER, saved alongside the data):
#   UCD ICD-10: X40-X44, X60-X64, X85, Y10-Y14  (drug overdose, all intents)
#   MCD ICD-10: T40.0-T40.4, T40.6              (opioid involvement)
#   Group by:   State; Year (2018-2024), 50 states + DC, totals off
#   Measures:   Deaths, Population, Crude Rate (95% CI), Age Adjusted Rate (95% CI)
#
# Input : the raw DOM dump JSON (table rows + notes + criteria), kept in
#         data/raw/cdc/wonder_dom_dump.json
# Output: data/raw/cdc/wonder_mcd_opioid_state_year_2018_2024.txt

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# defaults; pass <dump.json> <out.txt> as args for other pulls (e.g. the 2014-2020
# extension from database D77, "Multiple Cause of Death, 1999-2020")
DUMP = Path(sys.argv[1]) if len(sys.argv) > 2 else ROOT / "data" / "raw" / "cdc" / "wonder_dom_dump.json"
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "data" / "raw" / "cdc" / "wonder_mcd_opioid_state_year_2018_2024.txt"


def split_rate(cell):
    """'7.8\\n(7.0 - 8.6)' -> ('7.8', '7.0', '8.6')"""
    m = re.match(r"^([\d.,]+|Unreliable|Suppressed)\s*(?:\(([\d.,]+)\s*-\s*([\d.,]+)\))?", cell.strip())
    return (m.group(1), m.group(2) or "", m.group(3) or "") if m else (cell.strip(), "", "")


def main():
    payload = json.loads(DUMP.read_text(encoding="utf-8"))
    # tool-result wrapper: [{type, text}] with our JSON string inside,
    # possibly followed by non-JSON footer lines -> raw_decode takes the JSON
    if isinstance(payload, list):
        payload, _ = json.JSONDecoder().raw_decode(payload[0]["text"])
    if isinstance(payload, str):  # JS returned a JSON string -> decode once more
        payload = json.loads(payload)
    rows, notes, criteria = payload["rows"], payload["notes"], payload["criteria"]

    state_re = re.compile(r"^(.*)\s\((\d{2})\)$")  # 'Alabama (01)'
    out_rows, cur_state, cur_fips = [], None, None
    for cells in rows:
        cells = [c for c in cells]
        if len(cells) == 6 and state_re.match(cells[0]):
            m = state_re.match(cells[0])
            cur_state, cur_fips = m.group(1), m.group(2)
            data = cells[1:]
        elif len(cells) == 5 and cur_state and re.match(r"^20\d\d$", cells[0].strip()):
            data = cells
        else:
            continue  # header/blank/etc.
        year, deaths, pop = data[0].strip(), data[1].strip(), data[2].strip()
        crude, c_lo, c_hi = split_rate(data[3])
        aar, a_lo, a_hi = split_rate(data[4])
        out_rows.append([cur_state, cur_fips, year, year, deaths, pop,
                         crude, c_lo, c_hi, aar, a_lo, a_hi])

    # integrity: 51 states x 7 years
    states = {r[1] for r in out_rows}
    years = sorted({r[2] for r in out_rows})
    print(f"parsed {len(out_rows)} state-year rows | {len(states)} states | years {years[0]}-{years[-1]}")
    if len(out_rows) != len(states) * len(years):
        sys.exit(f"grain mismatch: {len(states)} x {len(years)} != {len(out_rows)}")

    # national sums per year for eyeball sanity vs published opioid death counts
    for y in years:
        tot = sum(int(r[4].replace(",", "")) for r in out_rows if r[2] == y)
        print(f"  national opioid deaths {y}: {tot:,}")

    header = ["Notes", "State", "State Code", "Year", "Year Code", "Deaths", "Population",
              "Crude Rate", "Crude Rate Lower 95% Confidence Interval",
              "Crude Rate Upper 95% Confidence Interval", "Age Adjusted Rate",
              "Age Adjusted Rate Lower 95% Confidence Interval",
              "Age Adjusted Rate Upper 95% Confidence Interval"]
    q = lambda v: '"' + str(v).replace('"', "'") + '"'
    lines = ["\t".join(q(h) for h in header)]
    for r in out_rows:
        lines.append("\t".join([""] + [q(v) for v in r]))
    lines.append(q("---"))
    for blk in ("criteria", "notes"):
        for ln in payload[blk].splitlines():
            ln = ln.replace("\t", " ").strip()
            if ln:
                lines.append(q(ln))

    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
