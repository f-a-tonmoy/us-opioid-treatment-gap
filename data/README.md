# Data

```
data/
├── raw/
│   ├── cdc/          # CDC WONDER mortality exports (manual; recipe below)
│   ├── samhsa/       # facility locator pull (scripts/01)
│   └── census/       # ACS pulls (scripts/02)
└── processed/        # outputs of scripts/03
```

**No data lives in this repo.** This file documents how to obtain and
regenerate all of it: two sources are one script each, and the CDC mortality
data is a short manual export.

---

## 1. CDC WONDER: opioid overdose deaths by state × year (2014–2024)

CDC WONDER has **no public API for sub-national mortality**; queries must go
through the web interface. To export:

Two databases at [wonder.cdc.gov](https://wonder.cdc.gov) (Multiple Cause of
Death): **D157** ("2018–2024, Single Race") and **D77** ("1999–2020", used for
2014–2017). Identical query in both:

- **Group by:** State, then Year. Measures: Deaths, Population, Crude Rate,
  Age-adjusted Rate. Uncheck "Show Totals".
- **Underlying cause (ICD-10):** `X40–X44, X60–X64, X85, Y10–Y14`
  (drug overdose, all intents).
- **Multiple cause (ICD-10):** `T40.0, T40.1, T40.2, T40.3, T40.4, T40.6`
  (any opioid involved).
- Export, and save the `.txt` into `data/raw/cdc/`. **Name files so the newer
  database sorts last alphabetically** (e.g. `..._2014_2020.txt`,
  `..._2018_2024.txt`); the merge keeps the last file's values for overlap
  years.

Sanity anchors (national totals must match published NCHS figures exactly):
80,411 (2021) · 81,806 (2022) · 79,358 (2023) · 54,045 (2024).

Three further WONDER queries feed the corrections, all from `D157`:

- **Drug specificity**, state x year, 2018-2024, two runs on the same UCD
  codes: one with no multiple-cause filter (every drug overdose death) and one
  with MCD restricted to `T36-T50.8` (deaths naming a specific drug). The
  difference is deaths with no drug named. Save the two results together as
  `data/raw/cdc/drug_specificity_capture.csv` with columns
  `fips,year,overdose_total,overdose_specified`; `scripts/03a_drug_specificity.py`
  reads that file.
- **Race**, national, `Single Race 6` x Year, 2018-2024, same opioid
  definition. Save as `data/raw/cdc/opioid_by_race_2018_2024.csv` with columns
  `race,year,deaths,pop,rate`. Do not splice this with `D77`: that database
  uses bridged-race categories, which are not comparable.

Mortality exports must keep the `wonder_mcd_*` filename prefix. The merge step
globs on that prefix so the capture files above are not mistaken for mortality
exports.

## 2. SAMHSA FindTreatment.gov: treatment facilities

```
python scripts/01_fetch_samhsa.py
```

No key needed. Pulls the full national substance-use facility list
(~12,000 rows, ~55 MB JSON) from the locator's `exportsAsJson/v2` endpoint and
prints the service taxonomy used to define the OTP / MOUD filters. Note: this
is a **snapshot**; the locator changes as facilities update their listings,
so a re-pull will differ slightly from the September 2026 snapshot behind the
article.

## 3. SAMHSA certified OTP directory

The authoritative list of certified opioid treatment programs, exported as CSV
by `scripts/01_fetch_samhsa.py` alongside the locator pull. Use this for any
OTP count. The locator's own self-reported certification flag agrees at the
national level (2,037 vs 2,100) but is wrong state by state: threefold too high
in West Virginia, a third too low in Massachusetts and New York, and it claims
a program in Wyoming, which has none.

    https://www.samhsa.gov/find-help/locators/opioid-treatment-program-directory

## 4. U.S. Census: ACS 2023 1-year estimates

**Requires a free API key** (data queries reject keyless requests):

1. Request one at <https://api.census.gov/data/key_signup.html> (arrives by
   email in minutes).
2. Create a `.env` file at the repo root (gitignored; never commit it):

   ```
   CENSUS_API_KEY=your_key_here
   ```

3. ```
   python scripts/02_fetch_census.py
   ```

Pulls `B01003` (population), `B19013` (median household income), and `DP03`
(insurance coverage) for all states, verifying each variable's label against
the API's metadata first.

## 5. Processed outputs

```
python scripts/03_clean_merge.py
```

- **`processed/state_merged.csv`**: one row per state (50 + DC). Key columns:
  `fips`, `state`, `abbrev`, latest-year deaths and crude/age-adjusted rates,
  `deaths_3yr`/`rate_3yr` (2022–2024 combined), `otp_count`, `moud_count`,
  per-100k densities, `deaths_per_moud_facility`, ACS population / median
  income / insurance coverage.
- **`processed/mortality_by_state_year.csv`**: long table, 51 × 11 years
  (2014–2024): deaths, population, crude and age-adjusted rates.

The merge halts unless every source matches a fixed 51-jurisdiction FIPS table
exactly, and it verifies that the two CDC databases agree on all overlapping
state-year cells before splicing them.
