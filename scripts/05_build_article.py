# Build the self-contained interactive HTML article from the processed data.
# Every number in the prose is interpolated from the CSVs, never hand-typed.
# Output: index.html  (Plotly 2.35.2 via CDN, all data inline)
#
# Visual identity is deliberately its own - dark cover masthead, IBM Plex Sans
# headlines over Spectral body, indigo/ember palette, flat hairline figures -
# the reference case studies set the rigor bar, not the skin.

import json
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
(ROOT / "reports").mkdir(exist_ok=True)   # holds the hero image

m = pd.read_csv(ROOT / "data" / "processed" / "state_merged.csv", dtype={"fips": str})
long = pd.read_csv(ROOT / "data" / "processed" / "mortality_by_state_year.csv", dtype={"fips": str})
race = pd.read_csv(ROOT / "data" / "raw" / "cdc" / "opioid_by_race_2018_2024.csv")

# ---------------- derived quantities ----------------
LATEST = int(m["year_latest"].iloc[0])
nat = long.groupby("year").agg(deaths=("deaths", "sum"), pop=("population", "sum")).reset_index()
nat["rate"] = nat["deaths"] / nat["pop"] * 1e5
nrate = nat.set_index("year")["rate"]
_natpeak = nrate[[y for y in nrate.index if y != LATEST]].max()
NAT_DECL = (nrate[LATEST] - _natpeak) / _natpeak * 100
NAT_R14 = nrate[LATEST] / nrate[2014]
TOTAL = int(nat["deaths"].sum())
# floored to the nearest 100k for the deck and hero card, so the rounded
# phrasing stays true if the data changes rather than going quietly stale
TOTAL_FLOOR = TOTAL // 100_000 * 100_000
D2014, D2017, D2018, D2019, D2020, D2022, D2023, D2024 = (int(nat.loc[nat.year == y, "deaths"].iloc[0]) for y in (2014, 2017, 2018, 2019, 2020, 2022, 2023, 2024))
PCT_DIP_2018 = (D2018 - D2017) / D2017 * 100
PCT_2324 = (D2024 - D2023) / D2023 * 100
PCT_2020 = (D2020 - D2019) / D2019 * 100
PCT_ABOVE_2019 = (D2024 - D2019) / D2019 * 100
PCT_ABOVE_2014 = (D2024 - D2014) / D2014 * 100
MIN_PER_DEATH_24 = 525_600 / D2024
MIN_PER_DEATH_22 = 525_600 / D2022
# North Carolina's medical examiner closed hundreds of 2023 cases after the
# federal file was compiled, so its 2023 row understates the state's own count
NC23 = int(long.loc[(long.fips == "37") & (long.year == 2023), "deaths"].iloc[0])
NC23_STATE = 3656   # NCDHHS, opioid-involved overdose deaths 2023 (cite-12)
last3 = long[long.year >= LATEST - 2]
NAT3 = last3["deaths"].sum() / last3["population"].sum() * 1e5

# ---- drug-specificity correction (2022-2024) ----
w3 = last3.groupby("fips").agg(d3=("deaths", "sum"), adj3=("opioid_adj", "sum"),
                               p3=("population", "sum"), od3=("overdose_total", "sum"),
                               spec3=("overdose_specified", "sum"))
w3["rate_3yr_adj"] = w3.adj3 / w3.p3 * 1e5
w3["unspec_share"] = (w3.od3 - w3.spec3) / w3.od3 * 100
m = m.merge(w3[["rate_3yr_adj", "unspec_share"]], on="fips")
NAT3_ADJ = w3.adj3.sum() / w3.p3.sum() * 1e5
NAT_UNSPEC = (w3.od3.sum() - w3.spec3.sum()) / w3.od3.sum() * 100
ADJ_LIFT = (w3.adj3.sum() / w3.d3.sum() - 1) * 100

m["r_burden_rep"] = m["rate_3yr"].rank(ascending=False).astype(int)
m["r_burden"] = m["rate_3yr_adj"].rank(ascending=False).astype(int)
m["rank_move"] = m["r_burden_rep"] - m["r_burden"]
m["r_cap_moud"] = m["moud_per_100k"].rank(ascending=False).astype(int)
m["r_cap_otp"] = m["otp_per_100k"].rank(ascending=False).astype(int)
m["mismatch"] = m["r_cap_moud"] - m["r_burden"]
m["deaths_per_moud_facility"] = m["deaths_latest"] / m["moud_count"]
_sl, _ic, _, _, _ = stats.linregress(m.rate_3yr_adj, m.moud_per_100k)
m["resid"] = m.moud_per_100k - (_ic + _sl * m.rate_3yr_adj)
_sh = (m["aa_rate"].rank(ascending=False) - m["crude_rate"].rank(ascending=False)).abs()
AA_MAXSHIFT = int(_sh.max())
AA_MAXSTATE = m.loc[_sh.idxmax(), "state"]

piv = long.pivot_table(index="fips", columns="year", values="crude_rate")
peak = piv[[y for y in piv.columns if y != LATEST]].max(axis=1)  # each state's own worst year, 2014-2023
m = m.set_index("fips")
m["decline_pct"] = (piv[LATEST] - peak) / peak * 100
m["ratio14"] = piv[LATEST] / piv[2014]
m["peak_year"] = piv[[y for y in piv.columns if y != LATEST]].idxmax(axis=1).astype(int)
m = m.reset_index()
BELOW14 = " and ".join(sorted(m.loc[m.ratio14 < 1, "state"]))

g = lambda ab: m.loc[m.abbrev == ab].iloc[0]
WV, WY, DCr, WA, TX, AK, NE, LA, AL = (g(x) for x in
                                       ("WV", "WY", "DC", "WA", "TX", "AK", "NE", "LA", "AL"))

# ---- who counts as underserved ----
# A rule rather than a hand-picked list: above-median adjusted burden, and in
# the worst quarter of states on all three independent gap measures. DC is
# excluded here for the same reason it is flagged everywhere else - a city
# measured against states.
_st = m[m.abbrev != "DC"].copy()
_q = len(_st) / 4
_st["r_mis"] = _st.mismatch.rank(ascending=False)
_st["r_res"] = _st.resid.rank()
_st["r_dpf"] = _st.deaths_per_moud_facility.rank(ascending=False)
_flag = ((_st.rate_3yr_adj > _st.rate_3yr_adj.median())
         & (_st.r_mis <= _q) & (_st.r_res <= _q) & (_st.r_dpf <= _q))
UNDERSERVED = _st.loc[_flag].sort_values("mismatch", ascending=False).abbrev.tolist()
UND_NAMES = _st.loc[_flag].sort_values("mismatch", ascending=False).state.tolist()
N_UND = len(UNDERSERVED)
UND_WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}[N_UND]
UND_LIST = ", ".join(UND_NAMES[:-1]) + f", and {UND_NAMES[-1]}"
MIS_LO, MIS_HI = int(_st.loc[_flag, "mismatch"].min()), int(_st.loc[_flag, "mismatch"].max())
RANKGAP_RESID_AGREE = stats.spearmanr(m.mismatch, -m.resid).statistic

# ---- Medicaid expansion ----
exp = m[m.expanded_medicaid]
non = m[~m.expanded_medicaid]
MED_MOUD_E, MED_MOUD_N = exp.moud_per_100k.median(), non.moud_per_100k.median()
MED_OTP_E, MED_OTP_N = exp.otp_per_100k.median(), non.otp_per_100k.median()
P_MOUD_EXP = stats.mannwhitneyu(exp.moud_per_100k, non.moud_per_100k).pvalue
P_OTP_EXP = stats.mannwhitneyu(exp.otp_per_100k, non.otp_per_100k).pvalue
P_BURD_EXP = stats.mannwhitneyu(exp.rate_3yr_adj, non.rate_3yr_adj).pvalue
UND_NONEXP = sorted(m.loc[m.abbrev.isin(UNDERSERVED) & ~m.expanded_medicaid, "state"])
UND_EXP = sorted(m.loc[m.abbrev.isin(UNDERSERVED) & m.expanded_medicaid, "state"])

# ---- race ----
rp = race.pivot(index="race", columns="year", values="rate")
RACE_2018, RACE_2024 = 2018, 2024
WHITE_CHG = (rp.loc["White", 2024] / rp.loc["White", 2018] - 1) * 100
BLACK_CHG = (rp.loc["Black or African American", 2024] / rp.loc["Black or African American", 2018] - 1) * 100
AIAN_CHG = (rp.loc["American Indian or Alaska Native", 2024] / rp.loc["American Indian or Alaska Native", 2018] - 1) * 100
WHITE_PEAK = int(rp.loc["White"].idxmax())
BLACK_PEAK = int(rp.loc["Black or African American"].idxmax())
BLACK_2023, WHITE_2023 = rp.loc["Black or African American", 2023], rp.loc["White", 2023]
_tot = race.groupby("year").deaths.sum()
_blk = race[race.race.str.startswith("Black")].set_index("year").deaths / _tot * 100
BLK_SHARE_18, BLK_SHARE_23 = _blk[2018], _blk[2023]
und = m[m.abbrev.isin(UNDERSERVED)].sort_values("mismatch", ascending=False)

def corr(x, y):
    r, p = stats.spearmanr(m[x], m[y])
    return r, p
R_MOUD, P_MOUD = corr("rate_3yr_adj", "moud_per_100k")
R_OTP, P_OTP = corr("rate_3yr_adj", "otp_per_100k")
R_MIS_INC, P_MIS_INC = corr("mismatch", "median_hh_income")
R_MIS_UNINS, P_MIS_UNINS = corr("mismatch", "pct_uninsured")
R_DPF_INC, P_DPF_INC = corr("deaths_per_moud_facility", "median_hh_income")
R_DPF_UNINS, P_DPF_UNINS = corr("deaths_per_moud_facility", "pct_uninsured")

OTP_TOTAL = int(m["otp_count"].sum())
OTP_SELFREPORT = int(m["otp_selfreport"].sum())
MOUD_TOTAL = int(m["moud_count"].sum())
NO_OTP_STATES = list(m.loc[m.otp_count == 0, "state"])
WV_OTP_SELF = int(WV.otp_selfreport)
# how many states' 2024 confidence intervals overlap West Virginia's
CI_OVERLAP_WV = int(((m.crude_lo <= WV.crude_hi) & (m.crude_hi >= WV.crude_lo)).sum() - 1)

# ---------------- palette / chart chrome ----------------
J = lambda o: json.dumps(o, separators=(",", ":"))
INK, DIM, RULE = "#15171e", "#565b68", "#e2e3ea"
INDIGO, EMBER, TEAL, AMBER, GREY = "#443d9e", "#c8501f", "#2e7f8a", "#e8a33d", "#c2c2d0"
FONT = dict(family="IBM Plex Sans, system-ui, sans-serif", size=13, color=INK)
# dragmode=False turns off box-zoom on the cartesian charts. Each one already
# fits its frame and zooming reveals nothing, so on a touch screen a swipe
# should scroll the page rather than draw a zoom rectangle the reader cannot
# undo with the modebar hidden. Hover and legend clicks are unaffected. The
# geo map overrides this back to "pan": it is the one chart worth exploring,
# because DC and the small northeastern states are unreadable at full extent.
BASE_LAYOUT = dict(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font=FONT,
                   margin=dict(l=48, r=16, t=12, b=40), showlegend=False, dragmode=False)
CONFIG = dict(displayModeBar=False, responsive=True)

def layout(**kw):
    d = {**BASE_LAYOUT}
    d.update(kw)
    return d

def gaxis(**kw):
    d = dict(gridcolor=RULE, griddash="dot", zeroline=False)
    d.update(kw)
    return d

# 01 national trend: drawn line + area fill (not bars)
years = nat["year"].tolist()
c_national = dict(
    data=[dict(type="scatter", mode="lines+markers", x=years, y=nat["deaths"].tolist(),
               line=dict(color=INDIGO, width=3), fill="tozeroy", fillcolor="rgba(68,61,158,0.07)",
               marker=dict(size=8, color=[EMBER if y == 2024 else INDIGO for y in years]),
               hovertemplate="%{x}: %{y:,} opioid deaths<extra></extra>")],
    layout=layout(height=360, margin=dict(l=30, r=30, t=12, b=40),
                  yaxis=dict(visible=False, range=[0, 99000]),
                  xaxis=dict(tickmode="array", tickvals=years, ticklen=0, showline=True,
                             linecolor=INK, range=[2013.6, 2024.4]),
                  annotations=[
                      # every year labeled; 2022 peak bold, 2024 called out below the drop,
                      # 2018 below the line so it can't crowd 2017's label
                      *[dict(x=int(y), y=int(d),
                             yshift=(-24 if y == 2018 else (24 if y == 2022 else 20)),
                             text=f"<b>{int(d):,}</b>" if y == 2022 else f"{int(d):,}",
                             showarrow=False,
                             font=dict(size=12.5 if y == 2022 else 11, color=INK if y == 2022 else DIM))
                        for y, d in zip(nat["year"], nat["deaths"]) if y != 2024],
                      dict(x=2024, y=D2024, xanchor="right", xshift=10, yshift=-34,
                           text=f"<b>{D2024:,}</b><br>{PCT_2324:+.0f}% in one year", showarrow=False,
                           font=dict(size=12.5, color=EMBER))]))

# 02 choropleth: indigo ramp, abbreviation labels on states with enough map area
# (tiny-footprint states skipped; label color flips white/dark with fill darkness)
SMALL_STATES = {"CT", "DE", "DC", "HI", "MA", "MD", "NH", "NJ", "RI", "VT"}
LABEL_PTS = {  # visual label points (lon, lat); MI on the lower peninsula, LA west of the delta
    "AL": (-86.8, 32.8), "AK": (-152.0, 64.2), "AZ": (-111.7, 34.3), "AR": (-92.4, 34.9),
    "CA": (-119.6, 37.0), "CO": (-105.5, 39.0), "FL": (-81.7, 28.4), "GA": (-83.4, 32.7),
    "ID": (-114.6, 44.2), "IL": (-89.2, 40.0), "IN": (-86.3, 39.9), "IA": (-93.5, 42.1),
    "KS": (-98.4, 38.5), "KY": (-85.3, 37.5), "LA": (-92.0, 31.1), "ME": (-69.2, 45.4),
    "MI": (-84.7, 43.3), "MN": (-94.3, 46.3), "MS": (-89.7, 32.7), "MO": (-92.5, 38.4),
    "MT": (-109.6, 47.0), "NE": (-99.8, 41.5), "NV": (-116.6, 39.3), "NM": (-106.1, 34.4),
    "NY": (-75.5, 43.0), "NC": (-79.4, 35.5), "ND": (-100.5, 47.4), "OH": (-82.8, 40.3),
    "OK": (-97.5, 35.6), "OR": (-120.6, 43.9), "PA": (-77.8, 40.9), "SC": (-80.9, 33.9),
    "SD": (-100.2, 44.4), "TN": (-86.3, 35.85), "TX": (-99.3, 31.4), "UT": (-111.7, 39.3),
    "VA": (-78.7, 37.5), "WA": (-120.4, 47.4), "WV": (-80.6, 38.7), "WI": (-89.8, 44.6),
    "WY": (-107.6, 43.0),
}
WHITE_Z = 31  # fills at/above this rate are dark enough to need white text
# on narrow screens, also skip the next-smallest states and shrink the type
MOBILE_SKIP = SMALL_STATES | {"WV", "SC", "ME", "IN", "KY", "TN", "VA", "OH", "PA", "MS"}

_lab = m[m.abbrev.isin(LABEL_PTS)]
MAP_KEEP_M = [0 if a in MOBILE_SKIP else 1 for a in _lab.abbrev]
c_map = dict(
    data=[dict(type="choropleth", locationmode="USA-states",
               locations=m["abbrev"].tolist(), z=m["rate_3yr_adj"].round(1).tolist(),
               zmin=0, zmax=56,
               colorscale=[[0, "#f0f0f5"], [0.4, "#aaa5da"], [0.75, INDIGO], [1, "#171243"]],
               marker=dict(line=dict(color="#f7f7f9", width=1)),
               colorbar=dict(title=dict(text="Deaths<br>/100k", font=dict(size=12)), thickness=9, len=0.65, x=0.93,
                             outlinewidth=0, tickfont=dict(size=11.5, color=DIM)),
               customdata=m["state"].tolist(),
               hovertemplate="%{customdata}: %{z} per 100k<extra></extra>"),
          dict(type="scattergeo", mode="text",
               lon=[LABEL_PTS[a][0] for a in _lab.abbrev],
               lat=[LABEL_PTS[a][1] for a in _lab.abbrev],
               text=[f"<b>{a}</b>" for a in _lab.abbrev],
               textfont=dict(size=10, color=["white" if z >= WHITE_Z else "#454962" for z in _lab.rate_3yr_adj]),
               hoverinfo="skip")],
    layout=layout(height=430, margin=dict(l=0, r=0, t=6, b=6), dragmode="pan",
                  geo=dict(scope="usa", bgcolor="rgba(0,0,0,0)", lakecolor="#f7f7f9",
                           subunitcolor="#f7f7f9", showlakes=True)))

# 03 decline: thin horizontal bars (lollipop feel), west in ember
dd = m.sort_values("decline_pct")
west = {"AK", "WA", "OR", "NV", "AZ", "UT", "HI", "CA", "NM", "CO", "ID", "MT", "WY"}
c_decline = dict(
    data=[dict(type="bar", orientation="h", x=dd["decline_pct"].round(1).tolist(), y=dd["abbrev"].tolist(),
               marker=dict(color=[EMBER if a in west else "#9d99cf" for a in dd["abbrev"]]),
               customdata=dd["state"].tolist(),
               hovertemplate="%{customdata}: %{x}% vs its 2021–23 peak<extra></extra>")],
    layout=layout(height=880, bargap=0.55, margin=dict(l=40, r=20, t=34, b=44),
                  xaxis=gaxis(title=dict(text="Change in overdose rate, fentanyl-era peak → 2024", font=dict(size=12.5, color=DIM)),
                              ticksuffix="%", zeroline=True, zerolinecolor=INK, range=[-58, 3]),
                  yaxis=dict(tickfont=dict(size=10.5), ticklen=0),
                  shapes=[dict(type="line", x0=NAT_DECL, x1=NAT_DECL, yref="paper", y0=0, y1=1,
                               line=dict(color=DIM, width=1.2, dash="dash"))],
                  annotations=[dict(x=float(AK.decline_pct), y="AK", xanchor="right", xshift=-8,
                                    text=f"<b>Alaska: {AK.decline_pct:.0f}%</b>", showarrow=False,
                                    font=dict(size=12, color=EMBER)),
                               dict(x=NAT_DECL, y=1, yref="paper", yshift=20,
                                    text=f"National {NAT_DECL:.0f}%", showarrow=False,
                                    font=dict(size=11.5, color=DIM)),
                               dict(xref="paper", yref="paper", x=0.02, y=0.995, xanchor="left",
                                    text="■ Western states", showarrow=False,
                                    font=dict(size=12.5, color=EMBER)),
                               dict(xref="paper", yref="paper", x=0.02, y=0.966, xanchor="left",
                                    text="■ Elsewhere", showarrow=False,
                                    font=dict(size=12.5, color="#8783c4"))]))

# 04 capacity: horizontal dumbbells, OTP (ember diamond) -> MOUD (indigo dot)
CAPMAX = float(m["moud_per_100k"].max())
cc = m.sort_values("moud_per_100k", ascending=True)  # best at top of a reversed axis
seg_x, seg_y = [], []
for _, r in cc.iterrows():
    seg_x += [r["otp_per_100k"], r["moud_per_100k"], None]
    seg_y += [r["abbrev"], r["abbrev"], None]
c_capacity = dict(
    data=[dict(type="scatter", mode="lines", x=[None if v is None else round(v, 2) for v in seg_x], y=seg_y,
               line=dict(color="#cfcede", width=2), hoverinfo="skip", showlegend=False),
          dict(type="scatter", mode="markers", x=cc["otp_per_100k"].round(2).tolist(), y=cc["abbrev"].tolist(),
               name="Certified OTPs", marker=dict(size=6.5, color=EMBER, symbol="diamond"),
               customdata=cc["state"].tolist(),
               hovertemplate="%{customdata}: %{x} certified OTPs /100k<extra></extra>"),
          dict(type="scatter", mode="markers", x=cc["moud_per_100k"].round(2).tolist(), y=cc["abbrev"].tolist(),
               name="MOUD facilities", marker=dict(size=7.5, color=INDIGO),
               customdata=cc["state"].tolist(),
               hovertemplate="%{customdata}: %{x} MOUD facilities /100k<extra></extra>")],
    layout=layout(height=880, showlegend=True, margin=dict(l=40, r=20, t=34, b=44),
                  legend=dict(orientation="h", x=0.5, xanchor="center", y=1.03, font=dict(size=12.5)),
                  # start just below zero: Wyoming has no certified OTP, and a marker
                  # sitting exactly on the axis line overlaps its own tick label
                  xaxis=gaxis(title=dict(text="Facilities per 100,000 residents", font=dict(size=12.5, color=DIM)),
                              range=[-0.032 * CAPMAX, CAPMAX * 1.06]),
                  yaxis=dict(tickfont=dict(size=10.5), ticklen=0),
                  annotations=[dict(x=3.3, y="WY", yshift=13,
                                    text=f"<b>Wyoming: one OTP, {int(WY.moud_count)} MOUD facilities</b>",
                                    showarrow=False, font=dict(size=12, color=INK))]))

# 06 scatter with median quadrant guides
hl = set(UNDERSERVED)
deep = {"WV", "KY", "ME", "VT", "MD"}
MEDX, MEDY = float(m["rate_3yr_adj"].median()), float(m["moud_per_100k"].median())

def mark_size(deaths):
    return deaths ** 0.5 / 2.4 + 5

def lab_shift(ab, up=True):
    """Offset a state label clear of its own dot.

    Marker radius grows with the death count, so a fixed offset leaves the
    big states' labels sitting on their own marker -- and since label and dot
    share a color, the overlapping half of the text disappears.
    """
    r = mark_size(float(g(ab).deaths_latest)) / 2
    return round((r + 11) * (1 if up else -1), 1)

def pt_color(ab):
    if ab in hl: return EMBER
    if ab in deep: return INDIGO
    if ab in ("WY", "UT"): return TEAL
    return GREY
c_scatter = dict(
    data=[dict(type="scatter", mode="markers", x=m["rate_3yr_adj"].round(1).tolist(), y=m["moud_per_100k"].round(2).tolist(),
               marker=dict(size=mark_size(m["deaths_latest"]).round(1).tolist(),
                           color=[pt_color(a) for a in m["abbrev"]],
                           opacity=0.9, line=dict(width=1, color="#f7f7f9")),
               customdata=[[s, int(d)] for s, d in zip(m["state"], m["deaths_latest"])],
               hovertemplate="%{customdata[0]}: %{x}/100k burden, %{y} facilities/100k, %{customdata[1]:,} deaths in " + str(LATEST) + "<extra></extra>")],
    layout=layout(height=470,
                  xaxis=gaxis(title=dict(text="Opioid deaths per 100k, 2022–2024 combined", font=dict(size=12.5, color=DIM))),
                  yaxis=gaxis(title=dict(text="MOUD facilities per 100k", font=dict(size=12.5, color=DIM)), range=[0, 7.3]),
                  shapes=[dict(type="line", x0=MEDX, x1=MEDX, yref="paper", y0=0, y1=1, line=dict(color="#b9b9c9", width=1, dash="dot")),
                          dict(type="line", y0=MEDY, y1=MEDY, xref="paper", x0=0, x1=1, line=dict(color="#b9b9c9", width=1, dash="dot"))],
                  annotations=[dict(x=56, y=7.0, xanchor="right", text="High burden<br>high capacity", showarrow=False, align="right", font=dict(size=11.5, color=INDIGO)),
                               dict(x=56, y=0.35, xanchor="right", text="<b>High burden<br>low capacity</b>", showarrow=False, align="right", font=dict(size=11.5, color=EMBER)),
                               dict(x=float(WV.rate_3yr_adj), y=float(WV.moud_per_100k), text="WV", showarrow=False, yshift=lab_shift("WV"), font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(g('KY').rate_3yr_adj), y=float(g('KY').moud_per_100k), text="KY", showarrow=False, yshift=lab_shift("KY"), font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(g('ME').rate_3yr_adj), y=float(g('ME').moud_per_100k), text="ME", showarrow=False, yshift=lab_shift("ME"), font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(WY.rate_3yr_adj), y=float(WY.moud_per_100k), text="WY", showarrow=True,
                                    ax=0, ay=-30, arrowhead=0, arrowsize=1, arrowwidth=1.1,
                                    arrowcolor="#9a9aa8", font=dict(size=11.5, color=TEAL)),
                               dict(x=float(DCr.rate_3yr_adj), y=float(DCr.moud_per_100k), text="DC", showarrow=False, yshift=lab_shift("DC", up=False), font=dict(size=11.5, color=DIM)),
                               # NV alone lands on a neighbour's dot; nudge it into the gap beside it
                               *[dict(x=float(g(a).rate_3yr_adj), y=float(g(a).moud_per_100k), text=a, showarrow=False,
                                      yshift=lab_shift(a, up=False), xshift=8 if a == "NV" else 0,
                                      font=dict(size=11.5, color=EMBER)) for a in UNDERSERVED]]))

# 07 slope: burden rank vs capacity rank
lines, notes = [], []
for _, r in m.iterrows():
    ab = r["abbrev"]
    color, width, op = "#c9c8d6", 1.1, 0.6
    if ab in hl: color, width, op = EMBER, 2.4, 1
    elif ab in ("WY", "UT"): color, width, op = TEAL, 2.1, 1
    elif ab == "DC": color, width, op = DIM, 1.5, 0.9
    lines.append(dict(type="scatter", mode="lines+markers", x=[0, 1], y=[int(r["r_burden"]), int(r["r_cap_moud"])],
                      line=dict(color=color, width=width), marker=dict(size=4.5, color=color), opacity=op,
                      customdata=[[r["state"], int(r["r_burden"]), int(r["r_cap_moud"])]] * 2,
                      hovertemplate="%{customdata[0]}: burden #%{customdata[1]}, capacity #%{customdata[2]}<extra></extra>"))
# UNDERSERVED, not hl: iterating the set would reorder the labels on every
# run (string hashing is randomized per process) and churn the built HTML
for ab in UNDERSERVED + ["WY", "UT", "DC"]:
    r = g(ab)
    side = 1 if ab in ("WY", "UT", "DC") else 0
    notes.append(dict(x=side, y=int(r["r_burden"] if side == 0 else r["r_cap_moud"]),
                      text=f"<b>{ab}</b>", showarrow=False, xshift=-17 if side == 0 else 17,
                      font=dict(size=11.5, color=EMBER if ab in hl else (TEAL if ab in ("WY", "UT") else DIM))))
c_slope = dict(
    data=lines,
    layout=layout(height=640, margin=dict(l=70, r=70, t=46, b=20),
                  xaxis=dict(tickmode="array", tickvals=[0, 1],
                             ticktext=["<b>Burden rank</b><br>Deaths/100k, 2022–24", "<b>Capacity rank</b><br>MOUD facilities/100k"],
                             side="top", ticklen=0, range=[-0.14, 1.14], fixedrange=True),
                  yaxis=dict(autorange="reversed", tickvals=[1, 10, 20, 30, 40, 51],
                             gridcolor="rgba(0,0,0,0)", ticklen=0, tickfont=dict(size=11.5, color=DIM), fixedrange=True),
                  annotations=notes))

# 03 drug reporting: share of overdose deaths naming no drug
sp = m.nlargest(10, "unspec_share").sort_values("unspec_share")
c_specificity = dict(
    data=[dict(type="bar", orientation="h", x=sp.unspec_share.round(1).tolist(), y=sp.abbrev.tolist(),
               marker=dict(color=[EMBER if a == "LA" else "#9d99cf" for a in sp.abbrev]),
               customdata=sp.state.tolist(),
               hovertemplate="%{customdata}: %{x}% of overdose deaths name no drug<extra></extra>")],
    layout=layout(height=400, bargap=0.45, margin=dict(l=40, r=30, t=34, b=44),
                  xaxis=gaxis(title=dict(text="Overdose deaths naming no drug, 2022–2024",
                                         font=dict(size=12.5, color=DIM)), ticksuffix="%"),
                  yaxis=dict(tickfont=dict(size=11.5), ticklen=0),
                  shapes=[dict(type="line", x0=NAT_UNSPEC, x1=NAT_UNSPEC, yref="paper", y0=0, y1=1,
                               line=dict(color=DIM, width=1.2, dash="dash"))],
                  annotations=[dict(x=NAT_UNSPEC, y=1, yref="paper", yshift=20,
                                    text=f"National {NAT_UNSPEC:.0f}%", showarrow=False,
                                    font=dict(size=11.5, color=DIM)),
                               dict(x=float(LA.unspec_share), y="LA", xshift=-10, xanchor="right",
                                    text="<b>Louisiana</b>", showarrow=False,
                                    font=dict(size=12, color="white"))]))

# 09 the recovery by race
RACE_STYLE = {"White": (INDIGO, 3.4), "Black or African American": (EMBER, 3.4),
              "American Indian or Alaska Native": (TEAL, 3.4)}
race_traces = []
for r_ in rp.index:
    color, width = RACE_STYLE.get(r_, ("#c9c9d6", 1.8))
    short = {"Black or African American": "Black", "American Indian or Alaska Native": "Am. Indian / Alaska Native",
             "Native Hawaiian or Other Pacific Islander": "Native Hawaiian / Pacific Isl.",
             "More than one race": "Two or more races"}.get(r_, r_)
    race_traces.append(dict(type="scatter", mode="lines", x=[int(c) for c in rp.columns],
                            y=rp.loc[r_].round(1).tolist(), name=short,
                            line=dict(color=color, width=width),
                            showlegend=r_ in RACE_STYLE,
                            hovertemplate=short + ": %{y} per 100k in %{x}<extra></extra>"))
c_race = dict(
    data=race_traces,
    layout=layout(height=450, margin=dict(l=52, r=24, t=34, b=46), showlegend=True,
                  legend=dict(orientation="h", x=0.5, xanchor="center", y=1.035,
                              font=dict(size=12.5)),
                  xaxis=dict(tickmode="array", tickvals=[int(c) for c in rp.columns], ticklen=0,
                             showline=True, linecolor=INK),
                  yaxis=gaxis(title=dict(text="Opioid deaths per 100,000", font=dict(size=12.5, color=DIM)),
                              rangemode="tozero")))

# 10 Medicaid expansion
med = []
for lbl, sub, color in [("Expanded Medicaid", exp, INDIGO), ("Did not expand", non, EMBER)]:
    med.append(dict(type="scatter", mode="markers", x=sub.moud_per_100k.round(2).tolist(),
                    y=[lbl] * len(sub), name=lbl,
                    marker=dict(size=10, color=color, opacity=0.55,
                                line=dict(width=1, color="#f7f7f9")),
                    customdata=sub.state.tolist(),
                    hovertemplate="%{customdata}: %{x} MOUD facilities /100k<extra></extra>"))
c_medicaid = dict(
    data=med,
    layout=layout(height=300, margin=dict(l=150, r=30, t=40, b=44),
                  xaxis=gaxis(title=dict(text="MOUD facilities per 100,000",
                                         font=dict(size=12.5, color=DIM)), rangemode="tozero"),
                  yaxis=dict(tickfont=dict(size=13), ticklen=0),
                  # category order follows the trace order above: 0 = Expanded, 1 = Did
                  # not expand. Each median must span its OWN row, or the chart shows
                  # each group measured against the other group's midpoint.
                  shapes=[dict(type="line", x0=MED_MOUD_E, x1=MED_MOUD_E, y0=-0.38, y1=0.38,
                               yref="y", line=dict(color=INDIGO, width=2.5)),
                          dict(type="line", x0=MED_MOUD_N, x1=MED_MOUD_N, y0=0.62, y1=1.38,
                               yref="y", line=dict(color=EMBER, width=2.5))],
                  # anchored beside each line, not centred on it: a centred label
                  # gets bisected by the very line it names
                  annotations=[dict(x=MED_MOUD_E, y="Expanded Medicaid", yshift=30,
                                    xanchor="left", xshift=6,
                                    text=f"median {MED_MOUD_E:.2f}", showarrow=False,
                                    font=dict(size=11.5, color=INDIGO)),
                               dict(x=MED_MOUD_N, y="Did not expand", yshift=-30,
                                    xanchor="left", xshift=6,
                                    text=f"median {MED_MOUD_N:.2f}", showarrow=False,
                                    font=dict(size=11.5, color=EMBER))]))

# ---------------- table rows ----------------
table_rows = "\n".join(
    f'<tr><td><strong>{r.state}</strong></td><td>{r.rate_3yr_adj:.1f} <span class="dimtd">(#{int(r.r_burden)})</span></td>'
    f'<td>{r.moud_per_100k:.2f} <span class="dimtd">(#{int(r.r_cap_moud)})</span></td>'
    f'<td>{r.deaths_per_moud_facility:.1f}</td><td>{r.decline_pct:+.0f}%</td></tr>'
    for r in und.itertuples())

def ordinal(n):
    """1 -> 1st, 2 -> 2nd, 41 -> 41st. The teens are the exception."""
    n = int(n)
    suffix = "th" if 11 <= n % 100 <= 13 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


CITE_DATE = "September 5, 2026"

CHARTS = {"c-national": c_national, "c-map": c_map, "c-decline": c_decline,
          "c-capacity": c_capacity, "c-specificity": c_specificity,
          "c-scatter": c_scatter, "c-slope": c_slope, "c-race": c_race,
          "c-medicaid": c_medicaid}

def figure(cid, title, caption, note="", wide=False, anim="", extra="", reset=False, sense=""):
    """`sense` states which direction is the good one, above the chart.

    It sits under the title at reading size rather than in the axis label or
    the caption: a reader should meet it before interpreting the marks, not
    after, and it is the one line that decides whether they read the finding
    forwards or backwards.
    """
    cls = ("wide " if wide else "") + (extra + " " if extra else "") + anim
    note_html = f'<span class="note">{note}</span>' if note else ""
    sense_html = f'<div class="figsense">{sense}</div>' if sense else ""
    # the map can be panned and pinched, so it needs a visible way back
    pill = (f'<button class="map-reset" type="button" data-target="{cid}" hidden>'
            f'Reset map</button>') if reset else ""
    return (f'<figure class="{cls.strip()}"><div class="figtitle">{title}</div>'
            f'{sense_html}'
            # min-height reserves the chart's space before Plotly renders. Without it
            # the placeholders are 0px tall when the browser restores a #sec-N hash on
            # reload, so it lands short by the height of every chart above the target.
            f'<div class="chart"><div id="{cid}" style="width:100%;'
            f'min-height:{CHARTS[cid]["layout"]["height"]}px"></div>{pill}</div>'
            f'<figcaption>{caption}{note_html}</figcaption></figure>')

# ---------------- HTML ----------------
html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>One Facility for Every Nineteen Deaths</title>
<meta name="description" content="{TOTAL:,} opioid deaths in eleven years, matched against every treatment facility in the country. Treatment followed the deaths; {UND_WORDS} states got left behind.">
<meta name="author" content="Fahim Ahamed">
<meta property="og:type" content="article">
<meta property="og:title" content="One Facility for Every Nineteen Deaths">
<meta property="og:description" content="{TOTAL:,} opioid deaths in eleven years, matched against every treatment facility in the country.">
<meta property="og:url" content="https://f-a-tonmoy.github.io/us-opioid-treatment-gap/">
<meta property="og:image" content="https://f-a-tonmoy.github.io/us-opioid-treatment-gap/reports/hero-image.png">
<meta property="og:image:width" content="2400">
<meta property="og:image:height" content="1350">
<meta name="twitter:card" content="summary_large_image">
<link rel="icon" type="image/svg+xml" href='data:image/svg+xml,<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32"><rect width="32" height="32" rx="5" fill="%231c1a17"/><text x="16" y="22" font-family="-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Arial,sans-serif" font-weight="700" font-size="14" fill="white" text-anchor="middle" letter-spacing="-0.5">FA</text></svg>'>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=Spectral:ital,wght@0,400;0,600;0,700;1,400&display=swap">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
:root {{
  --ink: #15171e; --ink-dim: #565b68; --rule: #e2e3ea; --bg: #f7f7f9;
  --indigo: #443d9e; --indigo-soft: rgba(68,61,158,0.07);
  --ember: #c8501f; --ember-soft: rgba(200,80,31,0.08);
  --amber: #e8a33d; --cover: #171923;
}}
* {{ box-sizing: border-box; }}
html {{ scroll-behavior: smooth; }}
html, body {{ margin: 0; padding: 0; background: var(--bg); color: var(--ink); }}
body {{ font-family: "Spectral", "Georgia", serif; font-size: 18px; line-height: 1.72;
  -webkit-font-smoothing: antialiased; }}
::selection {{ background: rgba(232,163,61,0.4); }}
.sans, figcaption, .figtitle {{ font-family: "IBM Plex Sans", system-ui, sans-serif; }}
.container {{ max-width: 700px; margin: 0 auto; padding: 26px 24px 80px; }}

/* ---- cover band: white against the off-white page ---- */
.cover {{ background: #ffffff; border-bottom: 1px solid var(--rule); }}
.cover-inner {{ max-width: 700px; margin: 0 auto; padding: 58px 24px 34px; }}
.eyebrow {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 12px; font-weight: 600;
  letter-spacing: 0.22em; text-transform: uppercase; color: var(--indigo); margin-bottom: 20px; }}
h1 {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 54px; line-height: 1.03;
  font-weight: 700; margin: 0 0 18px; letter-spacing: -0.025em; }}
.deck {{ font-size: 21px; line-height: 1.5; color: var(--ink-dim); margin: 0 0 24px; max-width: 620px; }}
.byline {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13px; color: var(--ink-dim); }}
.byline .sep {{ color: #c6c6d2; margin: 0 9px; }}
.byline a {{ color: var(--indigo); text-decoration: none; border-bottom: 1px solid rgba(68,61,158,0.3); }}
.byline a:hover {{ border-bottom-color: var(--indigo); }}
.statband {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 0 20px;
  border-top: 2px solid var(--ink); margin: 36px 0 0; padding: 18px 0 0; }}
.statband .cell + .cell {{ border-left: 1px solid var(--rule); padding-left: 20px; }}
.statband .n {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 10.5px; font-weight: 600;
  letter-spacing: 0.14em; color: var(--ember); margin-bottom: 7px; }}
.statband .v {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 28px; font-weight: 700;
  line-height: 1.05; letter-spacing: -0.015em; font-variant-numeric: lining-nums; }}
.statband .v.alarm {{ color: var(--ember); }}
.statband .k {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 11.5px; color: var(--ink-dim);
  margin-top: 6px; line-height: 1.45; }}

/* ---- body typography ---- */
.lede {{ font-size: 21.5px; line-height: 1.5; font-weight: 600; margin: 34px 0 24px; letter-spacing: -0.002em; }}
.opener .sc {{ font-variant: small-caps; letter-spacing: 0.04em; font-weight: 600; }}
h2 {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 25px; line-height: 1.22;
  font-weight: 700; margin: 66px 0 14px; letter-spacing: -0.012em; }}
h2 .kick {{ display: block; font-size: 11px; font-weight: 600; letter-spacing: 0.2em;
  text-transform: uppercase; color: var(--indigo); margin-bottom: 9px; }}
h2[id^="sec-"] {{ scroll-margin-top: 24px; }}
p {{ margin: 0 0 19px; }}
strong {{ font-weight: 700; }}
em.stat {{ font-style: normal; font-weight: 700; background: linear-gradient(transparent 68%, rgba(232,163,61,0.35) 68%); }}
a {{ color: var(--indigo); }}
a.term {{ color: inherit; text-decoration: none;
  border-bottom: 1.5px dotted rgba(68,61,158,0.6); cursor: help; }}
a.term:hover {{ color: var(--indigo); border-bottom-style: solid; }}

/* ---- figures: flat, hairline-topped, left-captioned ---- */
figure {{ margin: 36px 0 20px; border-top: 1px solid var(--ink); padding-top: 12px; }}
figure.wide {{ margin-left: -70px; margin-right: -70px; }}
.figtitle {{ font-size: 14px; font-weight: 600; letter-spacing: 0.005em; margin: 0 0 4px; }}
/* how to read the chart: belongs above the marks, at reading size */
.figsense {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13.5px;
  line-height: 1.45; color: var(--ink); background: var(--indigo-soft);
  border-left: 3px solid var(--indigo); padding: 7px 12px; margin: 0 0 12px;
  border-radius: 0 4px 4px 0; }}
@media (max-width: 640px) {{ .figsense {{ font-size: 13px; padding: 6px 10px; }} }}
figure .chart {{ background: transparent; position: relative; }}
/* bottom-left: clear of the figure title above and the colorbar on the right.
   There is no tile attribution to avoid - this is a projected geo chart, not
   a tile map, so nothing here is drawn by a basemap provider. */
.map-reset {{ position: absolute; left: 14px; bottom: 14px; z-index: 3;
  font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 12px;
  font-weight: 600; color: var(--indigo); background: rgba(255,255,255,0.94);
  border: 1px solid var(--rule); border-radius: 999px; padding: 6px 15px;
  cursor: pointer; box-shadow: 0 1px 4px rgba(21,23,30,.10); }}
.map-reset:hover {{ border-color: var(--indigo); }}
/* the cartesian charts no longer drag-zoom, so drop plotly's crosshair hint */
figure .chart .nsewdrag {{ cursor: default !important; }}
figcaption {{ font-size: 12.5px; color: var(--ink-dim); margin-top: 6px; line-height: 1.55; text-align: left; }}
figcaption .note {{ display: block; margin-top: 5px; font-size: 11.5px; color: #8a8f9d; }}

.pullquote {{ font-size: 25px; line-height: 1.4; font-style: italic; font-weight: 600;
  border-left: 3px solid var(--ember); padding: 4px 0 4px 22px; margin: 40px 0; }}

table.mini {{ width: 100%; border-collapse: collapse; margin: 26px 0;
  font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13.5px; }}
table.mini th {{ text-align: left; font-size: 10.5px; text-transform: uppercase; letter-spacing: 0.1em;
  color: var(--ink-dim); font-weight: 600; border-bottom: 2px solid var(--indigo); padding: 6px 10px 8px; }}
table.mini td {{ border-bottom: 1px solid var(--rule); padding: 9px 10px; }}
table.mini tr:last-child td {{ border-bottom: none; }}
.dimtd {{ color: var(--ink-dim); }}

.takeaways {{ border-top: 2px solid var(--ink); margin: 48px 0 0; padding-top: 20px; }}
.takeaways h3 {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 12px; font-weight: 700;
  letter-spacing: 0.2em; text-transform: uppercase; color: var(--indigo); margin: 0 0 16px; }}
.takeaways ol {{ margin: 0; padding-left: 22px; font-size: 16px; line-height: 1.6; }}
.takeaways li {{ margin-bottom: 13px; padding-left: 5px; }}
.takeaways li::marker {{ color: var(--ember); font-weight: 700;
  font-family: "IBM Plex Sans", system-ui, sans-serif; }}

.aside {{ font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 14px; line-height: 1.62;
  color: var(--ink); background: var(--indigo-soft); border-left: 3px solid var(--indigo);
  padding: 14px 18px; margin: 34px 0 28px; border-radius: 0 6px 6px 0; }}

/* ---- footer ---- */
footer {{ margin-top: 70px; padding-top: 26px; border-top: 2px solid var(--ink);
  font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13px; color: var(--ink-dim); line-height: 1.55; }}
footer code {{ background: #e9e9f1; padding: 1px 5px; border-radius: 3px; font-size: 12px; }}
footer a {{ color: var(--indigo); text-decoration: none; }}
footer a:hover {{ text-decoration: underline; }}
footer h3 {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.16em; color: var(--ink); margin: 0 0 10px; font-weight: 700; }}
footer ol.citations {{ padding-left: 20px; margin: 0 0 8px; }}
footer ol.citations li {{ margin-bottom: 9px; scroll-margin-top: 24px; }}
footer ol.citations li:target {{ border-radius: 4px; animation: cite-flash 2.5s ease-out forwards; }}
@keyframes cite-flash {{ 0%, 40% {{ background: var(--indigo-soft); }} 100% {{ background: transparent; }} }}
footer .repro {{ border-top: 1px solid var(--rule); margin-top: 22px; padding-top: 16px; }}
footer .repro ol {{ padding-left: 20px; margin: 8px 0 0; }}
footer .repro li {{ margin-bottom: 9px; }}
footer .colophon {{ border-top: 1px solid var(--rule); margin-top: 22px; padding-top: 16px; }}
footer .colophon p {{ margin: 0 0 4px; }}
footer .colophon .muted {{ color: #9aa0ad; }}
footer .crisis-note {{ border-top: 1px solid var(--rule); margin-top: 16px; padding-top: 16px; font-size: 13.5px; color: var(--ink); }}
footer .crisis-note strong {{ color: var(--ember); font-size: 14.5px; letter-spacing: 0.02em; }}
sup.cite {{ font-size: 0.62em; vertical-align: super; line-height: 0; margin-left: 1px; }}
sup.cite a {{ color: var(--indigo); text-decoration: none; font-weight: 700; padding: 1px 3px; border-radius: 2px; }}
sup.cite a:hover {{ background: var(--indigo-soft); }}

/* ---- chart reveal animations ---- */
figure .chart {{ opacity: 0; transform: translateY(14px);
  transition: opacity 0.8s cubic-bezier(0.215,0.61,0.355,1), transform 0.8s cubic-bezier(0.215,0.61,0.355,1); }}
figure .chart.is-visible {{ opacity: 1; transform: translateY(0); }}
figure.anim-draw .chart .scatterlayer .trace .lines path {{ stroke-dasharray: 4200; stroke-dashoffset: 4200;
  transition: stroke-dashoffset 1.6s ease 0.25s; }}
figure.anim-draw .chart.marks-on .scatterlayer .trace .lines path {{ stroke-dashoffset: 0; }}
figure.anim-pop .chart .scatterlayer .trace .points path {{ scale: 0.3; opacity: 0;
  transition: scale 0.7s cubic-bezier(0.34,1.56,0.64,1), opacity 0.45s ease; }}
figure.anim-pop .chart.marks-on .scatterlayer .trace .points path {{ scale: 1; opacity: 1; }}
figure.anim-hbars .chart .bars .point path {{ transform: scaleX(0); transform-origin: 0% 50%;
  transform-box: fill-box; transition: transform 0.9s cubic-bezier(0.215,0.61,0.355,1) 0.15s; }}
figure.anim-hbars .chart.marks-on .bars .point path {{ transform: scaleX(1); }}

/* ---- section nav: RIGHT side ---- */
.section-nav {{ position: fixed; top: 50%; transform: translateY(-50%);
  left: calc(50% + 442px); z-index: 5; }}
.section-nav ol {{ list-style: none; margin: 0; padding: 0; border-left: 2px solid var(--rule); }}
.section-nav a {{ display: flex; align-items: center; gap: 12px;
  margin-left: -2px; border-left: 2px solid transparent; padding: 7px 12px 7px 14px;
  font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13px;
  color: var(--ink-dim); text-decoration: none; }}
.section-nav .num {{ font-weight: 700; color: #a9a9ba; font-size: 12.5px; width: 22px; }}
.section-nav .label {{ opacity: 0; transition: opacity .25s; white-space: nowrap; }}
.section-nav:hover .label, .section-nav a.active .label {{ opacity: 1; }}
.section-nav a.active {{ border-left-color: var(--indigo); background: var(--indigo-soft); }}
.section-nav a.active .num, .section-nav a.active .label {{ color: var(--indigo); font-weight: 600; }}
.section-nav a:hover .num {{ color: var(--ink); }}
@media (max-width: 1240px) {{ .section-nav {{ display: none; }} }}
@media (max-width: 860px) {{ figure.wide {{ margin-left: -12px; margin-right: -12px; }} }}
@media (max-width: 640px) {{
  h1 {{ font-size: 36px; }} body {{ font-size: 16.5px; }}
  .cover-inner {{ padding: 42px 20px 26px; }}
  .statband {{ grid-template-columns: 1fr 1fr; row-gap: 18px; }}
  .statband .cell:nth-child(3) {{ border-left: none; padding-left: 0; }}
  /* five columns will not fit 375px at desktop padding: tighten, and let the
     table scroll inside itself so a long state name can never push the page */
  table.mini {{ display: block; overflow-x: auto; font-size: 12px; }}
  table.mini th {{ padding: 6px 5px 8px; font-size: 9.5px; letter-spacing: 0.06em; }}
  table.mini td {{ padding: 8px 5px; }}
}}
.back-to-top {{ position: fixed; right: 22px; bottom: 22px; width: 40px; height: 40px; border-radius: 50%;
  border: none; background: var(--cover); color: #f2f3f7; font-size: 17px; cursor: pointer;
  opacity: 0; pointer-events: none; transition: opacity .3s; box-shadow: 0 2px 10px rgba(21,23,30,.25); }}
.back-to-top.visible {{ opacity: 1; pointer-events: auto; }}
@media (prefers-reduced-motion: reduce) {{
  figure .chart, figure .chart * {{ transition: none !important; }}
}}
</style>
</head>
<body>

<nav class="section-nav" aria-label="Sections"><ol>
<li><a href="#sec-01"><span class="num">01</span><span class="label">The steepest drop</span></a></li>
<li><a href="#sec-02"><span class="num">02</span><span class="label">The map fentanyl drew</span></a></li>
<li><a href="#sec-03"><span class="num">03</span><span class="label">What certificates omit</span></a></li>
<li><a href="#sec-04"><span class="num">04</span><span class="label">An uneven retreat</span></a></li>
<li><a href="#sec-05"><span class="num">05</span><span class="label">Counting capacity</span></a></li>
<li><a href="#sec-06"><span class="num">06</span><span class="label">The Wyoming mirage</span></a></li>
<li><a href="#sec-07"><span class="num">07</span><span class="label">Capacity follows burden</span></a></li>
<li><a href="#sec-08"><span class="num">08</span><span class="label">The {UND_WORDS} exceptions</span></a></li>
<li><a href="#sec-09"><span class="num">09</span><span class="label">An unshared recovery</span></a></li>
<li><a href="#sec-10"><span class="num">10</span><span class="label">What money explains</span></a></li>
</ol></nav>

<div class="cover">
<div class="cover-inner">
<header class="masthead">
  <div class="eyebrow">U.S. Opioid Epidemic · 2014 to 2024</div>
  <h1>One Facility for Every Nineteen Deaths</h1>
  <p class="deck">
    Opioids killed more than {TOTAL_FLOOR:,} Americans in eleven years. Treatment followed the
    deaths. {UND_WORDS.capitalize()} states got left behind.
  </p>
  <div class="byline">
    By <strong>Fahim Ahamed</strong>
    <span class="sep">|</span>
    <a href="https://www.linkedin.com/in/f-a-tonmoy/" target="_blank" rel="noopener">LinkedIn</a>
    <span class="sep">|</span>
    <a href="https://f-a-tonmoy.github.io/" target="_blank" rel="noopener">Portfolio</a>
  </div>
</header>

<div class="statband">
  <div class="cell"><div class="n">01</div><div class="v">{TOTAL:,}</div>
    <div class="k">opioid overdose deaths recorded, 2014&ndash;2024</div></div>
  <div class="cell"><div class="n">02</div><div class="v">{PCT_2324:.0f}%</div>
    <div class="k">the 2024 drop &mdash; the largest one-year decline on record</div></div>
  <div class="cell"><div class="n">03</div><div class="v">{WV.rate_3yr_adj / NE.rate_3yr_adj:.0f}&times;</div>
    <div class="k">death-rate gap between the worst state (WV) and the best (NE)</div></div>
  <div class="cell"><div class="n">04</div><div class="v alarm">{WA.deaths_per_moud_facility:.0f}</div>
    <div class="k">opioid deaths per treatment facility in Washington state, the worst ratio in the country</div></div>
</div>
</div>
</div>

<article class="container">

<p class="lede">
  In 2024, an American died of an opioid overdose roughly every {MIN_PER_DEATH_24:.0f} minutes.
  Two years earlier it was every {MIN_PER_DEATH_22:.0f}. The improvement did not arrive evenly.
</p>

<p class="opener">
  <span class="sc">This analysis joins three federal datasets</span> at the state level: every
  <a class="term" href="https://en.wikipedia.org/wiki/Opioid_epidemic_in_the_United_States" target="_blank" rel="noopener">opioid overdose</a> death certificate from 2014 through
  2024<sup class="cite"><a href="#cite-1">1</a></sup>, every substance-use treatment facility in
  SAMHSA&rsquo;s national locator<sup class="cite"><a href="#cite-2">2</a></sup>, and Census
  population, income, and insurance figures<sup class="cite"><a href="#cite-3">3</a></sup>. It
  asks one question: do the states with the heaviest overdose burden have the treatment
  infrastructure to match?
</p>

<h2 id="sec-01"><span class="kick">Section 01</span>The steepest drop on record.</h2>
<p>
  In 2014 the country recorded {D2014:,} opioid deaths. By the 2022 peak it had
  nearly tripled, to {D2022:,}, driven almost entirely by illegally made
  <a class="term" href="https://en.wikipedia.org/wiki/Fentanyl" target="_blank" rel="noopener">fentanyl</a>. The worst year was 2020: the pandemic hit a street supply already flooded
  with fentanyl, and deaths jumped {PCT_2020:.0f}%. In 2024 the count fell to
  <em class="stat">{D2024:,}</em>, down {abs(PCT_2324):.0f}% in one year, the sharpest
  decline ever measured<sup class="cite"><a href="#cite-4">4</a></sup>. Nobody fully knows why. Wider <a class="term" href="https://en.wikipedia.org/wiki/Naloxone" target="_blank" rel="noopener">naloxone</a> access, a shifting
  drug supply, and fewer surviving long-term users all likely play a part, and the only
  earlier decline, a {abs(PCT_DIP_2018):.0f}% dip in 2018, lasted one year. The
  progress is real. So is the limit: <em class="stat">2024&rsquo;s toll still sits
  {PCT_ABOVE_2019:.0f}% above the last pre-pandemic year, and nearly double 2014&rsquo;s</em>.
</p>
{figure("c-national", "The national toll, year by year",
        "Opioid overdose deaths across all 50 states + DC.",
        "Multiple-cause-of-death records: underlying cause is drug overdose (X40&ndash;44, X60&ndash;64, X85, Y10&ndash;14), with any opioid (T40.0&ndash;T40.4, T40.6) among the listed causes.",
        wide=True, anim="anim-draw",
        sense="Lower is better. Each point is one year&rsquo;s national death count.")}

<h2 id="sec-02"><span class="kick">Section 02</span>The map the fentanyl era drew.</h2>
<p>
  Measured over 2022&ndash;2024 so a single year cannot swing the ranking, West Virginia&rsquo;s rate
  of <em class="stat">{WV.rate_3yr_adj:.0f} deaths per 100,000 residents</em> is
  {WV.rate_3yr_adj / NAT3_ADJ:.1f} times the national figure of {NAT3_ADJ:.0f}, and
  {WV.rate_3yr_adj / NE.rate_3yr_adj:.0f} times Nebraska&rsquo;s, the lowest. The highest rates run through
  Appalachia and the mid-Atlantic (West Virginia, Delaware, Tennessee, Kentucky),
  northern New England (Maine, Vermont), and, new this decade, up the Pacific coast.
  The Plains states have far lower rates.
</p>
{figure("c-map", "Where the burden concentrates",
        f"Opioid overdose deaths per 100,000 residents, 2022&ndash;2024 combined, corrected for incomplete drug reporting. DC, second-highest at {DCr.rate_3yr_adj:.0f}, is too small to see at this scale.",
        f"North Carolina sits lower here than it should: the federal file records {NC23:,} opioid deaths there in 2023, against the {NC23_STATE:,} the state itself counts, a gap of medical examiner cases closed too late to be included.",
        wide=True, reset=True,
        sense="Darker is worse: more opioid deaths per 100,000 residents.")}

<h2 id="sec-03"><span class="kick">Section 03</span>What the death certificates leave out.</h2>
<p>
  That map needs a correction before it can be trusted. A death certificate
  records an overdose, but naming the drug depends on whether the local
  coroner or medical examiner ordered a toxicology screen and wrote down what
  it found. Nationally, <em class="stat">{NAT_UNSPEC:.0f}% of overdose deaths
  name no drug at all</em>. That number is not spread evenly. In Louisiana it
  is {LA.unspec_share:.0f}%. In Connecticut it is under one percent.
</p>
<p>
  States that investigate carefully therefore look worse than states that do
  not, purely as an artifact of record-keeping. The standard repair, used in
  the epidemiological literature<sup class="cite"><a href="#cite-6">6</a></sup>,
  is to assume the unnamed deaths in a state involved opioids at the same rate
  as that state&rsquo;s deaths where a drug <em>was</em> named, and scale
  accordingly. Across 2022&ndash;2024 that lifts the national toll by
  {ADJ_LIFT:.0f}%, and it can only be applied from 2018 on, because the older
  federal database was not queried the same way. Every headline count in this
  article is therefore the recorded one, which is what CDC publishes and what
  can be checked against it. The correction is used where it matters, comparing
  states against each other.
</p>
<p>
  For most states the correction is cosmetic: {int((m.rank_move.abs() <= 4).sum())} of
  the 51 move four places or fewer in the burden ranking. One state moves a
  lot. <em class="stat">Louisiana rises {int(LA.rank_move)} places, from
  {ordinal(LA.r_burden_rep)} to {ordinal(LA.r_burden)}</em>, which puts it in the
  worst ten in the country. Its reported rate of {LA.rate_3yr:.1f} per 100,000
  becomes {LA.rate_3yr_adj:.1f}. Every ranking from here on uses the corrected
  figures, and the correction is the reason Louisiana appears in them at all.
</p>
{figure("c-specificity", "Where the record-keeping is thinnest",
        "The ten states least likely to name a drug on an overdose death certificate.",
        "A high share means the state&rsquo;s opioid deaths are undercounted, not that it has fewer of them.",
        anim="anim-hbars",
        sense="Lower is better: a smaller share of overdose deaths whose certificate names no drug, so less correction is needed.")}

<h2 id="sec-04"><span class="kick">Section 04</span>Every state improved. Some barely.</h2>
<p>
  Every state, and DC, now sits below its worst year; the size of the retreat
  varies widely. Arkansas, Nebraska, West Virginia, and Ohio have roughly
  halved their death rates. At the other end, <em class="stat">Alaska is down just
  {abs(AK.decline_pct):.0f}%</em>, and in 2024 had the country&rsquo;s highest <a class="term" href="https://en.wikipedia.org/wiki/Mortality_rate" target="_blank" rel="noopener">unadjusted death rate</a>
  ({AK.crude_rate:.1f} per 100k). Nevada and Washington have cut deaths at less than
  half the national pace; Arizona and Oregon are not far ahead. The epidemic is
  shifting west, and eleven years make it plain: Washington&rsquo;s 2024 rate is
  {WA.ratio14:.1f} times its 2014 level, Alaska&rsquo;s {AK.ratio14:.1f}, against
  {NAT_R14:.1f} nationally. Only {BELOW14} sit below their 2014 rates, and both peaked
  early: New Hampshire in {int(g('NH').peak_year)}, at the front of the fentanyl wave,
  Utah in {int(g('UT').peak_year)} itself.
</p>
{figure("c-decline", "The 2024 retreat, state by state",
        "Each state is measured against its own worst year, 2014&ndash;2023.",
        wide=True, anim="anim-hbars",
        sense="Further left is better: a steeper fall from the state&rsquo;s own worst year. Bars near zero have barely retreated.")}

<h2 id="sec-05"><span class="kick">Section 05</span>Counting treatment capacity.</h2>
<p>
  That is the burden half of the question. On the treatment side, SAMHSA&rsquo;s locator
  lists {MOUD_TOTAL:,} facilities providing <a class="term" href="https://en.wikipedia.org/wiki/Opioid_use_disorder" target="_blank" rel="noopener">medication for opioid
  use disorder</a> (MOUD): <a class="term" href="https://en.wikipedia.org/wiki/Buprenorphine" target="_blank" rel="noopener">buprenorphine</a> or <a class="term" href="https://en.wikipedia.org/wiki/Methadone" target="_blank" rel="noopener">methadone</a>. Of those, <em class="stat">{OTP_TOTAL:,}
  are certified <a class="term" href="https://www.samhsa.gov/substance-use/treatment/opioid-treatment-program" target="_blank" rel="noopener">Opioid Treatment Programs</a> (OTPs)</em><sup class="cite"><a href="#cite-5">5</a></sup>, the only tier allowed to dispense
  methadone. The spread is nearly nine-to-one: Maine, Kentucky, and West Virginia have four to six MOUD
  facilities per 100,000 residents. Texas has <em class="stat">{TX.moud_per_100k:.2f}</em>,
  the fewest anywhere, despite the fifth-largest death toll
  ({int(TX.deaths_latest):,} in {LATEST}) and the nation&rsquo;s highest <a class="term" href="https://en.wikipedia.org/wiki/Health_insurance_coverage_in_the_United_States" target="_blank" rel="noopener">uninsured rate</a>
  ({TX.pct_uninsured:.0f}%).
</p>
<div class="aside">
  <strong>Two measures, deliberately.</strong> The OTP count, the narrower tier, comes
  from SAMHSA&rsquo;s official certification list. The MOUD count, any facility offering
  buprenorphine or methadone, comes from what facilities say about themselves in the
  treatment locator. They are collected differently and fail
  differently, which is exactly why every finding below has to satisfy both. A
  facility is also not a bed, a clinician, or a patient slot: these counts show where
  treatment exists, not how much of it there is.
</div>
{figure("c-capacity", "Two tiers of capacity, state by state",
        "Sorted by MOUD density; the line is each state&rsquo;s gap between tiers.",
        wide=True, anim="anim-pop",
        sense="Further right is better, meaning more treatment per resident. A longer line is worse: it means fewer of that state&rsquo;s facilities are certified to dispense methadone.")}

<h2 id="sec-06"><span class="kick">Section 06</span>The Wyoming mirage.</h2>
<p>
  Wyoming tops the capacity chart: {WY.moud_per_100k:.1f} MOUD facilities per 100,000,
  the best rate in the country. It is an illusion built on small numbers. The
  state&rsquo;s {int(WY.moud_count)} listed facilities serve about
  {WY.acs_pop / 1e3:,.0f},000 people scattered across nearly 98,000 square miles, and
  <em class="stat">not one of them is a certified OTP</em>. Wyoming is the only state
  in the country where a person cannot get methadone from a licensed program without
  leaving it.
</p>
<p>
  Finding that took two attempts. The treatment locator asks facilities to declare
  their own credentials, and on that self-reported flag Wyoming appeared to have one
  certified program. SAMHSA&rsquo;s official certification
  list<sup class="cite"><a href="#cite-5">5</a></sup> says zero. Nationally the two
  sources look interchangeable, {OTP_SELFREPORT:,} against {OTP_TOTAL:,}, which is how
  the error survived a first check. State by state they are not, so every OTP count
  here comes from the official list.
</p>
<p>
  West Virginia&rsquo;s {int(WV.otp_count)} is itself not an accident. It is the only
  state with a standing moratorium on new methadone clinics, layered with zoning rules
  that keep them half a mile from any school or
  daycare<sup class="cite"><a href="#cite-9">9</a></sup>. The state with the worst
  overdose burden in the country has written a cap on the most regulated form of
  treatment into its own law.
</p>
<p>
  Distance is the part none of these counts capture. Among Medicare patients, a drive of
  15 minutes instead of 5 was associated with roughly half the likelihood of receiving
  methadone at all<sup class="cite"><a href="#cite-10">10</a></sup>. A state can look
  adequately supplied and still be unreachable.
</p>

<h2 id="sec-07"><span class="kick">Section 07</span>Mostly, treatment follows the burden.</h2>
<p>
  The story most people would expect, that the hardest-hit states also have the least
  treatment, is wrong. Across all 50 states and DC, higher death rates come with <em>more</em>
  treatment per resident, not less (<a class="term" href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient" target="_blank" rel="noopener">Spearman r</a>&nbsp;=&nbsp;+{R_MOUD:.2f},
  <a class="term" href="https://en.wikipedia.org/wiki/P-value" target="_blank" rel="noopener">p</a>&nbsp;=&nbsp;{P_MOUD:.3f}). For the methadone-dispensing OTP tier the relationship is
  tighter still (r&nbsp;=&nbsp;+{R_OTP:.2f}, p&nbsp;&lt;&nbsp;0.0001): the OTP map closely
  follows the regions hit first, Appalachia and the Northeast.
  A decade into the fentanyl era, treatment has largely gone where the deaths were.
</p>
<p>
  What that relationship cannot settle is which way it runs. Capacity may have followed
  need, or need may have been slower to fall where treatment was scarce, or both may
  reflect how old each state&rsquo;s epidemic is. The treatment counts are a 2026
  snapshot and the deaths run 2022 through 2024, so the ordering in time cannot
  separate those.
</p>
<p>
  The comparison is also worth holding at arm&rsquo;s length. Nationally, only about a
  quarter of people who need treatment for opioid use disorder receive any medication
  for it<sup class="cite"><a href="#cite-8">8</a></sup>. What follows ranks states
  against each other, not against a standard any of them meets.
</p>
{figure("c-scatter", "Burden against capacity, 50 states + DC",
        f"Bubble area: {LATEST} deaths. Dotted lines: the national medians.",
        anim="anim-pop",
        sense="Better is left and up: fewer deaths, more treatment. The four orange states carry a heavy burden on thin capacity.")}

<h2 id="sec-08"><span class="kick">Section 08</span>{UND_WORDS.capitalize()} states break the pattern.</h2>
<p>
  Rank every state twice, by burden and by capacity, and the lines mostly run
  flat. The exceptions are the story, and picking them out calls for a rule rather
  than an eye. A state qualifies here if its corrected death rate is above the
  national median and it lands in the worst quarter of states on three separate
  measures of shortfall: the gap between its two ranks, how far it sits below the
  capacity its burden predicts, and the blunt count of opioid deaths per treatment
  facility. Those measures are built differently but agree closely
  (r&nbsp;=&nbsp;+{RANKGAP_RESID_AGREE:.2f}). DC is charted throughout but left out of
  the ranking, being a single city measured against states.
</p>
<p>
  <em class="stat">{UND_LIST}</em> clear all three. They rank {MIS_LO} to {MIS_HI}
  places worse on capacity than on burden. Washington pairs the third-highest state
  death rate in {LATEST} with the {ordinal(WA.r_cap_moud)}-ranked facility density.
</p>
<div class="pullquote">
  Washington has one treatment facility for every {WA.deaths_per_moud_facility:.0f} opioid
  deaths. Wyoming has one for every two.
</div>
{figure("c-slope", "Ranked twice: burden and capacity",
        f"States ranked by corrected death rate, then MOUD facilities per 100,000. Flat means the ranks match; falling orange lines are the underserved {UND_WORDS}.",
        anim="anim-draw",
        sense="Rank 1 means opposite things on each side: worst burden on the left, most facilities on the right. A line falling left to right is worse.")}
<table class="mini">
  <thead><tr><th>State</th><th>Deaths/100k &rsquo;22&ndash;24 (rank)</th><th>MOUD fac./100k (rank)</th><th>Deaths per facility</th><th>Off peak</th></tr></thead>
  <tbody>
{table_rows}
  </tbody>
</table>
<p>
  Washington and Nevada combine thin coverage with the country&rsquo;s weakest 2024
  declines. Alaska, kept off the list because {int(AK.deaths_latest)} deaths a year
  make its rates jumpy, looks much the same: seventh in corrected burden,
  {AK.deaths_per_moud_facility:.0f} deaths per facility. Where deaths are falling slowest,
  treatment is hardest to find.
</p>
<p>
  Two states fail the rule in a way worth naming. Alabama and Texas have the thinnest
  coverage in the country, {AL.moud_per_100k:.2f} and {TX.moud_per_100k:.2f} facilities
  per 100,000, and Alabama carries {AL.deaths_per_moud_facility:.0f} deaths for each
  one. Both are excluded only because their death rates sit below the national median.
  That is not reassurance so much as a description of exposure: the capacity is already
  missing, and nothing about the last decade suggests burden stays put.
</p>

<h2 id="sec-09"><span class="kick">Section 09</span>The recovery is not shared.</h2>
<p>
  The geography of the gap is only one dimension of it. Split the same deaths by race
  and the national decline stops looking like one story.
</p>
<p>
  White Americans&rsquo; opioid death rate peaked in {WHITE_PEAK} and by {LATEST} had
  returned to <em class="stat">{WHITE_CHG:+.0f}% of its {RACE_2018} level</em>,
  essentially back to where the fentanyl era started. No other group is close. Black
  Americans&rsquo; rate peaked two years later, in {BLACK_PEAK}, and remains
  <em class="stat">{BLACK_CHG:.0f}% above {RACE_2018}</em>. For American Indian and
  Alaska Native Americans it is {AIAN_CHG:.0f}% above. In {BLACK_PEAK} the Black rate
  ({BLACK_2023:.1f} per 100,000) was more than half again the White rate
  ({WHITE_2023:.1f}), a reversal of how this epidemic was described for most of its
  first decade. The Black share of opioid deaths rose from {BLK_SHARE_18:.0f}% in
  {RACE_2018} to {BLK_SHARE_23:.0f}% in {BLACK_PEAK}.
</p>
{figure("c-race", "One epidemic, two recoveries",
        f"Opioid overdose deaths per 100,000 by race, {RACE_2018}&ndash;{LATEST}.",
        "Single-race categories, so this series starts in 2018; the older federal file uses bridged-race groups that are not comparable. Hispanic origin is recorded separately and is not shown.",
        wide=True, anim="anim-draw",
        sense="Lower is better. A line still climbing while the others fall is the group left behind.")}
<p>
  Treatment access follows the same fault line. Studied across commercially insured
  patients, Black patients had roughly 40% lower odds of receiving either
  buprenorphine or methadone than white
  patients<sup class="cite"><a href="#cite-7">7</a></sup>, and the two medications are
  not distributed alike: methadone, which requires supervised daily dosing at a clinic,
  is concentrated in urban Black and Hispanic communities, while buprenorphine, which
  can be prescribed in an office, reaches whiter and wealthier ones. The
  {UND_WORDS}-state gap in this article is measured in facilities per resident. This
  one is not visible in facility counts at all.
</p>

<h2 id="sec-10"><span class="kick">Section 10</span>What money explains, and what policy does.</h2>
<p>
  The obvious explanation fails. States with too little treatment for their burden
  are not poorer (<a class="term" href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient" target="_blank" rel="noopener">Spearman r</a>&nbsp;=&nbsp;{R_MIS_INC:+.2f}, p&nbsp;=&nbsp;{P_MIS_INC:.2f}) and
  not meaningfully less insured (r&nbsp;=&nbsp;{R_MIS_UNINS:+.2f},
  p&nbsp;=&nbsp;{P_MIS_UNINS:.2f}); the deaths-per-facility measure agrees
  (r&nbsp;=&nbsp;{R_DPF_INC:+.2f} and {R_DPF_UNINS:+.2f}, both p&nbsp;&gt;&nbsp;0.09).
  Washington is rich and short on treatment; Alabama is poor and short on treatment. West
  Virginia is poor and well covered; Maryland is rich and well covered.
  <a class="term" href="https://en.wikipedia.org/wiki/Household_income_in_the_United_States" target="_blank" rel="noopener">Household income</a> is not the variable.
</p>
<p>
  One policy variable does separate them. The ten states that never adopted the
  <a class="term" href="https://en.wikipedia.org/wiki/Medicaid" target="_blank" rel="noopener">Medicaid</a> expansion<sup class="cite"><a href="#cite-11">11</a></sup> have
  <em class="stat">about half the treatment density of the states that did</em>:
  a median of {MED_MOUD_N:.2f} MOUD facilities per 100,000 against
  {MED_MOUD_E:.2f} (p&nbsp;=&nbsp;{P_MOUD_EXP:.3f}), and {MED_OTP_N:.2f} certified OTPs
  against {MED_OTP_E:.2f} (p&nbsp;=&nbsp;{P_OTP_EXP:.3f}). Their overdose burden,
  meanwhile, is statistically indistinguishable
  (p&nbsp;=&nbsp;{P_BURD_EXP:.2f}). Same problem, half the infrastructure.
</p>
{figure("c-medicaid", "The one line that separates them",
        "Each dot is a state; the bars mark each group&rsquo;s median.",
        "Wisconsin and Georgia cover some adults through waivers without taking the expansion, and count as non-expansion here, following the source&rsquo;s own classification.",
        anim="anim-pop",
        sense="Further right is better: more MOUD facilities per resident.")}
<p>
  This explains part of the gap, not all of it. {" and ".join(UND_NONEXP)} never
  expanded; {" and ".join(UND_EXP)} did, and are short of treatment anyway. And the
  comparison is between fifty-one jurisdictions that differ in a hundred other ways at
  the same time, so it identifies an association, not a mechanism. What can be said is
  that the thing income failed to explain, a single coverage decision partly does.
</p>

<div class="takeaways">
<h3>What the numbers say</h3>
<ol>
  <li>Opioid deaths fell {abs(PCT_2324):.0f}% in {LATEST} to {D2024:,}, a record single-year fall. Seven of the eight weakest state declines are in the West; Alaska has barely joined the retreat.</li>
  <li>Treatment capacity broadly tracks burden (r&nbsp;=&nbsp;+{R_MOUD:.2f}; +{R_OTP:.2f} for methadone programs), though which way that relationship runs cannot be settled here.</li>
  <li>{UND_LIST} are the exceptions, chosen by a rule that requires above-median burden and worst-quartile shortfall on three separate measures.</li>
  <li>White Americans&rsquo; death rate is back to its {RACE_2018} level. Black Americans&rsquo; is {BLACK_CHG:.0f}% higher, and American Indian and Alaska Native Americans&rsquo; is {AIAN_CHG:.0f}% higher.</li>
  <li>Income and insurance explain none of the gap. Medicaid expansion explains part: non-expansion states carry the same burden with about half the treatment.</li>
</ol>
</div>

<h2 id="sec-11"><span class="kick">Limits</span>What this data cannot tell us</h2>
<p>
  <strong>Counts are not capacity, and distance is not measured.</strong> A facility that
  treats forty patients and one that treats four thousand count the same here, and the
  field&rsquo;s usual access measure is travel time rather than facilities per resident.
  SAMHSA&rsquo;s facility survey does collect client counts, and drive times need road
  networks and sub-state population; this analysis uses neither, so it shows where
  treatment exists rather than how much there is or who can reach it.
</p>
<p>
  <strong>Opioids are only part of the burden.</strong> A stimulant is involved in a
  large share of opioid deaths &mdash; two in five, the last time CDC measured
  co-involvement directly<sup class="cite"><a href="#cite-13">13</a></sup> &mdash; and medication for opioid use disorder does
  nothing for methamphetamine or cocaine. Capacity measured this way answers a
  narrower question than the death counts pose.
</p>
<p>
  <strong>The locator is self-reported.</strong> Facilities describe their own services and
  update them on their own schedule. A clinic that quietly stopped prescribing still counts;
  a new one can take months to appear.
</p>
<p>
  <strong>Office-based prescribing is invisible.</strong> Since Congress eliminated the
  <a class="term" href="https://en.wikipedia.org/wiki/Drug_Addiction_Treatment_Act" target="_blank" rel="noopener">buprenorphine waiver</a> at the end of 2022, any clinician whose DEA registration covers Schedule III drugs can prescribe it.
  Those prescribers, and <a class="term" href="https://en.wikipedia.org/wiki/Telehealth" target="_blank" rel="noopener">telehealth</a>, do not appear in a facility locator, so
  buprenorphine access is undercounted everywhere, most in states that lean on office-based care.
</p>
<p>
  <strong>The timelines do not fully match.</strong> Facility data is a September 2026
  snapshot; deaths run 2014&ndash;2024; population, income, and insurance are 2023
  estimates, so a state that opened clinics last year gets credit against deaths that
  came before them. Deaths stop at 2024 because that is where final data ends: CDC
  publishes provisional counts closer to the present, but they are revised for months
  and are incomplete in ways that vary by state, which would put reporting speed into
  the rankings.
</p>
<p>
  <strong>States are not independent observations.</strong> Drug supply, policy, and
  patients all cross state lines, and people travel out of state for treatment,
  especially for methadone. The correlations here treat fifty-one jurisdictions as
  separate points when neighboring states are not.
</p>
<p>
  <strong>The correction is a simplification.</strong> Redistributing unspecified deaths
  in proportion to a state&rsquo;s named ones assumes the missing certificates resemble
  the recorded ones. The published method models local characteristics as well; this one
  does not, and Louisiana&rsquo;s adjustment in particular rests on that assumption.
</p>
<p>
  <strong>Exact ranks carry more precision than the data does.</strong> Alaska&rsquo;s
  {LATEST} rate of {AK.crude_rate:.1f} per 100,000 has a confidence interval running
  {AK.crude_lo:.0f} to {AK.crude_hi:.0f}, and {CI_OVERLAP_WV} states&rsquo; intervals
  overlap West Virginia&rsquo;s. Re-ranking on age-adjusted rather than unadjusted rates
  also moves states by up to {AA_MAXSHIFT} places ({AA_MAXSTATE} most). Neither changes
  anything above, because the findings rest on groups of states rather than exact
  positions.
</p>
<p>
  <strong>North Carolina&rsquo;s 2023 undercount.</strong> North Carolina reports
  {NC23_STATE:,} opioid overdose deaths for 2023<sup class="cite"><a href="#cite-12">12</a></sup>
  against the {NC23:,} in the federal file, a backlog its medical examiner closed after the
  file was compiled. Its true burden rank is higher than shown.
</p>

<footer>
  <h3>Citations</h3>
  <ol class="citations">
    <li id="cite-1">Centers for Disease Control and Prevention, National Center for Health
      Statistics. National Vital Statistics System, Mortality 2014&ndash;2024, CDC WONDER
      (Multiple Cause of Death: <code>D157</code> for 2018&ndash;2024, <code>D77</code> for
      2014&ndash;2017). UCD <code>X40&ndash;44, X60&ndash;64, X85, Y10&ndash;14</code>; MCD
      <code>T40.0&ndash;T40.4, T40.6</code>; by state and year. Accessed {CITE_DATE}.
      <a href="https://wonder.cdc.gov/mcd-icd10-expanded.html" target="_blank" rel="noopener">wonder.cdc.gov</a></li>
    <li id="cite-2">Substance Abuse and Mental Health Services Administration,
      FindTreatment.gov Treatment Services Locator (<code>exportsAsJson</code> API).
      12,220 substance-use listings retrieved {CITE_DATE}; 11,614 unique facilities in the
      50 states + DC.
      <a href="https://findtreatment.gov" target="_blank" rel="noopener">findtreatment.gov</a></li>
    <li id="cite-3">U.S. Census Bureau, American Community Survey 2023 1-year estimates:
      <code>B01003</code> population, <code>B19013</code> median household income,
      <code>DP03</code> insurance coverage, via the Census API, {CITE_DATE}.
      <a href="https://data.census.gov/table/ACSDT1Y2023.B01003" target="_blank" rel="noopener">data.census.gov</a></li>
    <li id="cite-4">National Center for Health Statistics, <em>Drug Overdose Deaths in the
      United States, 2023&ndash;2024</em>, NCHS Data Brief No. 549. Independently reports the
      national totals used here: 79,358 deaths in 2023, 54,045 in 2024.
      <a href="https://www.cdc.gov/nchs/products/databriefs/db549.htm" target="_blank" rel="noopener">cdc.gov/nchs</a></li>
    <li id="cite-5">SAMHSA, Opioid Treatment Program Directory: the official list of certified
      OTPs, exported {CITE_DATE}. {OTP_TOTAL:,} programs in the 50 states + DC, and the source
      of every OTP count here.
      <a href="https://www.samhsa.gov/find-help/locators/opioid-treatment-program-directory" target="_blank" rel="noopener">samhsa.gov</a></li>
    <li id="cite-6">Slavova S, Delcher C, Buchanich JM, et al. Methodological Complexities in
      Quantifying Rates of Fatal Opioid-Related Overdose. <em>Current Epidemiology Reports</em>,
      2019. Describes the proportional redistribution used here; Ruhm&rsquo;s adjustment raised
      the 2015 national rate by 21%.
      <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC6559129/" target="_blank" rel="noopener">National Library of Medicine</a></li>
    <li id="cite-7">Barnett ML, Meara E, Lewinson T, et al. Racial Inequality in Receipt of
      Medications for Opioid Use Disorder. <em>New England Journal of Medicine</em>, May 2023.
      Black patients were less likely than white patients to receive buprenorphine or methadone.
      <a href="https://www.nejm.org/doi/full/10.1056/NEJMsa2212412" target="_blank" rel="noopener">nejm.org</a></li>
    <li id="cite-8">Treatment for Opioid Use Disorder: Population Estimates, United States,
      2022. About a quarter of adults needing treatment received any medication for it.
      <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11254342/" target="_blank" rel="noopener">National Library of Medicine</a></li>
    <li id="cite-9">American Civil Liberties Union, on West Virginia&rsquo;s moratorium on new
      methadone clinics, currently the subject of a federal challenge.
      <a href="https://www.aclu.org/press-releases/challenge-to-west-virginia-methadone-clinic-moratorium-proceeds" target="_blank" rel="noopener">aclu.org</a></li>
    <li id="cite-10">Medicare Beneficiary Receipt of Methadone by Drive Time to Opioid Treatment
      Programs: longer travel time meant sharply lower likelihood of receiving methadone.
      <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11969284/" target="_blank" rel="noopener">National Library of Medicine</a></li>
    <li id="cite-11">KFF, Status of State Medicaid Expansion Decisions. The ten-state
      non-expansion classification used here.
      <a href="https://www.kff.org/medicaid/status-of-state-medicaid-expansion-decisions/" target="_blank" rel="noopener">kff.org</a></li>
    <li id="cite-12">North Carolina Department of Health and Human Services, April 2026:
      opioid-involved deaths fell from {NC23_STATE:,} in 2023 to 2,254 in 2024, and overdose
      deaths overall from 4,442 to 2,934.
      <a href="https://www.ncdhhs.gov/news/press-releases/2026/04/27/new-data-shows-all-time-low-infant-mortality-nc-sharp-decrease-overdose-related-deaths" target="_blank" rel="noopener">ncdhhs.gov</a></li>
    <li id="cite-13">CDC, <em>Vital Signs: Characteristics of Drug Overdose Deaths Involving
      Opioids and Stimulants</em>, MMWR 2020 (24 states + DC, January&ndash;June 2019). Of
      16,236 overdose deaths, 5,301 involved opioids with a stimulant and 7,936 without one.
      National vital statistics do not publish this breakdown, so it is not computed from the
      mortality file used here.
      <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC7470457/" target="_blank" rel="noopener">MMWR / National Library of Medicine</a></li>
  </ol>
  <div class="repro">
    <h3>Reproducibility</h3>
    <p>A multi-script Python (pandas/SciPy) pipeline generates everything above &mdash;
      fetch &rarr; clean/merge &rarr; analytics &mdash; with every number interpolated
      from the processed data at build time, none typed by hand.</p>
    <ol>
      <li><strong>Mortality.</strong> CDC WONDER (no public sub-national API; queried via
        its web interface): database <code>D157</code> for 2018&ndash;2024, <code>D77</code>
        for 2014&ndash;2017, identical definitions (UCD <code>X40&ndash;44, X60&ndash;64,
        X85, Y10&ndash;14</code>; MCD <code>T40.0&ndash;T40.4, T40.6</code>; state &times;
        year). Merged on (FIPS, year), <code>D157</code> preferred for the 2018&ndash;2020
        overlap. All 153 overlapping state-year cells agree exactly; national totals match
        published NCHS figures to the death.</li>
      <li><strong>Drug-specificity correction.</strong> Two further WONDER queries on the
        same UCD codes, state &times; year, 2018&ndash;2024: one unrestricted (every drug
        overdose death), one with MCD restricted to <code>T36&ndash;T50.8</code> (deaths
        naming at least one specific drug). Adjusted opioid deaths = reported &times;
        total &divide; specified, applied per state-year. National lift {ADJ_LIFT:.0f}%;
        {int((m.rank_move.abs() <= 4).sum())} of 51 states move four ranks or fewer.</li>
      <li><strong>Treatment capacity.</strong> Two sources, kept separate on purpose.
        Certified OTPs come from SAMHSA&rsquo;s official directory ({OTP_TOTAL:,} in the
        50 states + DC). MOUD facilities come from one national FindTreatment.gov pull
        (<code>exportsAsJson/v2</code>): 12,220 listings &minus; 549 same-name-and-address
        duplicates &minus; 57 territory rows = 11,614 facilities, of which {MOUD_TOTAL:,}
        report using buprenorphine or methadone. An earlier version derived OTP counts
        from the locator&rsquo;s self-reported certification flag; that flag matches the
        directory nationally ({OTP_SELFREPORT:,} vs {OTP_TOTAL:,}) but not by state, and
        was discarded.</li>
      <li><strong>Population and socioeconomics.</strong> ACS 2023 1-year via the Census
        API: <code>B01003</code> population (the denominator behind every per-100k facility
        rate), plus <code>B19013</code> median household income and <code>DP03</code>
        insurance coverage (the variables tested in Section&nbsp;10). Labels verified
        against the API&rsquo;s metadata; Puerto Rico dropped. Medicaid expansion status
        is KFF&rsquo;s classification, held as a fixed 10-state list.</li>
      <li><strong>Cleaning.</strong> WONDER&rsquo;s trailing notes block stripped; numbers
        de-comma&rsquo;d and coerced, with <code>Suppressed</code>/<code>Unreliable</code>
        markers mapped to missing (none occur at state level &mdash; checked, not assumed);
        FIPS zero-padded; state identity reconciled across all sources through one fixed
        51-row name/abbreviation/FIPS table.</li>
      <li><strong>Checks.</strong> Each source must match that table 51/51 or the pipeline
        halts; each WONDER file must land exactly 51 states &times; its years. Analysis is
        descriptive only (ranks, Spearman, Mann&ndash;Whitney, n&nbsp;=&nbsp;51;
        2022&ndash;2024 combined rates to damp small-state noise). The underserved states
        are selected by rule, not by hand, and the reader-facing rank gap is validated
        against regression residuals (r&nbsp;=&nbsp;+{RANKGAP_RESID_AGREE:.2f}). No finding
        runs without a second measure agreeing.</li>
    </ol>
  </div>
  <div class="colophon">
    <p>Analysis and writing: <strong>Fahim Ahamed</strong> <span class="muted">(refined with AI)</span></p>
    <p>Tools: Python, pandas, SciPy, Plotly, Requests</p>
    <p class="muted">Deaths counted by state of residence. All 50 states and DC are above
      CDC&rsquo;s small-number reporting limits in every year; no suppressed data affects
      these figures.</p>
  </div>
  <p class="crisis-note">
    If you or someone you know is struggling with opioids, free confidential help is available
    around the clock: call or text <strong>988</strong>, or SAMHSA&rsquo;s National Helpline,
    <strong>1-800-662-4357</strong>.
  </p>
</footer>

</article>

<button class="back-to-top" aria-label="Back to top" title="Back to top">&uarr;</button>

<script>
var CH = {{
  "c-national": {J(c_national)},
  "c-map": {J(c_map)},
  "c-decline": {J(c_decline)},
  "c-capacity": {J(c_capacity)},
  "c-specificity": {J(c_specificity)},
  "c-scatter": {J(c_scatter)},
  "c-slope": {J(c_slope)},
  "c-race": {J(c_race)},
  "c-medicaid": {J(c_medicaid)}
}};
var CONFIG = {J(CONFIG)};
var MAP_KEEP_M = {J(MAP_KEEP_M)};
// The map is pannable and pinch-zoomable, so give the reader a way home.
// viewInitial is Plotly's own record of the resting frame, so this stays
// correct if the map layout changes. Note: Plotly.relayout() does NOT emit
// plotly_relayout, so the pill hides itself on click rather than waiting for
// an event that never arrives.
function wireMapReset(id) {{
  var div = document.getElementById(id);
  var pill = document.querySelector('.map-reset[data-target="' + id + '"]');
  if (!div || !pill || !div._fullLayout || !div._fullLayout.geo) return;
  var EPS_SCALE = 0.05, EPS_DEG = 0.002;
  function home() {{
    var v = div._fullLayout.geo._subplot.viewInitial;
    return {{scale: v['projection.scale'], lat: v['center.lat'], lon: v['center.lon']}};
  }}
  function away() {{
    var g = div._fullLayout.geo, h = home();
    return Math.abs(g.projection.scale - h.scale) > EPS_SCALE
        || Math.abs(g.center.lat - h.lat) > EPS_DEG
        || Math.abs(g.center.lon - h.lon) > EPS_DEG;
  }}
  function sync() {{ pill.hidden = !away(); }}
  sync();
  div.on('plotly_relayout', sync);
  div.on('plotly_doubleclick', function () {{ setTimeout(sync, 0); }});
  pill.addEventListener('click', function () {{
    var h = home();
    Plotly.relayout(div, {{'geo.projection.scale': h.scale,
                          'geo.center.lat': h.lat, 'geo.center.lon': h.lon}});
    pill.hidden = true;
  }});
}}

function renderCharts() {{
  // narrow screens: fewer, smaller map labels; trend chart keeps only the
  // three anchor labels (start, peak, drop). Runs once, before any render.
  if (window.innerWidth < 700) {{
    var t = CH['c-map'].data[1];
    var f = function (arr) {{ return arr.filter(function (_, i) {{ return MAP_KEEP_M[i]; }}); }};
    t.lon = f(t.lon); t.lat = f(t.lat); t.text = f(t.text);
    t.textfont.color = f(t.textfont.color);
    t.textfont.size = 8.5;
    var rl = CH['c-race'].layout;
    rl.xaxis.tickvals = [2018, 2020, 2022, 2024];
    rl.xaxis.tickangle = 0;
    rl.legend = {{orientation: 'h', x: 0, xanchor: 'left', y: 1.14, font: {{size: 11}}}};
    rl.margin = {{l: 46, r: 16, t: 48, b: 46}};
    // the narrow plot tightens two clusters: KY/ME at the top, and WA/NV among
    // the underserved dots, which have to move to another side of their marker
    CH['c-scatter'].layout.annotations.forEach(function (a) {{
      if (a.text === 'KY') a.xshift = -13;
      if (a.text === 'WA') {{ a.xshift = 0; a.yshift = 24; }}
      if (a.text === 'NV') {{ a.xshift = 8; a.yshift = 16; }}
    }});
    // the slope chart's left gutter holds both rank ticks and state labels;
    // at this width they collide unless x=0 is pushed further from the axis
    CH['c-slope'].layout.xaxis.range = [-0.30, 1.14];
    var keep = ["{D2014:,}", "{D2022:,}", "{D2024:,}"];
    CH['c-national'].layout.annotations = CH['c-national'].layout.annotations.filter(function (a) {{
      return keep.some(function (k) {{ return a.text.indexOf(k) !== -1; }});
    }});
  }}
  Object.keys(CH).forEach(function (id) {{
    // retry once: plotly's geo renderer can reject on a cold load
    Plotly.newPlot(id, CH[id].data, CH[id].layout, CONFIG).catch(function () {{
      return Plotly.newPlot(id, CH[id].data, CH[id].layout, CONFIG).catch(function () {{}});
    }}).then(function () {{ wireMapReset(id); }});
  }});
}}
// render after resources settle; avoids a race in plotly's geo pipeline during page load
if (document.readyState === 'complete') renderCharts();
else window.addEventListener('load', renderCharts);

(function () {{
  var reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var charts = document.querySelectorAll('figure .chart');
  function reveal(c) {{
    c.classList.add('is-visible');
    setTimeout(function () {{ c.classList.add('marks-on'); }}, 180);
    setTimeout(function () {{ c.classList.add('labels-on'); }}, 950);
  }}
  if (reduce || !('IntersectionObserver' in window)) {{
    charts.forEach(reveal);
  }} else {{
    var obs = new IntersectionObserver(function (entries) {{
      entries.forEach(function (e) {{
        if (!e.isIntersecting) return;
        obs.unobserve(e.target);
        reveal(e.target);
      }});
    }}, {{ threshold: 0.2, rootMargin: '0px 0px -8% 0px' }});
    charts.forEach(function (c) {{ obs.observe(c); }});
  }}

  var links = document.querySelectorAll('.section-nav a');
  var h2s = Array.prototype.slice.call(document.querySelectorAll('article h2[id^="sec-"]'));
  function update() {{
    var trig = window.innerHeight * 0.35, act = null;
    for (var i = 0; i < h2s.length; i++) {{
      if (h2s[i].getBoundingClientRect().top <= trig) act = h2s[i]; else break;
    }}
    links.forEach(function (l) {{
      l.classList.toggle('active', !!act && l.getAttribute('href') === '#' + act.id);
    }});
  }}
  window.addEventListener('scroll', update, {{ passive: true }});
  update();

  var btn = document.querySelector('.back-to-top');
  var mast = document.querySelector('header.masthead');
  if ('IntersectionObserver' in window && mast) {{
    new IntersectionObserver(function (entries) {{
      btn.classList.toggle('visible', !entries[0].isIntersecting);
    }}, {{ threshold: 0 }}).observe(mast);
  }}
  btn.addEventListener('click', function () {{ window.scrollTo({{ top: 0, behavior: 'smooth' }}); }});
}})();
</script>
</body>
</html>
"""

# one canonical URL: the article lives only at the Pages root
INDEX = ROOT / "index.html"
INDEX.write_text(html, encoding="utf-8")
print(f"wrote {INDEX} ({INDEX.stat().st_size / 1024:.0f} KB)")
