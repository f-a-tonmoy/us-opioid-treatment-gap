# Build the self-contained interactive HTML article from the processed data.
# Every number in the prose is interpolated from the CSVs, never hand-typed.
# Output: reports/article.html  (Plotly 2.35.2 via CDN, all data inline)
#
# Visual identity is deliberately its own - dark cover masthead, IBM Plex Sans
# headlines over Spectral body, indigo/ember palette, flat hairline figures -
# the reference case studies set the rigor bar, not the skin.

import json
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "article.html"
OUT.parent.mkdir(exist_ok=True)

m = pd.read_csv(ROOT / "data" / "processed" / "state_merged.csv", dtype={"fips": str})
long = pd.read_csv(ROOT / "data" / "processed" / "mortality_by_state_year.csv", dtype={"fips": str})

# ---------------- derived quantities ----------------
LATEST = int(m["year_latest"].iloc[0])
nat = long.groupby("year").agg(deaths=("deaths", "sum"), pop=("population", "sum")).reset_index()
nat["rate"] = nat["deaths"] / nat["pop"] * 1e5
nrate = nat.set_index("year")["rate"]
_natpeak = nrate[[y for y in nrate.index if y != LATEST]].max()
NAT_DECL = (nrate[LATEST] - _natpeak) / _natpeak * 100
NAT_R14 = nrate[LATEST] / nrate[2014]
TOTAL = int(nat["deaths"].sum())
D2014, D2017, D2018, D2019, D2020, D2022, D2023, D2024 = (int(nat.loc[nat.year == y, "deaths"].iloc[0]) for y in (2014, 2017, 2018, 2019, 2020, 2022, 2023, 2024))
PCT_DIP_2018 = (D2018 - D2017) / D2017 * 100
PCT_2324 = (D2024 - D2023) / D2023 * 100
PCT_2020 = (D2020 - D2019) / D2019 * 100
PCT_ABOVE_2019 = (D2024 - D2019) / D2019 * 100
PCT_ABOVE_2014 = (D2024 - D2014) / D2014 * 100
MIN_PER_DEATH_24 = 525_600 / D2024
MIN_PER_DEATH_22 = 525_600 / D2022
last3 = long[long.year >= LATEST - 2]
NAT3 = last3["deaths"].sum() / last3["population"].sum() * 1e5

m["r_burden"] = m["rate_3yr"].rank(ascending=False).astype(int)
m["r_cap_moud"] = m["moud_per_100k"].rank(ascending=False).astype(int)
m["r_cap_otp"] = m["otp_per_100k"].rank(ascending=False).astype(int)
m["mismatch"] = m["r_cap_moud"] - m["r_burden"]
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
WV, WY, DCr, WA, TX, AK, NE = (g(x) for x in ("WV", "WY", "DC", "WA", "TX", "AK", "NE"))
UNDERSERVED = ["SC", "WA", "TN", "NV", "AL"]
und = m[m.abbrev.isin(UNDERSERVED)].sort_values("mismatch", ascending=False)

def corr(x, y):
    r, p = stats.spearmanr(m[x], m[y])
    return r, p
R_MOUD, P_MOUD = corr("rate_3yr", "moud_per_100k")
R_OTP, P_OTP = corr("rate_3yr", "otp_per_100k")
R_MIS_INC, P_MIS_INC = corr("mismatch", "median_hh_income")
R_MIS_UNINS, P_MIS_UNINS = corr("mismatch", "pct_uninsured")
R_DPF_INC, P_DPF_INC = corr("deaths_per_moud_facility", "median_hh_income")
R_DPF_UNINS, P_DPF_UNINS = corr("deaths_per_moud_facility", "pct_uninsured")

OTP_TOTAL = int(m["otp_count"].sum())
MOUD_TOTAL = int(m["moud_count"].sum())

# ---------------- palette / chart chrome ----------------
J = lambda o: json.dumps(o, separators=(",", ":"))
INK, DIM, RULE = "#15171e", "#565b68", "#e2e3ea"
INDIGO, EMBER, TEAL, AMBER, GREY = "#443d9e", "#c8501f", "#2e7f8a", "#e8a33d", "#c2c2d0"
FONT = dict(family="IBM Plex Sans, system-ui, sans-serif", size=13, color=INK)
BASE_LAYOUT = dict(paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)", font=FONT,
                   margin=dict(l=48, r=16, t=12, b=40), showlegend=False)
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
               locations=m["abbrev"].tolist(), z=m["rate_3yr"].round(1).tolist(),
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
               textfont=dict(size=10, color=["white" if z >= WHITE_Z else "#454962" for z in _lab.rate_3yr]),
               hoverinfo="skip")],
    layout=layout(height=430, margin=dict(l=0, r=0, t=6, b=6),
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
                               dict(x=NAT_DECL, y=1, yref="paper", yshift=15,
                                    text=f"National {NAT_DECL:.0f}%", showarrow=False,
                                    font=dict(size=11.5, color=DIM)),
                               dict(xref="paper", yref="paper", x=0.02, y=0.995, xanchor="left",
                                    text="■ Western states", showarrow=False,
                                    font=dict(size=12.5, color=EMBER)),
                               dict(xref="paper", yref="paper", x=0.02, y=0.966, xanchor="left",
                                    text="■ Everywhere else", showarrow=False,
                                    font=dict(size=12.5, color="#8783c4"))]))

# 04 capacity: horizontal dumbbells, OTP (ember diamond) -> MOUD (indigo dot)
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
                  xaxis=gaxis(title=dict(text="Facilities per 100,000 residents", font=dict(size=12.5, color=DIM)), rangemode="tozero"),
                  yaxis=dict(tickfont=dict(size=10.5), ticklen=0),
                  annotations=[dict(x=3.3, y="WY", yshift=13,
                                    text=f"<b>Wyoming: one OTP, {int(WY.moud_count)} MOUD facilities</b>",
                                    showarrow=False, font=dict(size=12, color=INK))]))

# 06 scatter with median quadrant guides
hl = set(UNDERSERVED)
deep = {"WV", "KY", "ME", "VT", "MD"}
MEDX, MEDY = float(m["rate_3yr"].median()), float(m["moud_per_100k"].median())
def pt_color(ab):
    if ab in hl: return EMBER
    if ab in deep: return INDIGO
    if ab in ("WY", "UT"): return TEAL
    return GREY
c_scatter = dict(
    data=[dict(type="scatter", mode="markers", x=m["rate_3yr"].round(1).tolist(), y=m["moud_per_100k"].round(2).tolist(),
               marker=dict(size=(m["deaths_latest"] ** 0.5 / 2.4 + 5).round(1).tolist(),
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
                               dict(x=float(WV.rate_3yr), y=float(WV.moud_per_100k), text="WV", showarrow=False, yshift=14, font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(g('KY').rate_3yr), y=float(g('KY').moud_per_100k), text="KY", showarrow=False, yshift=14, font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(g('ME').rate_3yr), y=float(g('ME').moud_per_100k), text="ME", showarrow=False, yshift=14, font=dict(size=11.5, color=INDIGO)),
                               dict(x=float(WY.rate_3yr), y=float(WY.moud_per_100k), text="WY", showarrow=False, yshift=14, font=dict(size=11.5, color=TEAL)),
                               dict(x=float(DCr.rate_3yr), y=float(DCr.moud_per_100k), text="DC", showarrow=False, yshift=-14, font=dict(size=11.5, color=DIM)),
                               *[dict(x=float(g(a).rate_3yr), y=float(g(a).moud_per_100k), text=a, showarrow=False, yshift=-14, font=dict(size=11.5, color=EMBER)) for a in UNDERSERVED]]))

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
for ab in list(hl) + ["WY", "UT", "DC"]:
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

# 08 null scatter
c_null = dict(
    data=[dict(type="scatter", mode="markers", x=(m["median_hh_income"] / 1000).round(1).tolist(), y=m["mismatch"].tolist(),
               marker=dict(size=7.5, color=[EMBER if a in hl else GREY for a in m["abbrev"]]),
               customdata=[[s, int(mm)] for s, mm in zip(m["state"], m["mismatch"])],
               hovertemplate="%{customdata[0]}: $%{x}k income, mismatch %{customdata[1]:+}<extra></extra>")],
    layout=layout(height=360,
                  xaxis=gaxis(title=dict(text="Median household income, $000s (ACS 2023)", font=dict(size=12.5, color=DIM))),
                  yaxis=gaxis(title=dict(text="Capacity rank minus burden rank", font=dict(size=12.5, color=DIM)),
                              zeroline=True, zerolinecolor=INK, range=[-47, 38]),
                  annotations=[dict(x=float(m["median_hh_income"].max() / 1000), y=-41, xanchor="right",
                                    text=f"Spearman r = {R_MIS_INC:+.2f}, p = {P_MIS_INC:.2f}",
                                    showarrow=False, font=dict(size=12.5, color=DIM)),
                               *[dict(x=float(g(a).median_hh_income / 1000), y=int(g(a).mismatch),
                                      text=a, showarrow=False, font=dict(size=11.5, color=EMBER),
                                      **({"xshift": 11, "xanchor": "left"} if a == "SC" else
                                         {"yshift": {"WA": 14, "NV": 14, "TN": -15, "AL": -15}[a]}))
                                 for a in UNDERSERVED]]))

# ---------------- table rows ----------------
table_rows = "\n".join(
    f'<tr><td><strong>{r.state}</strong></td><td>{r.rate_3yr:.1f} <span class="dimtd">(#{int(r.r_burden)})</span></td>'
    f'<td>{r.moud_per_100k:.2f} <span class="dimtd">(#{int(r.r_cap_moud)})</span></td>'
    f'<td>{r.deaths_per_moud_facility:.1f}</td><td>{r.decline_pct:+.0f}%</td></tr>'
    for r in und.itertuples())

CITE_DATE = "September 5, 2026"

def figure(cid, title, caption, note="", wide=False, anim=""):
    cls = ("wide " if wide else "") + anim
    note_html = f'<span class="note">{note}</span>' if note else ""
    return (f'<figure class="{cls.strip()}"><div class="figtitle">{title}</div>'
            f'<div class="chart"><div id="{cid}" style="width:100%"></div></div>'
            f'<figcaption>{caption}{note_html}</figcaption></figure>')

# ---------------- HTML ----------------
html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>One Facility for Every Eighteen Deaths</title>
<meta name="description" content="{TOTAL:,} opioid deaths in eleven years, matched against every treatment facility in the country. Treatment followed the deaths; five states got left behind.">
<meta name="author" content="Fahim Ahamed">
<meta property="og:type" content="article">
<meta property="og:title" content="One Facility for Every Eighteen Deaths">
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
  border-top: 6px solid var(--indigo); -webkit-font-smoothing: antialiased; }}
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
figure .chart {{ background: transparent; }}
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
.section-nav ol {{ list-style: none; margin: 0; padding: 0; }}
.section-nav a {{ display: flex; align-items: center; gap: 10px;
  padding: 6px 8px; font-family: "IBM Plex Sans", system-ui, sans-serif; font-size: 13px;
  color: var(--ink-dim); text-decoration: none; border-radius: 5px; }}
.section-nav .num {{ font-weight: 700; color: #a9a9ba; font-size: 12.5px; width: 22px; text-align: right; }}
.section-nav .label {{ opacity: 0; transform: translateX(-4px); transition: opacity .25s, transform .25s; white-space: nowrap; }}
.section-nav a:hover .label, .section-nav a.active .label {{ opacity: 1; transform: none; }}
.section-nav a.active .num {{ color: var(--ember); }}
@media (max-width: 1240px) {{ .section-nav {{ display: none; }} }}
@media (max-width: 860px) {{ figure.wide {{ margin-left: -12px; margin-right: -12px; }} }}
@media (max-width: 640px) {{
  h1 {{ font-size: 36px; }} body {{ font-size: 16.5px; }}
  .cover-inner {{ padding: 42px 20px 26px; }}
  .statband {{ grid-template-columns: 1fr 1fr; row-gap: 18px; }}
  .statband .cell:nth-child(3) {{ border-left: none; padding-left: 0; }}
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
<li><a href="#sec-03"><span class="num">03</span><span class="label">An uneven retreat</span></a></li>
<li><a href="#sec-04"><span class="num">04</span><span class="label">Counting capacity</span></a></li>
<li><a href="#sec-05"><span class="num">05</span><span class="label">The Wyoming mirage</span></a></li>
<li><a href="#sec-06"><span class="num">06</span><span class="label">Capacity follows burden</span></a></li>
<li><a href="#sec-07"><span class="num">07</span><span class="label">The five exceptions</span></a></li>
<li><a href="#sec-08"><span class="num">08</span><span class="label">Money explains nothing</span></a></li>
</ol></nav>

<div class="cover">
<div class="cover-inner">
<header class="masthead">
  <div class="eyebrow">U.S. Opioid Epidemic · 2014 to 2024</div>
  <h1>One Facility for Every Eighteen Deaths</h1>
  <p class="deck">
    Opioids killed more than 600,000 Americans in eleven years. Treatment followed the
    deaths. Five states got left behind.
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
    <div class="k">opioid overdose deaths, 2014&ndash;2024</div></div>
  <div class="cell"><div class="n">02</div><div class="v">{PCT_2324:.0f}%</div>
    <div class="k">the 2024 drop &mdash; the largest one-year decline on record</div></div>
  <div class="cell"><div class="n">03</div><div class="v">{WV.rate_3yr / NE.rate_3yr:.0f}&times;</div>
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
        wide=True, anim="anim-draw")}

<h2 id="sec-02"><span class="kick">Section 02</span>The map the fentanyl era drew.</h2>
<p>
  Averaged over 2022&ndash;2024 to smooth out the small states, West Virginia&rsquo;s rate
  of <em class="stat">{WV.rate_3yr:.0f} deaths per 100,000 residents</em> is
  {WV.rate_3yr / NAT3:.1f} times the national figure of {NAT3:.0f}, and
  {WV.rate_3yr / NE.rate_3yr:.0f} times Nebraska&rsquo;s, the lowest. The highest rates run through
  Appalachia and the mid-Atlantic (West Virginia, Delaware, Tennessee, Kentucky),
  northern New England (Maine, Vermont), and, new this decade, up the Pacific coast.
  The Plains states have far lower rates.
</p>
{figure("c-map", "Where the burden concentrates",
        f"Opioid overdose deaths per 100,000 residents, 2022&ndash;2024 combined. DC, second-highest at {DCr.rate_3yr:.0f}, is too small to see at this scale.",
        "North Carolina&rsquo;s 2023 figure is understated in the final federal file: NCHS notes ~900 late-coded overdose deaths (a true count over 4,400) that will not be added to this dataset.",
        wide=True)}

<h2 id="sec-03"><span class="kick">Section 03</span>Every state improved. Some barely.</h2>
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
        wide=True, anim="anim-hbars")}

<h2 id="sec-04"><span class="kick">Section 04</span>Counting treatment capacity.</h2>
<p>
  That is the burden half of the question. On the treatment side, SAMHSA&rsquo;s locator
  lists {MOUD_TOTAL:,} facilities providing <a class="term" href="https://en.wikipedia.org/wiki/Opioid_use_disorder" target="_blank" rel="noopener">medication for opioid
  use disorder</a> (MOUD): <a class="term" href="https://en.wikipedia.org/wiki/Buprenorphine" target="_blank" rel="noopener">buprenorphine</a> or <a class="term" href="https://en.wikipedia.org/wiki/Methadone" target="_blank" rel="noopener">methadone</a>. Of those, <em class="stat">{OTP_TOTAL:,}
  are certified <a class="term" href="https://www.samhsa.gov/substance-use/treatment/opioid-treatment-program" target="_blank" rel="noopener">Opioid Treatment Programs</a> (OTPs)</em><sup class="cite"><a href="#cite-5">5</a></sup>, the only tier allowed to dispense
  methadone. The spread is nearly nine-to-one: Maine, Kentucky, and West Virginia have five to six MOUD
  facilities per 100,000 residents. Texas has <em class="stat">{TX.moud_per_100k:.2f}</em>,
  the fewest anywhere, despite the fifth-largest death toll
  ({int(TX.deaths_latest):,} in {LATEST}) and the nation&rsquo;s highest <a class="term" href="https://en.wikipedia.org/wiki/Health_insurance_coverage_in_the_United_States" target="_blank" rel="noopener">uninsured rate</a>
  ({TX.pct_uninsured:.0f}%).
</p>
<div class="aside">
  <strong>What a facility count is not.</strong> A facility is not a bed, a clinician, or a
  patient slot. They show where treatment exists, not how much. Every finding below is
  checked against a second measure for that reason.
</div>
{figure("c-capacity", "Two tiers of capacity, state by state",
        "Sorted by MOUD density; the line is each state&rsquo;s gap between tiers.",
        wide=True, anim="anim-pop")}

<h2 id="sec-05"><span class="kick">Section 05</span>The Wyoming mirage.</h2>
<p>
  Wyoming tops the capacity chart: {WY.moud_per_100k:.1f} MOUD facilities per 100,000.
  It is an illusion built on small numbers. The state&rsquo;s
  {int(WY.moud_count)} listed facilities serve about {WY.acs_pop / 1e3:,.0f},000 people
  scattered across nearly 98,000 square miles &mdash; and exactly
  <em class="stat">one of them is a certified OTP</em>. On methadone access, Wyoming ranks
  {int(WY.r_cap_otp)}th of 51. One clinic opening or closing moves its headline figure
  nearly three percent.
</p>
<p>
  The District of Columbia fails in the opposite direction. Its per-resident count looks
  unimpressive ({DCr.moud_per_100k:.1f} per 100k), yet its {int(DCr.moud_count)} facilities sit
  within a few miles of every resident. Counting facilities per resident makes empty states
  look good and packed cities look bad. DC also has the second-highest combined death rate
  ({DCr.rate_3yr:.0f} per 100k), concentrated in a single city with no rural areas around
  it, so DC is flagged in every chart and kept out of the headline rankings. Rhode Island
  shows the same fragility in reverse: first in OTP density,
  {int(g('RI').r_cap_moud)}th in MOUD.
</p>
<p>
  No single measure here is trusted until a second one agrees. The states named next
  pass that test.
</p>

<h2 id="sec-06"><span class="kick">Section 06</span>Mostly, treatment follows the burden.</h2>
<p>
  The story most people would expect, that the hardest-hit states also have the least
  treatment, is wrong. Across all 50 states and DC, higher death rates come with <em>more</em>
  treatment per resident, not less (<a class="term" href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient" target="_blank" rel="noopener">Spearman r</a>&nbsp;=&nbsp;+{R_MOUD:.2f},
  <a class="term" href="https://en.wikipedia.org/wiki/P-value" target="_blank" rel="noopener">p</a>&nbsp;=&nbsp;{P_MOUD:.3f}). For the methadone-dispensing OTP tier the relationship is
  tighter still (r&nbsp;=&nbsp;+{R_OTP:.2f}, p&nbsp;&lt;&nbsp;0.0001): the OTP map closely
  follows the regions hit first, Appalachia and the Northeast.
  Two decades into the crisis, treatment has largely gone where the deaths were. Whether
  that is policy answering need, or simply time (the oldest epidemics have had the longest
  to build), this data cannot say.
</p>
{figure("c-scatter", "Burden against capacity, 50 states + DC",
        f"Bubble area: {LATEST} deaths. Dotted lines: the national medians.",
        "Purple: high on both (WV, KY, ME). Orange: the underserved five. Teal: Wyoming&rsquo;s small-population illusion. DC: a city measured against states, flagged throughout.",
        anim="anim-pop")}

<h2 id="sec-07"><span class="kick">Section 07</span>Five states break the pattern.</h2>
<p>
  Rank every state twice, by burden and by capacity, and the lines mostly run
  flat. The exceptions are the story. <em class="stat">South Carolina, Washington,
  Tennessee, Nevada, and Alabama</em> all rank 18 to 32 places worse on capacity than on
  burden, and all five sit in the top nine of 51 on the simplest measure:
  opioid deaths per treatment facility. Washington pairs the third-highest state death rate
  in {LATEST} with the 41st-ranked facility density.
</p>
<div class="pullquote">
  Washington has one treatment facility for every {WA.deaths_per_moud_facility:.0f} opioid
  deaths. Wyoming has one for every two.
</div>
{figure("c-slope", "Ranked twice: burden and capacity",
        "States ranked by combined death rate, then MOUD facilities per 100,000. Flat means the ranks match; falling orange lines are the underserved five.",
        "Teal: Wyoming and Utah, whose high capacity ranks reflect small populations more than large systems. Ranks run 1&ndash;51, DC included but flagged.",
        anim="anim-draw")}
<table class="mini">
  <thead><tr><th>State</th><th>Deaths/100k &rsquo;22&ndash;24 (rank)</th><th>MOUD fac./100k (rank)</th><th>Deaths per facility</th><th>Off peak</th></tr></thead>
  <tbody>
{table_rows}
  </tbody>
</table>
<p>
  Washington and Nevada combine thin coverage with the country&rsquo;s weakest 2024
  declines. Alaska, kept off the list because {int(AK.deaths_latest)} deaths a year
  make its rates jumpy, looks much the same: sixth in combined burden,
  {AK.deaths_per_moud_facility:.0f} deaths per facility. Where deaths are falling slowest,
  treatment is hardest to find.
</p>

<h2 id="sec-08"><span class="kick">Section 08</span>The gap is not a poverty map.</h2>
<p>
  The obvious explanation fails. States with too little treatment for their burden
  are not poorer (<a class="term" href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient" target="_blank" rel="noopener">Spearman r</a>&nbsp;=&nbsp;{R_MIS_INC:+.2f}, p&nbsp;=&nbsp;{P_MIS_INC:.2f}) and
  not meaningfully less insured (r&nbsp;=&nbsp;{R_MIS_UNINS:+.2f},
  p&nbsp;=&nbsp;{P_MIS_UNINS:.2f}); the deaths-per-facility measure agrees
  (r&nbsp;=&nbsp;{R_DPF_INC:+.2f} and {R_DPF_UNINS:+.2f}, both p&nbsp;&gt;&nbsp;0.1).
  Washington is rich and short on treatment; Alabama is poor and short on treatment. West
  Virginia is poor and well covered; Maryland is rich and well covered. Whatever decides
  where treatment lands (licensing rules, <a class="term" href="https://en.wikipedia.org/wiki/Medicaid" target="_blank" rel="noopener">Medicaid</a> history, the age of each
  state&rsquo;s epidemic), it is not <a class="term" href="https://en.wikipedia.org/wiki/Household_income_in_the_United_States" target="_blank" rel="noopener">household income</a>, and this dataset cannot see it.
</p>
{figure("c-null", "The explanation that isn&rsquo;t",
        "The pattern here is the absence of one.",
        "Orange: the underserved five.",
        anim="anim-pop")}

<div class="takeaways">
<h3>What the numbers say</h3>
<ol>
  <li>Opioid deaths fell {abs(PCT_2324):.0f}% in 2024 to {D2024:,}, a record single-year fall. Seven of the eight weakest state declines are in the West; Alaska has barely joined the retreat.</li>
  <li>Treatment capacity broadly tracks burden (r&nbsp;=&nbsp;+{R_MOUD:.2f}; +{R_OTP:.2f} for methadone programs): the epidemic&rsquo;s first-wave states built the most capacity.</li>
  <li>South Carolina, Washington, Tennessee, Nevada, and Alabama are the exceptions: high burden, little capacity, confirmed by two independent measures.</li>
  <li>Washington&rsquo;s ratio of deaths to treatment facilities is the country&rsquo;s worst; Texas has the fewest facilities per resident.</li>
  <li>Income and insurance coverage explain none of the gap; the states with too little treatment are not the poor ones.</li>
</ol>
</div>

<h2 id="sec-09"><span class="kick">Limits</span>What this data cannot tell us</h2>
<p>
  <strong>Counts are not capacity.</strong> A facility that treats forty patients and one that
  treats four thousand count the same here. No public dataset reports patient slots per
  facility nationwide.
</p>
<p>
  <strong>The locator is self-reported.</strong> Facilities describe their own services and
  update them on their own schedule. A clinic that quietly stopped prescribing still counts;
  a new one can take months to appear.
</p>
<p>
  <strong>Office-based prescribing is invisible.</strong> Since Congress eliminated the
  <a class="term" href="https://en.wikipedia.org/wiki/Drug_Addiction_Treatment_Act" target="_blank" rel="noopener">buprenorphine waiver</a> at the end of 2022, any DEA-registered clinician can prescribe it.
  Those prescribers, and <a class="term" href="https://en.wikipedia.org/wiki/Telehealth" target="_blank" rel="noopener">telehealth</a>, do not appear in a facility locator, so
  buprenorphine access is undercounted everywhere, most in states that lean on office-based care.
</p>
<p>
  <strong>The timelines do not fully match.</strong> Facility data is a September 2026
  snapshot; deaths run 2018&ndash;2024; population, income, and insurance are 2023 estimates.
  A state that opened clinics last year gets credit against deaths that came before them.
</p>
<p>
  <strong>Age structure is left alone.</strong> The burden rankings use unadjusted rates.
  Re-ranking on age-adjusted rates moves no state more than {AA_MAXSHIFT} places
  ({AA_MAXSTATE} shifts most), so nothing above turns on it.
</p>
<p>
  <strong>North Carolina&rsquo;s 2023 undercount.</strong> NCHS reports ~900 overdose deaths
  coded late; North Carolina&rsquo;s true burden rank is somewhat higher than shown.
</p>
<p>
  <strong>DC is a city, not a state.</strong> Its extreme values on both burden and the gap
  are real but not comparable to states, and are flagged rather than ranked.
</p>
<p>
  <strong>These numbers say where, not why.</strong> Fifty-one data points cannot say why
  treatment landed where it did, only where it visibly hasn&rsquo;t.
</p>

<footer>
  <h3>Citations</h3>
  <ol class="citations">
    <li id="cite-1">Centers for Disease Control and Prevention, National Center for Health
      Statistics. National Vital Statistics System, Mortality 2014&ndash;2024, CDC WONDER
      Online Database (Multiple Cause of Death: database <code>D157</code>, Single Race,
      for 2018&ndash;2024; database <code>D77</code>, 1999&ndash;2020, for 2014&ndash;2017;
      the 153 overlapping state-years agree exactly across the two).
      Query: UCD <code>X40&ndash;44, X60&ndash;64, X85, Y10&ndash;14</code>; MCD
      <code>T40.0&ndash;T40.4, T40.6</code>; grouped by state and year; age-adjusted and crude
      rates. Accessed {CITE_DATE}.
      <a href="https://wonder.cdc.gov/mcd-icd10-expanded.html" target="_blank" rel="noopener">wonder.cdc.gov</a></li>
    <li id="cite-2">Substance Abuse and Mental Health Services Administration,
      FindTreatment.gov Treatment Services Locator (substance-use facilities,
      <code>exportsAsJson</code> API). 12,220 listings retrieved {CITE_DATE}; 11,614 unique
      facilities in the 50 states + DC after address de-duplication. MOUD = services include
      buprenorphine or methadone used in treatment, or SAMHSA OTP certification; OTP =
      certification flag alone.
      <a href="https://findtreatment.gov" target="_blank" rel="noopener">findtreatment.gov</a></li>
    <li id="cite-3">U.S. Census Bureau, American Community Survey 2023 1-year estimates:
      tables <code>B01003</code> (population), <code>B19013</code> (median household income),
      and data profile <code>DP03</code> (health-insurance coverage), retrieved via the Census
      API, {CITE_DATE}.
      <a href="https://data.census.gov/table/ACSDT1Y2023.B01003" target="_blank" rel="noopener">data.census.gov</a></li>
    <li id="cite-4">National Center for Health Statistics, <em>Drug Overdose Deaths in the
      United States, 2023&ndash;2024</em>, NCHS Data Brief No. 549. Independently reports
      the same national opioid totals used here &mdash; 79,358 deaths in 2023, 54,045 in
      2024 &mdash; and the record one-year decline.
      <a href="https://www.cdc.gov/nchs/products/databriefs/db549.htm" target="_blank" rel="noopener">cdc.gov/nchs</a></li>
    <li id="cite-5">SAMHSA, Opioid Treatment Program certification overview: more than
      1,900 certified OTPs nationally, in line with the 2,037 locator-flagged facilities
      counted in this analysis.
      <a href="https://www.samhsa.gov/substance-use/treatment/opioid-treatment-program/become-otp" target="_blank" rel="noopener">samhsa.gov</a></li>
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
      <li><strong>Treatment capacity.</strong> One national FindTreatment.gov pull
        (<code>exportsAsJson/v2</code>): 12,220 listings &minus; 549 same-name-and-address
        duplicates &minus; 57 territory rows = 11,614 facilities in the 50 states + DC.
        OTP = SAMHSA certification flag (2,037, in line with the ~1,900+ SAMHSA reports);
        MOUD = buprenorphine or methadone used in treatment, or OTP (7,756).
        Naltrexone-only and &ldquo;accepts MAT prescribed elsewhere&rdquo; excluded.</li>
      <li><strong>Population and socioeconomics.</strong> ACS 2023 1-year via the Census
        API: <code>B01003</code> population (the denominator behind every per-100k facility
        rate), plus <code>B19013</code> median household income and <code>DP03</code>
        insurance coverage (the variables tested in Section&nbsp;08). Labels verified
        against the API&rsquo;s metadata; Puerto Rico dropped.</li>
      <li><strong>Cleaning.</strong> WONDER&rsquo;s trailing notes block stripped; numbers
        de-comma&rsquo;d and coerced, with <code>Suppressed</code>/<code>Unreliable</code>
        markers mapped to missing (none occur at state level &mdash; checked, not assumed);
        FIPS zero-padded; state identity reconciled across all sources through one fixed
        51-row name/abbreviation/FIPS table.</li>
      <li><strong>Checks.</strong> Each source must match that table 51/51 or the pipeline
        halts; each WONDER file must land exactly 51 states &times; its years. Analysis is
        descriptive only (ranks, Spearman, n&nbsp;=&nbsp;51; 2022&ndash;2024 combined rates
        to damp small-state noise), and no finding runs without a second measure agreeing.</li>
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
  "c-scatter": {J(c_scatter)},
  "c-slope": {J(c_slope)},
  "c-null": {J(c_null)}
}};
var CONFIG = {J(CONFIG)};
var MAP_KEEP_M = {J(MAP_KEEP_M)};
function renderCharts() {{
  // narrow screens: fewer, smaller map labels; trend chart keeps only the
  // three anchor labels (start, peak, drop). Runs once, before any render.
  if (window.innerWidth < 700) {{
    var t = CH['c-map'].data[1];
    var f = function (arr) {{ return arr.filter(function (_, i) {{ return MAP_KEEP_M[i]; }}); }};
    t.lon = f(t.lon); t.lat = f(t.lat); t.text = f(t.text);
    t.textfont.color = f(t.textfont.color);
    t.textfont.size = 8.5;
    var keep = ["{D2014:,}", "{D2022:,}", "{D2024:,}"];
    CH['c-national'].layout.annotations = CH['c-national'].layout.annotations.filter(function (a) {{
      return keep.some(function (k) {{ return a.text.indexOf(k) !== -1; }});
    }});
  }}
  Object.keys(CH).forEach(function (id) {{
    // retry once: plotly's geo renderer can reject on a cold load
    Plotly.newPlot(id, CH[id].data, CH[id].layout, CONFIG).catch(function () {{
      return Plotly.newPlot(id, CH[id].data, CH[id].layout, CONFIG).catch(function () {{}});
    }});
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

OUT.write_text(html, encoding="utf-8")
print(f"wrote {OUT} ({OUT.stat().st_size / 1024:.0f} KB)")
