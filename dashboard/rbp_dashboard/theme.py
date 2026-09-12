"""Dark editorial theme, and the Plotly template every figure uses.

Colour was measured, not chosen. Three data colours, checked for pairwise OKLab separation
under normal vision and under simulated protanopia, deuteranopia and tritanopia (Machado 2009
matrices at full severity), and for contrast against the surface. Worst normal-vision separation
is dE 18.8 against a floor of 15, worst CVD separation is dE 17.6 against a target of 8, and the
weakest contrast is 3.26:1.

Separation is carried by LIGHTNESS, not saturation, and that is the whole reason this palette
looks calm. Four desaturated sets were measured first and all four failed: muting a hue removes
chroma, and chroma is what CVD separation depends on, so a soft blue and a soft clay collapse
into each other under deuteranopia. Spreading the three across a wide lightness range instead
buys separation that survives every simulation while keeping every colour muted.

NUCLEOTIDE is decorative and must never encode data. Its four colours fail CVD separation
against each other, which is acceptable for a hairline motif where nothing depends on telling
them apart. `tests/test_schemas.py` asserts no figure references it.
"""

SURFACE = "#10141C"
PANEL = "#151A24"
RAISED = "#1A212D"
GRID = "#1E2530"
BORDER = "#252D3A"
BORDER_SOFT = "#1C232E"
INK = "#E6EAF0"
INK_MUTED = "#93A0B3"
INK_FAINT = "#616D7E"

# Fixed order. A series keeps its colour when a filter removes its neighbours.
STEEL, BRASS, SLATE = "#A8C8E0", "#BE9752", "#55697F"
CATEGORICAL = (STEEL, BRASS, SLATE)
ACCENT = STEEL

# Decorative only. See the module docstring.
NUCLEOTIDE = {"A": "#5E7F6E", "C": "#4E6E86", "G": "#8A7449", "U": "#7E5A60"}

PROTOCOL_COLOUR = {"gc": STEEL, "dn": BRASS, "neg2": SLATE}
# Pre-mixed translucent fills. Written out rather than derived, because deriving them needs
# str.replace on the hex, and `replace` is on the read-only guard's banned list.
FILL_TINT = {
    "gc": "rgba(168,200,224,0.07)",
    "dn": "rgba(190,151,82,0.07)",
    "neg2": "rgba(85,105,127,0.09)",
}
MODEL_COLOUR = {"kmer": STEEL, "cnn": BRASS, "splicebert": SLATE}

PROTOCOL_LABEL = {"gc": "GC-matched", "dn": "Dinucleotide-matched", "neg2": "Bias-aware"}
PROTOCOL_SHORT = {"gc": "GC", "dn": "DINUC", "neg2": "BIAS-AWARE"}
MODEL_LABEL = {
    "kmer": "4-mer logistic regression",
    "cnn": "DeepBind-style CNN",
    "splicebert": "SpliceBERT",
}
MODEL_SHORT = {"kmer": "4-MER", "cnn": "CNN", "splicebert": "SPLICEBERT"}

FONT = "'Inter', ui-sans-serif, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif"
MONO = "'JetBrains Mono', ui-monospace, SFMono-Regular, Menlo, monospace"

# A CSS url() cannot contain a newline. This is built from parts so that no source line runs
# long while the URL handed to the browser stays on one piece.
_FONTS = (
    "https://fonts.googleapis.com/css2"
    "?family=Inter:wght@300;400;500;600"
    "&family=JetBrains+Mono:wght@400;500"
    "&display=swap"
)


def plotly_template():
    return {
        "layout": {
            "font": {"family": FONT, "size": 13, "color": INK_MUTED},
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)",
            "colorway": list(CATEGORICAL),
            "margin": {"l": 70, "r": 28, "t": 40, "b": 58},
            "transition": {"duration": 450, "easing": "cubic-in-out"},
            "hoverlabel": {
                "bgcolor": RAISED, "bordercolor": BORDER,
                "font": {"family": MONO, "size": 12, "color": INK},
            },
            "xaxis": {
                "gridcolor": BORDER_SOFT, "zeroline": False, "linecolor": BORDER,
                "ticks": "outside", "tickcolor": BORDER,
                "tickfont": {"color": INK_FAINT, "family": MONO, "size": 10},
                "title": {"font": {"color": INK_FAINT, "size": 11}},
            },
            "yaxis": {
                "gridcolor": BORDER_SOFT, "zeroline": False, "linecolor": BORDER,
                "ticks": "outside", "tickcolor": BORDER,
                "tickfont": {"color": INK_FAINT, "family": MONO, "size": 10},
                "title": {"font": {"color": INK_FAINT, "size": 11}},
            },
            "legend": {
                "orientation": "h", "yanchor": "bottom", "y": 1.04,
                "xanchor": "left", "x": 0,
                "font": {"size": 11, "color": INK_MUTED},
                "title": {"text": ""}, "bgcolor": "rgba(0,0,0,0)",
            },
            "title": {"font": {"size": 13, "color": INK_MUTED}, "x": 0, "xanchor": "left"},
        }
    }


_CSS_BODY = """
  /* The white band. Streamlit's chrome ships light and sits above everything. */
  header[data-testid="stHeader"] { background: transparent !important; height: 0 !important; }
  header[data-testid="stHeader"] * { display: none !important; }
  #MainMenu, footer, [data-testid="stToolbar"],
  [data-testid="stDecoration"] { display: none !important; }
  [data-testid="stAppViewContainer"] > .main { background: #10141C; }

  .stApp, .main { background: #10141C; color: #E6EAF0; }

  /* The background field. Fixed, behind everything, pointer-transparent. Two helix layers
     drifting at different rates so the depth reads without either being fast enough to
     distract, plus a vignette so the centre of the page stays quiet where the text is. */
  .bg { position: fixed; inset: 0; z-index: 0; pointer-events: none; overflow: hidden; }
  .bg-far, .bg-near { position: absolute; left: -50%; width: 200%; height: 360px;
      background-repeat: repeat-x; background-size: 1600px auto; }
  .bg-far  { top: 2%;  opacity: 0.50; animation: drift 170s linear infinite; }
  .bg-near { top: 60%; opacity: 0.40; animation: drift 112s linear infinite reverse; }

  /* Depth. Four layers at four rates; a single speed reads as a loop. */
  .bg-dust { position: absolute; inset: -10%; opacity: 0.42;
      background-repeat: repeat; background-size: 1200px 800px;
      animation: sink 240s linear infinite; }
  .bg-ring { position: absolute; top: 50%; left: 50%; width: 620px; height: 620px;
      margin: -310px 0 0 -310px; opacity: 0.34; background-repeat: no-repeat;
      background-size: contain; animation: turn 300s linear infinite; }
  @keyframes drift { from { transform: translateX(0); }
                     to   { transform: translateX(-1600px); } }
  @keyframes sink  { from { transform: translate3d(0, 0, 0); }
                     to   { transform: translate3d(-60px, -800px, 0); } }
  @keyframes turn  { from { transform: rotate(0deg); }
                     to   { transform: rotate(360deg); } }

  /* Keeps the centre of the page quiet where the text sits, without erasing the field at the
     edges. An earlier veil at 0.94 in the middle made the whole background invisible, and a
     CSS rewrite then deleted this rule outright, which left the field running unshaded under
     the body text. Both failures were silent in the browser. */
  .bg-veil { position: absolute; inset: 0; background:
      radial-gradient(1080px 600px at 50% 32%, rgba(16,20,28,0.84) 0%,
                      rgba(16,20,28,0.44) 54%, rgba(16,20,28,0.04) 100%); }

  @media (prefers-reduced-motion: reduce) {
    .bg-far, .bg-near, .bg-dust, .bg-ring { animation: none !important; }
  }
  .block-container {
      position: relative; z-index: 1; padding-top: 3rem; padding-bottom: 5rem; max-width: 1240px; }

  /* Typography leads. Colour is used sparingly and never for emphasis in text. */
  h1 { font-size: 1.72rem !important; font-weight: 300 !important; letter-spacing: -0.02em;
       color: #E6EAF0; margin-bottom: 0.35rem; }
  h2 { font-size: 0.72rem !important; font-weight: 500 !important; text-transform: uppercase;
       letter-spacing: 0.13em; color: #616D7E; margin: 2.6rem 0 0.9rem 0;
       padding-bottom: 0.5rem; border-bottom: 1px solid #1C232E; }
  h3 { font-size: 1.0rem !important; font-weight: 500 !important; color: #E6EAF0;
       margin-top: 1.6rem; letter-spacing: -0.008em; }
  h5 { font-size: 0.78rem !important; text-transform: uppercase; letter-spacing: 0.1em;
       color: #616D7E; font-weight: 500 !important; }
  p, li { color: #93A0B3; line-height: 1.62; font-weight: 400; }
  strong { color: #E6EAF0; font-weight: 500; }
  em { color: #A8C8E0; font-style: normal; }
  a { color: #A8C8E0; text-decoration: none; border-bottom: 1px solid rgba(168,200,224,0.22);
      transition: border-color 180ms ease; }
  a:hover { border-color: #A8C8E0; }

  /* Staggered entrance. Applied to the page, not to every element, so it reads as one motion. */
  @keyframes rise { from { opacity: 0; transform: translateY(9px); }
                    to   { opacity: 1; transform: none; } }
  /* Scoped deliberately. An animation on an ancestor creates a containing block, which makes
     every position:fixed descendant behave like position:absolute. The background field was
     invisible for exactly that reason: it was inside an animated subtree and got clipped. */
  .guide { animation: rise 560ms cubic-bezier(.16,.8,.3,1) both; }
  [data-testid="stMetric"] { animation: rise 520ms cubic-bezier(.16,.8,.3,1) both; }
  [data-testid="column"]:nth-child(1) [data-testid="stMetric"] { animation-delay: 40ms; }
  [data-testid="column"]:nth-child(2) [data-testid="stMetric"] { animation-delay: 110ms; }
  [data-testid="column"]:nth-child(3) [data-testid="stMetric"] { animation-delay: 180ms; }
  [data-testid="column"]:nth-child(4) [data-testid="stMetric"] { animation-delay: 250ms; }
  .stPlotlyChart { animation: rise 620ms cubic-bezier(.16,.8,.3,1) both; animation-delay: 90ms; }
  /* Names the elements that actually animate. It previously named
     `.block-container > div > div > div > div`, which no longer animates at all, so the rule
     was dead and the entrance motion was unsuppressible. */
  @media (prefers-reduced-motion: reduce) {
    .guide, [data-testid="stMetric"], .stPlotlyChart, .fig, .key-v {
        animation: none !important; }
    .rule i { animation: none !important; }
  }

  section[data-testid="stSidebar"] {
      position: relative; z-index: 2; background: #0D1116; border-right: 1px solid #1C232E; }
  section[data-testid="stSidebar"] .block-container { padding-top: 1.4rem; }
  section[data-testid="stSidebar"] * { color: #93A0B3; }
  .sidebar-title { font-family: 'JetBrains Mono', monospace; font-size: 0.66rem;
      letter-spacing: 0.19em; text-transform: uppercase; color: #616D7E; margin-bottom: 0.1rem; }

  /* Brand block */
  .brand { display: flex; align-items: center; gap: 0.65rem; padding: 0.1rem 0 0.2rem 0; }
  .brand-k { font-family: 'JetBrains Mono', monospace; font-size: 0.82rem; font-weight: 500;
      letter-spacing: 0.13em; color: #10141C; background: #A8C8E0;
      padding: 0.3rem 0.42rem; border-radius: 3px; line-height: 1; }
  .brand-t { font-size: 0.76rem; line-height: 1.28; color: #93A0B3; letter-spacing: 0.02em;
      text-transform: uppercase; font-weight: 500; }

  /* The index. One radio, restyled into a numbered navigation list. */
  section[data-testid="stSidebar"] div[role="radiogroup"] { gap: 0 !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label {
      display: flex; align-items: center; width: 100%;
      margin: 0 0 1px 0 !important; padding: 0.52rem 0.7rem 0.52rem 0.62rem;
      border-radius: 3px; border-left: 2px solid transparent; cursor: pointer;
      transition: background 170ms ease, border-color 170ms ease, padding-left 170ms ease; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:hover {
      background: #131A24; border-left-color: #3A4657; padding-left: 0.78rem; }
  /* The dot itself is redundant once the row is the target. */
  section[data-testid="stSidebar"] div[role="radiogroup"] > label > div:first-child {
      display: none !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label p {
      font-size: 0.83rem !important; color: #7E8B9C !important; margin: 0 !important;
      font-variant-numeric: tabular-nums; letter-spacing: 0.005em;
      transition: color 170ms ease; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:hover p {
      color: #C2CCDA !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:has(input:checked) {
      background: #16202C; border-left-color: #A8C8E0; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:has(input:checked) p {
      color: #E6EAF0 !important; font-weight: 500; }

  /* Section headers, injected above specific rows. A single radio cannot have real markdown
     interleaved between its options, and splitting it would let two views be selected. The
     nth-of-type numbers track the order of VIEWS in app.py and must move when a view is added. */
  section[data-testid="stSidebar"] div[role="radiogroup"] > label::before {
      font-family: 'JetBrains Mono', monospace; font-size: 0.56rem; letter-spacing: 0.15em;
      color: #4C5666; position: absolute; margin-top: -1.55rem; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(1) {
      margin-top: 1.5rem !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(1)::before {
      content: "START HERE"; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(3) {
      margin-top: 1.6rem !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(3)::before {
      content: "THE METHOD"; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(4) {
      margin-top: 1.6rem !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(4)::before {
      content: "THE FINDING"; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(8) {
      margin-top: 1.6rem !important; }
  section[data-testid="stSidebar"] div[role="radiogroup"] > label:nth-of-type(8)::before {
      content: "DOES IT HOLD UP"; }

  .navfoot { display: flex; justify-content: space-between; align-items: center;
      font-family: 'JetBrains Mono', monospace; font-size: 0.6rem; letter-spacing: 0.1em;
      text-transform: uppercase; color: #4C5666; margin: 1.4rem 0 0 0.62rem; }
  .navfoot span { color: #616D7E; }

  /* Hairline motif. Decorative: encodes nothing. */
  .rule { display: flex; gap: 2px; margin: 0.55rem 0 1.9rem 0; height: 2px; }
  .rule i { flex: 1; border-radius: 1px; animation: breathe 7s ease-in-out infinite; }
  @keyframes breathe { 0%,100% { opacity: 0.22; } 50% { opacity: 0.6; } }

  /* The reading guide. One per view, always at the top. */
  .guide { background: linear-gradient(180deg, #151A24 0%, #12171F 100%);
      border: 1px solid #252D3A; border-left: 2px solid #A8C8E0; border-radius: 3px;
      padding: 1.15rem 1.35rem; margin: 0.4rem 0 2rem 0; }
  .guide .gl { font-family: 'JetBrains Mono', monospace; font-size: 0.62rem;
      letter-spacing: 0.16em; text-transform: uppercase; color: #616D7E;
      margin-bottom: 0.55rem; }
  .guide p { color: #C2CCDA; margin: 0 0 0.7rem 0; font-size: 0.93rem; }
  .guide ol { margin: 0; padding-left: 1.15rem; }
  .guide li { color: #93A0B3; font-size: 0.875rem; margin-bottom: 0.34rem; }
  .guide li::marker { color: #616D7E; font-family: 'JetBrains Mono', monospace;
      font-size: 0.78rem; }

  [data-testid="stMetric"] { background: #151A24; border: 1px solid #252D3A;
      border-radius: 3px; padding: 0.9rem 1.05rem;
      transition: border-color 240ms ease, transform 240ms cubic-bezier(.16,.8,.3,1); }
  [data-testid="stMetric"]:hover { border-color: #3A4657; transform: translateY(-2px); }
  [data-testid="stMetricValue"] { font-family: 'JetBrains Mono', monospace;
      font-size: 1.3rem !important; font-weight: 400; color: #E6EAF0;
      font-variant-numeric: tabular-nums; letter-spacing: -0.01em; }
  [data-testid="stMetricLabel"] { color: #616D7E !important; font-size: 0.64rem !important;
      text-transform: uppercase; letter-spacing: 0.11em; font-weight: 500; }

  /* Schematic diagrams. Content, not decoration, so they get a frame and breathing room. */
  .fig { background: linear-gradient(180deg, #151A24 0%, #12171F 100%);
      border: 1px solid #1F2735; border-radius: 4px; padding: 1.1rem 1.3rem 0.9rem 1.3rem;
      margin: 0.6rem 0 0.9rem 0; overflow: hidden;
      animation: rise 660ms cubic-bezier(.16,.8,.3,1) both; }
  .fig svg { display: block; width: 100%; height: auto; }

  /* One plain sentence under every chart, for a reader outside the field. */
  .plain { color: #8E9CAF; font-size: 0.845rem; line-height: 1.55; margin: 0.15rem 0 0.3rem 0;
      padding-left: 0.75rem; border-left: 1px solid #252D3A; }
  .plain b { color: #C2CCDA; font-weight: 500; }

  .src { color: #4C5666; font-size: 0.68rem; font-family: 'JetBrains Mono', monospace;
      margin: -0.35rem 0 1.7rem 0; }

  .note, .warn, .good, .bad { border-radius: 3px; padding: 0.85rem 1.1rem; margin: 1.1rem 0;
      font-size: 0.875rem; border-left: 2px solid; background: #151A24; line-height: 1.6; }
  .note { border-color: #55697F; color: #B0BCCB; }
  .warn { border-color: #BE9752; color: #C9BCA3; }
  .good { border-color: #5E7F6E; color: #AFC3B7; }
  .bad  { border-color: #7E5A60; color: #C7AEB2; }

  .pub, .live, .static { display: inline-block; font-family: 'JetBrains Mono', monospace;
      font-size: 0.58rem; letter-spacing: 0.11em; text-transform: uppercase;
      padding: 0.2rem 0.48rem; border-radius: 2px; margin-left: 0.5rem;
      vertical-align: middle; font-weight: 500; }
  .pub    { background: rgba(168,200,224,0.09); color: #A8C8E0;
            border: 1px solid rgba(168,200,224,0.26); }
  .live   { background: rgba(190,151,82,0.10);  color: #BE9752;
            border: 1px solid rgba(190,151,82,0.3); }
  .static { background: rgba(97,109,126,0.13);  color: #7E8B9C;
            border: 1px solid rgba(97,109,126,0.3); }

  .stTabs [data-baseweb="tab-list"] { gap: 1.6rem; border-bottom: 1px solid #1C232E; }
  .stTabs [data-baseweb="tab"] { font-size: 0.83rem; color: #616D7E; padding: 0.4rem 0;
      background: transparent; transition: color 180ms ease; }
  .stTabs [data-baseweb="tab"]:hover { color: #93A0B3; }
  .stTabs [aria-selected="true"] { color: #E6EAF0 !important; }
  .stTabs [data-baseweb="tab-highlight"] { background: #A8C8E0; }

  div[data-testid="stExpander"] { border: 1px solid #1C232E; border-radius: 3px;
      background: #12171F; transition: border-color 200ms ease; }
  div[data-testid="stExpander"]:hover { border-color: #252D3A; }
  div[data-testid="stExpander"] summary { font-size: 0.85rem; color: #93A0B3; }

  /* Charts sit in the same frame the diagrams do, so a page reads as one column of panels
     rather than as figures floating on a background. */
  .stPlotlyChart { background: linear-gradient(180deg, #141A23 0%, #11161E 100%);
      border: 1px solid #1F2735; border-radius: 4px; padding: 0.9rem 0.5rem 0.4rem 0.5rem;
      transition: border-color 260ms ease; }
  .stPlotlyChart:hover { border-color: #2B3546; }

  /* Numbered section rules. The number is the reading order, which is otherwise invisible. */
  .step { display: flex; align-items: baseline; gap: 0.8rem; margin: 2.6rem 0 0.5rem 0;
      padding-bottom: 0.55rem; border-bottom: 1px solid #1C232E; }
  .step-n { font-family: 'JetBrains Mono', monospace; font-size: 0.68rem; color: #55697F;
      letter-spacing: 0.1em; }
  .step-t { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.13em;
      color: #616D7E; font-weight: 500; }

  /* A key-number band: the figures that carry the argument, set larger than body text. */
  .keys { display: flex; gap: 2.4rem; flex-wrap: wrap; margin: 0.4rem 0 1.6rem 0;
      padding: 1.05rem 1.25rem; background: #141A23; border: 1px solid #1F2735;
      border-radius: 4px; }
  .key-v { font-family: 'JetBrains Mono', monospace; font-size: 1.48rem; color: #E6EAF0;
      font-variant-numeric: tabular-nums; line-height: 1.15;
      animation: rise 620ms cubic-bezier(.16,.8,.3,1) both; }
  .key-l { font-size: 0.66rem; text-transform: uppercase; letter-spacing: 0.1em;
      color: #616D7E; margin-top: 0.22rem; }
  .key-s { color: #55697F; font-size: 0.95rem; }

  /* Pull quote, for the narrator's sharpest line on a view. */
  .pull { border-left: 2px solid #BE9752; padding: 0.2rem 0 0.2rem 1.15rem; margin: 1.5rem 0;
      font-size: 1.02rem; line-height: 1.6; color: #C2CCDA; font-weight: 300; }

  /* Narrow viewports. The desktop rhythm is built for a 1240px column; below roughly a
     tablet it needs less padding, smaller display type and single-column metric rows, or the
     guide box and the metric tiles both overflow. */
  @media (max-width: 640px) {
    .block-container { padding-top: 1.6rem; padding-left: 1rem; padding-right: 1rem; }
    h1 { font-size: 1.24rem !important; }
    .guide { padding: 0.9rem 1rem; }
    .guide p { font-size: 0.86rem; }
    .guide li { font-size: 0.82rem; }
    .keys { gap: 1.2rem; padding: 0.85rem 1rem; }
    .key-v { font-size: 1.18rem; }
    .pull { font-size: 0.94rem; padding-left: 0.9rem; }
    [data-testid="stMetricValue"] { font-size: 1.05rem !important; }
    [data-testid="column"] { min-width: 100% !important; }
    .stPlotlyChart { padding: 0.5rem 0.2rem 0.2rem 0.2rem; }
    .src { font-size: 0.62rem; }
    .bg-ring { display: none; }
  }

  .stDataFrame { border: 1px solid #1C232E; border-radius: 3px; }
  hr { border-color: #1C232E; margin: 2rem 0; }
  .stSlider [data-baseweb="slider"] div[role="slider"] { background: #A8C8E0 !important; }
  .stMultiSelect [data-baseweb="tag"] { background: #252D3A !important;
      color: #C2CCDA !important; }
</style>
"""

# Concatenated, not formatted and not .replace()d. The stylesheet is full of percent signs in
# @keyframes so a format string would need every one escaped, and `replace` is on the read-only
# guard's banned list because `Path.replace` is a destructive move and the guard cannot tell
# `str.replace` from it by reading the syntax tree. Avoiding the name costs one line.
CSS = "<style>\n  @import url('" + _FONTS + "');\n" + _CSS_BODY

RULE = (
    '<div class="rule">'
    + "".join(
        f'<i style="background:{NUCLEOTIDE[b]};animation-delay:{i * 0.21:.2f}s"></i>'
        for i, b in enumerate("ACGU" * 10)
    )
    + "</div>"
)

# ----------------------------------------------------------------- the background

def _helix_svg(width: int = 1600, height: int = 300, turns: int = 7,
               amplitude: int = 54, rungs: int = 84) -> str:
    """Two antiparallel strands and their base pairs, computed rather than hand-drawn.

    Computed because a hand-typed path cannot be re-tuned: changing the number of turns means
    retyping every control point. This takes the wave from sin() at a fixed step, so the shape
    follows the parameters.

    It is a background. It carries no data, it is not to scale, and it is not a structural claim
    about any molecule. Purely a field to sit behind the page.
    """
    import math

    mid = height / 2
    step = width / 220
    front, back = [], []
    x = 0.0
    while x <= width:
        theta = (x / width) * turns * 2 * math.pi
        front.append(f"{x:.1f},{mid + amplitude * math.sin(theta):.1f}")
        back.append(f"{x:.1f},{mid + amplitude * math.sin(theta + math.pi):.1f}")
        x += step

    bars = []
    for i in range(rungs):
        bx = width * i / (rungs - 1)
        theta = (bx / width) * turns * 2 * math.pi
        y1 = mid + amplitude * math.sin(theta)
        y2 = mid + amplitude * math.sin(theta + math.pi)
        base = "ACGU"[i % 4]
        # Rungs fade where the strands cross, which is where a real helix occludes itself.
        near = abs(y1 - y2) / (2 * amplitude)
        bars.append(
            f'<line x1="{bx:.1f}" y1="{y1:.1f}" x2="{bx:.1f}" y2="{y2:.1f}" '
            f'stroke="{NUCLEOTIDE[base]}" stroke-width="1.1" opacity="{0.1 + 0.5 * near:.2f}"/>'
        )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'width="{width}" height="{height}" fill="none">'
        f'<g>{"".join(bars)}</g>'
        f'<polyline points="{" ".join(front)}" stroke="{STEEL}" stroke-width="1.6" '
        f'opacity="0.5" stroke-linecap="round"/>'
        f'<polyline points="{" ".join(back)}" stroke="{SLATE}" stroke-width="1.6" '
        f'opacity="0.42" stroke-linecap="round"/>'
        f"</svg>"
    )


def _particle_field(width: int = 1200, height: int = 800, n: int = 46) -> str:
    """Loose nucleotide glyphs, placed on a deterministic pseudo-random lattice.

    Seeded arithmetic rather than `random` so the field is identical on every page load. A
    background that reshuffles on each rerun is a background a reader notices, which is the one
    thing it must not be.
    """
    parts = []
    for i in range(n):
        # cheap deterministic hash; no import, no seed to thread through
        x = (i * 7919 % 997) / 997 * width
        y = (i * 6551 % 991) / 991 * height
        size = 9 + (i * 13 % 7)
        base = "ACGU"[i % 4]
        dur = 16 + (i * 11 % 19)
        parts.append(
            f'<text x="{x:.0f}" y="{y:.0f}" fill="{NUCLEOTIDE[base]}" '
            f'font-family="ui-monospace,monospace" font-size="{size}" opacity="0.30">{base}'
            f'<animate attributeName="opacity" values="0.10;0.42;0.10" dur="{dur}s" '
            f'begin="{i * 0.37:.2f}s" repeatCount="indefinite"/></text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'width="{width}" height="{height}" fill="none">{"".join(parts)}</svg>')


def _ring(size: int = 620) -> str:
    """A helix seen end-on: concentric arcs with radial base pairs. Rotates very slowly.

    Gives the field a centre of gravity that the drifting strands do not, and reads as
    structural rather than decorative at the opacity it is used.
    """
    import math

    mid = size / 2
    parts = []
    for ring, radius in enumerate((mid * 0.34, mid * 0.55, mid * 0.78)):
        parts.append(
            f'<circle cx="{mid}" cy="{mid}" r="{radius:.0f}" stroke="{STEEL}" '
            f'stroke-width="1" opacity="{0.30 - ring * 0.07:.2f}" '
            f'stroke-dasharray="{18 + ring * 9} {26 + ring * 7}"/>')
    for i in range(36):
        angle = i * math.pi / 18
        inner, outer = mid * 0.34, mid * 0.78
        x1, y1 = mid + inner * math.cos(angle), mid + inner * math.sin(angle)
        x2, y2 = mid + outer * math.cos(angle), mid + outer * math.sin(angle)
        base = "ACGU"[i % 4]
        parts.append(
            f'<line x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" y2="{y2:.0f}" '
            f'stroke="{NUCLEOTIDE[base]}" stroke-width="1" opacity="0.16"/>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" '
            f'width="{size}" height="{size}" fill="none">{"".join(parts)}</svg>')


def background() -> str:
    """Four layers behind the page, moving at four rates.

    Parallax rather than a single drifting strip: a lone layer at a constant speed reads as a
    loop, and three at different rates reads as depth. The ring turns slowly enough that it is
    never the thing being looked at.

    Every layer is fixed, pointer-transparent and `aria-hidden`. None of it encodes data, none
    of it is to scale, and all of it stops under `prefers-reduced-motion`.
    """
    import base64

    def uri(svg: str) -> str:
        return "data:image/svg+xml;base64," + base64.b64encode(svg.encode()).decode()

    far = uri(_helix_svg(1600, 240, turns=11, amplitude=32, rungs=110))
    near = uri(_helix_svg(1600, 360, turns=4, amplitude=74, rungs=52))
    dust = uri(_particle_field())
    ring = uri(_ring())
    return (
        '<div class="bg" aria-hidden="true">'
        f'<div class="bg-ring" style="background-image:url({ring})"></div>'
        f'<div class="bg-dust" style="background-image:url({dust})"></div>'
        f'<div class="bg-far"  style="background-image:url({far})"></div>'
        f'<div class="bg-near" style="background-image:url({near})"></div>'
        '<div class="bg-veil"></div>'
        "</div>"
    )
