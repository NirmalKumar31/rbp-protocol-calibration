"""Colour and layout constants, and the Plotly template every figure uses.

The three categorical colours were checked, not chosen by eye: pairwise OKLab separation under
normal vision and under simulated protanopia, deuteranopia and tritanopia (Machado 2009
matrices at full severity), plus contrast against the chart surface. Worst normal-vision
separation is dE 20.4 against a floor of 15, worst CVD separation is dE 14.4 against a target
of 8, and the weakest contrast is 3.19:1. An earlier blue/amber/teal set failed at dE 6.4
because blue and teal collapse together under tritanopia.

The surface is pinned light in `.streamlit/config.toml` for the same reason: the palette was
validated against one surface, so that is the surface it is shown on.
"""

SURFACE = "#FFFFFF"
PANEL = "#F8FAFC"
GRID = "#E5E9F0"
INK = "#0F172A"
INK_MUTED = "#5A6B84"

# Fixed order. A series keeps its colour when a filter removes its neighbours.
BLUE, AMBER, SLATE = "#2563EB", "#D97706", "#475569"
CATEGORICAL = (BLUE, AMBER, SLATE)

PROTOCOL_COLOUR = {"gc": BLUE, "dn": AMBER, "neg2": SLATE}
MODEL_COLOUR = {"kmer": BLUE, "cnn": AMBER, "splicebert": SLATE}
ESTIMATOR_COLOUR = {"two_stage": AMBER, "cross_fitted": BLUE}

PROTOCOL_LABEL = {
    "gc": "GC-matched",
    "dn": "Dinucleotide-matched",
    "neg2": "Bias-aware",
}
MODEL_LABEL = {
    "kmer": "4-mer logistic regression",
    "cnn": "DeepBind-style CNN",
    "splicebert": "SpliceBERT",
}

FONT = "ui-sans-serif, -apple-system, 'Segoe UI', Inter, Helvetica, Arial, sans-serif"


def plotly_template():
    """One template, applied to every figure, so nothing is styled ad hoc."""
    return {
        "layout": {
            "font": {"family": FONT, "size": 13, "color": INK},
            "paper_bgcolor": SURFACE,
            "plot_bgcolor": SURFACE,
            "colorway": list(CATEGORICAL),
            "margin": {"l": 64, "r": 24, "t": 48, "b": 56},
            "hoverlabel": {
                "bgcolor": "#FFFFFF",
                "bordercolor": GRID,
                "font": {"family": FONT, "size": 12, "color": INK},
            },
            "xaxis": {
                "gridcolor": GRID, "zeroline": False, "linecolor": GRID,
                "ticks": "outside", "tickcolor": GRID, "tickfont": {"color": INK_MUTED},
                "title": {"font": {"color": INK_MUTED, "size": 12}},
            },
            "yaxis": {
                "gridcolor": GRID, "zeroline": False, "linecolor": GRID,
                "ticks": "outside", "tickcolor": GRID, "tickfont": {"color": INK_MUTED},
                "title": {"font": {"color": INK_MUTED, "size": 12}},
            },
            "legend": {
                "orientation": "h", "yanchor": "bottom", "y": 1.02,
                "xanchor": "left", "x": 0,
                "font": {"size": 12, "color": INK_MUTED},
                "title": {"text": ""},
            },
            "title": {"font": {"size": 15, "color": INK}, "x": 0, "xanchor": "left"},
        }
    }


CSS = """
<style>
  .stApp { background: #FFFFFF; }
  section[data-testid="stSidebar"] { background: #F8FAFC; border-right: 1px solid #E5E9F0; }
  h1, h2, h3 { letter-spacing: -0.011em; color: #0F172A; font-weight: 600; }
  h1 { font-size: 1.55rem !important; }
  h2 { font-size: 1.16rem !important; margin-top: 1.5rem; }
  h3 { font-size: 0.98rem !important; }
  .block-container { padding-top: 2.2rem; max-width: 1180px; }
  [data-testid="stMetricValue"] { font-size: 1.5rem; font-weight: 600;
      font-variant-numeric: tabular-nums; color: #0F172A; }
  [data-testid="stMetricLabel"] { color: #5A6B84; font-size: 0.78rem;
      text-transform: uppercase; letter-spacing: 0.045em; }
  .src { color: #5A6B84; font-size: 0.74rem; font-family: ui-monospace, SFMono-Regular,
      Menlo, monospace; margin: -0.4rem 0 1.4rem 0; }
  .note { background: #F8FAFC; border-left: 2px solid #475569; padding: 0.7rem 0.95rem;
      margin: 0.9rem 0; font-size: 0.87rem; color: #0F172A; }
  .warn { background: #FFFBEB; border-left: 2px solid #D97706; padding: 0.7rem 0.95rem;
      margin: 0.9rem 0; font-size: 0.87rem; color: #0F172A; }
  hr { margin: 1.6rem 0; border-color: #E5E9F0; }
  .stTabs [data-baseweb="tab"] { font-size: 0.9rem; }
</style>
"""
