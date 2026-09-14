"""Name the chart graphics that Plotly leaves unnamed.

Plotly builds its SVG in the browser and gives it no `<title>`, and Streamlit exposes no
argument to add one, so every chart arrived as an unnamed graphic. A screen reader meeting one
either skips it or reads its loose text nodes, which are axis ticks and bare numbers in no
useful order.

The fix pairs each chart with the sentence already written for it. `chart()` refuses to draw a
figure without a plain-language caption, so the text exists in the DOM immediately below every
chart; this attaches it. The chart container becomes `role="img"` with that sentence as its
accessible name, and the SVG underneath is marked `aria-hidden` so its fragments are not read
on top of the name.

`st.markdown` strips `<script>`, so this runs through `components.html`, whose iframe is
same-origin and can reach `window.parent.document`. A MutationObserver reapplies it because
Streamlit replaces the DOM on every rerun.

THE COST, STATED. That iframe makes Chromium log nine warnings per page: eight deprecated
feature-policy names in the `allow` attribute Streamlit writes, and one standard notice about a
sandbox carrying both allow-scripts and allow-same-origin. They come from Streamlit's iframe
configuration, not from this file, and none is an error. The trade is nine benign framework
warnings against two hundred and five unnamed graphics, and the graphics matter more. If
Streamlit ever exposes a chart-accessibility hook, this module should go.
"""

from __future__ import annotations

import streamlit.components.v1 as components

_SCRIPT = """
<script>
(function () {
  const doc = window.parent.document;

  function captionFor(plot) {
    // The caption chart() writes sits in the block immediately after the chart.
    let node = plot.closest('[data-testid="stElementContainer"]') || plot.parentElement;
    for (let hops = 0; node && hops < 4; hops++) {
      node = node.nextElementSibling;
      if (!node) break;
      const p = node.querySelector && node.querySelector('p.plain');
      if (p) return p.textContent.trim();
    }
    return null;
  }

  function patch() {
    doc.querySelectorAll('.js-plotly-plot').forEach(function (plot) {
      if (plot.dataset.a11yDone === '1') return;
      const text = captionFor(plot);
      if (!text) return;
      plot.setAttribute('role', 'img');
      plot.setAttribute('aria-label', text);
      plot.querySelectorAll('svg').forEach(function (svg) {
        svg.setAttribute('aria-hidden', 'true');
        svg.setAttribute('focusable', 'false');
      });
      plot.dataset.a11yDone = '1';
    });

    // Framework icons. Streamlit's tooltip triggers, its chrome glyphs and one offscreen
    // measurement node all arrive as unnamed SVGs. They carry no information a reader needs:
    // a tooltip's meaning is its text, which is exposed separately. Anything this project
    // draws itself lives inside .fig and already has a <title>, so it is excluded.
    doc.querySelectorAll('svg:not([aria-hidden]):not([aria-label])').forEach(function (svg) {
      if (svg.closest('.fig')) return;
      if (svg.querySelector('title')) return;
      svg.setAttribute('aria-hidden', 'true');
      svg.setAttribute('focusable', 'false');
    });
  }

  patch();
  new MutationObserver(function () { window.requestAnimationFrame(patch); })
    .observe(doc.body, { childList: true, subtree: true });
})();
</script>
"""


def attach_chart_labels() -> None:
    """Call once per run, after the charts are written."""
    components.html(_SCRIPT, height=0, width=0)
