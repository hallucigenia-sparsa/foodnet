"""The legend: what each node shape, arc color, dash and shade in a foodnet network means.

It is drawn from the same constants the Cytoscape style uses (`brand`), and `tests/test_legend.py` requires it
to name every direction, evidence, phase and caution of the model, so a new value cannot ship unexplained
(grownet's rule).
"""
from __future__ import annotations

import html

from . import brand
from .model import CAUTIONS

WIDTH = 900
ROW = 34

CAUTION_TEXT = {
    "single_replicate": "one replicate (by default inconclusive; decided only when judging by the mean alone)",
    "short_record": "metabolites recorded for less than 24 h",
    "window_beyond_data": "the phase or window starts before the first or ends after the last sample; it stands in",
    "stationary_not_reached": "the culture was still growing at its last sample: no stationary phase",
    "conflict": "experiments (or taxa, merged to genus) disagree: inconclusive, no arc",
    "inconclusive": "the replicates pin down neither a change nor no change: NA, no arc (cautions.csv)",
    "not_detected_in_value_medium": "seen in another medium, and measured without a change in the value medium",
    "phase_from_other_replicates": "no growth curve in this replicate; its experiment's median boundary is used",
    "coarse_sampling": "the phase boundary rests on fewer than three growth samples",
    "boundaries_differ": "the replicates end exponential growth further apart than a sampling interval",
    "no_variance": "identical replicates (rounding, or one series twice): no test",
    "amounts_differ": "the experiments agree in direction, with amounts more than twofold apart",
    "experiment_left_out": "an inconclusive experiment that does not contradict the others is left out",
    "still_changing": "the compound kept changing after the growth rate fell, before the 90% rule would have ended "
                      "growth (a slow-down, or a second substrate)",
    "growth_rate_boundary": "growth ended where its rate fell, more than a sample before 90% of its maximum",
    "whole_run": "no end of growth was found, so the change is over the whole run, not a phase",
    "not_grown": "the culture did not grow (neither 1.5-fold nor, in OD, by 0.1) nor metabolize: no value, no arc",
    "growth_unclear": "its curve shows no growth (often an OD read without its blank), but it metabolized",
    "growth_unknown": "no growth curve at all: whole-run change, growth never checked",
    "within_scatter": "a change within its window's own scatter: it can stop a call, not make one",
    "within_evaporation": "no growth: a rise evaporation could explain (a set share): it can stop a call, not make one",
    "start_differs": "the replicates' metabolite samples start apart, so they cover different stretches of the phase",
    "pair_decided": "decided on two replicates (a change, or no change); a pair errs more readily than three",
}


def _arrow(y: float, color: str, dash: str = "", opacity: float = 1.0, width: float = 3.0) -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<g opacity="{opacity}"><line x1="40" y1="{y}" x2="150" y2="{y}" stroke="{color}" '
            f'stroke-width="{width}"{dash_attr}/>'
            f'<path d="M 162 {y} L 148 {y - 7} L 148 {y + 7} z" fill="{color}"/></g>')


def _text(x: float, y: float, text: str, weight: str = "normal") -> str:
    return (f'<text x="{x}" y="{y + 5}" font-size="15" font-weight="{weight}" fill="{brand.INK}" '
            f'font-family="system-ui, sans-serif">{html.escape(text)}</text>')


def legend_svg() -> str:
    rows = []
    y = 30
    rows.append(_text(20, y, "Nodes", "600"))
    y += ROW
    rows.append(f'<circle cx="95" cy="{y}" r="12" fill="{brand.GENUS_COLORS[0]}" stroke="{brand.MUTED}"/>')
    rows.append(_text(190, y, "taxon (a strain, or a genus after merging), colored by genus"))
    y += ROW
    rows.append(f'<rect x="83" y="{y - 12}" width="24" height="24" rx="5" fill="{brand.METABOLITE_NODE}" '
                f'stroke="{brand.MUTED}"/>')
    rows.append(_text(190, y, "metabolite, by ChEBI id"))
    y += ROW + 6
    rows.append(_text(20, y, "Arcs (width: the amount in mM)", "600"))
    y += ROW
    rows.append(_arrow(y, brand.PRODUCED))
    rows.append(_text(190, y, "produced: from the taxon to the metabolite, which rose in its monoculture"))
    y += ROW
    rows.append(_arrow(y, brand.CONSUMED))
    rows.append(_text(190, y, "consumed: from the metabolite to the taxon, the metabolite fell"))
    y += ROW
    rows.append(_arrow(y, brand.MUTED, "12 7"))
    rows.append(_text(190, y, "presence_only: seen in another medium, so its direction counts and not its size"))
    y += ROW
    rows.append(_arrow(y, brand.MUTED, "2 5"))
    rows.append(_text(190, y, "measured in one replicate (single_replicate; only when judging by the mean alone)"))
    y += ROW
    rows.append(_arrow(y, brand.MUTED))
    rows.append(_text(190, y, "measured: a mean change beyond the detection limit, in the medium values come from"))
    y += ROW
    rows.append(_arrow(y, brand.PRODUCED, opacity=0.45))
    rows.append(_text(190, y, "stationary phase (fainter); exponential phase and window full; whole run between"))
    y += ROW + 6
    rows.append(_text(20, y, "Cautions (an arc column; unless one says no arc, the arc is shown)", "600"))
    for flag in CAUTIONS:
        y += 26
        rows.append(_text(40, y, f"{flag}: {CAUTION_TEXT[flag]}"))
    height = y + 30
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {height}" width="{WIDTH}" '
            f'height="{height}" role="img" aria-label="foodnet legend">'
            f'<rect width="{WIDTH}" height="{height}" fill="#ffffff"/>' + "".join(rows) + "</svg>\n")
