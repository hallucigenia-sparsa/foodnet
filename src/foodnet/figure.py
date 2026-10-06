"""The consumed and produced matrices as one image, in the style of the hand-checked reference matrices.

Karoline, 2026-10-04: "as the first result, generate and show a downloadable image of the 2 matrices". An SVG drawn with the standard library, so the tool keeps no runtime dependencies; it opens in
a browser and imports into Inkscape, Illustrator and PowerPoint.

The cell states are the figure's, and there are exactly these:

  * **a number on a gray**: a change beyond the detection limit in the value medium, in mM, two decimals
    below 1 and one above; the gray darkens with the amount on a square-root scale, since the largest
    turnover exceeds the smallest by two orders of magnitude;
  * **white**: measured in the value medium, and no change beyond the limit in this direction;
  * **pale orange**: never assayed for this taxon in the value medium, a color outside the gray scale, so the
    absence of a measurement never reads as a measured zero;
  * **an open circle**: the change was seen in another medium (presence only), drawn on white or orange.

With booleans, a measured 1 is dark gray with no number. With the phase choice Both, each metabolite has a
column per phase.
"""
from __future__ import annotations

import html
import math

from . import brand
from .matrix import columns, pair_rows, taxa_rows

CELL_W, CELL_H = 44, 24
NOT_ASSAYED = "#F6D9B8"
GAP = 36              # between the two panels
MAX_SIDE_BY_SIDE = 1000   # px; wider, the panels stack
FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif"


def _text_width(text: str, size: float) -> float:
    return len(text) * size * 0.56


def _taxon_label(name: str) -> str:
    """Genus and species in italics, the strain designation upright, as the figure writes names."""
    words = name.split()
    if len(words) < 2 or name.startswith("genus:"):
        return f"<tspan font-style=\"italic\">{html.escape(name)}</tspan>"
    italic = " ".join(words[:2])
    rest = " ".join(words[2:])
    return (f"<tspan font-style=\"italic\">{html.escape(italic)}</tspan>"
            + (f" {html.escape(rest)}" if rest else ""))


def _wrap(text: str, width: float, size: float) -> list:
    """Lines of `text` that fit `width` pixels at font `size` (estimated, as `_text_width`)."""
    lines, line = [], ""
    for word in text.split():
        trial = f"{line} {word}".strip()
        if line and _text_width(trial, size) > width:
            lines.append(line)
            line = word
        else:
            line = trial
    return lines + ([line] if line else [])


def _gray(share: float) -> str:
    level = round(255 - share * (255 - 45))
    return f"#{level:02x}{level:02x}{level:02x}"


def _number(v: float) -> str:
    return f"{v:.1f}" if v >= 1 else f"{v:.2f}"


def matrices_svg(result: dict) -> str:
    """Both matrices, side by side or stacked when that would be too wide, with their key, as an SVG."""
    net = result["network"]
    pair = pair_rows(net, result)
    taxa = taxa_rows(net, result)
    cols = columns(net, result)
    booleans = bool(result["settings"].get("booleans"))
    limit = result["settings"].get("detection_limit", 0.2)
    values = [v for d in ("consumed", "produced") for row in pair[d] for v in row if isinstance(v, (int, float))]
    vmax = max(values, default=0) or 1.0
    names = pair["taxa"]
    label_w = max((_text_width(n, 12) for n in names), default=60) + 12
    # short column labels: the metabolite, its phase only when the image holds several phases, and a star
    # for the second time window, which the caption explains; the CSV files keep the full labels
    second = set((result.get("second_window") or {}).get("metabolites") or ())
    shown_phases = {ph for m, ph, _ in cols if m.id not in second}
    short = []
    for (m, ph, _), full in zip(cols, pair["columns"], strict=True):
        base = full.split(" (")[0] if " (" in full and not full.startswith("(") else full
        if m.id in second:
            short.append(f"{base} *")
        elif len(shown_phases) > 1:
            short.append(f"{base} ({ph})")
        else:
            short.append(base)
    col_label_h = max((_text_width(c, 11) for c in short), default=40) * 0.72 + 24
    # the slanted labels lean right beyond the last column
    right = max((_text_width(c, 11) * 0.77 - (len(cols) - 1 - j) * CELL_W for j, c in enumerate(short)),
                default=0)
    right = max(20, right)
    panel_w = len(cols) * CELL_W
    grid_h = len(taxa) * CELL_H
    # side by side, as in the figure, while that stays readable; one above the other when it would not
    stacked = label_w + 2 * panel_w + GAP > MAX_SIDE_BY_SIDE
    block_h = 34 + col_label_h + grid_h
    if stacked:
        width = label_w + panel_w + right
        origins = [(label_w, 0), (label_w, block_h + 18)]
    else:
        width = label_w + 2 * panel_w + GAP + right
        origins = [(label_w, 0), (label_w + panel_w + GAP, 0)]
    out = []
    unit = "" if booleans else " (mM)"
    for p, direction in enumerate(("consumed", "produced")):
        x0, y0 = origins[p]
        top = y0 + 34 + col_label_h
        if p == 0 or stacked:
            for i, name in enumerate(names):
                y = top + i * CELL_H + CELL_H / 2 + 4
                out.append(f'<text x="{label_w - 8:.0f}" y="{y:.1f}" font-size="12" text-anchor="end" '
                           f'fill="{brand.INK}">{_taxon_label(name)}</text>')
        out.append(f'<text x="{x0:.0f}" y="{y0 + 20:.0f}" font-size="14" font-weight="600" fill="{brand.INK}">'
                   f'{direction.capitalize()}{unit}</text>')
        for j, label in enumerate(short):
            cx = x0 + j * CELL_W + CELL_W / 2
            out.append(f'<text transform="translate({cx:.1f},{top - 6:.1f}) rotate(-40)" font-size="11" '
                       f'fill="{brand.INK}">{html.escape(label)}</text>')
        for i, t in enumerate(taxa):
            for j, (m, ph, _) in enumerate(cols):
                x, y = x0 + j * CELL_W, top + i * CELL_H
                v, evidence = pair[direction][i][j], pair[f"evidence_{direction}"][i][j]
                if evidence == "measured":
                    share = 1.0 if booleans else math.sqrt(min(1.0, v / vmax))
                    fill = _gray(0.75 if booleans else share)
                    out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{fill}"/>')
                    if not booleans:
                        ink = "#ffffff" if share > 0.55 else brand.INK
                        out.append(f'<text x="{x + CELL_W / 2:.1f}" y="{y + CELL_H / 2 + 4:.1f}" font-size="10.5" '
                                   f'text-anchor="middle" fill="{ink}">{_number(v)}</text>')
                    continue
                cell = result["cells"].get((t.id, m.id, ph))
                assayed = cell is not None and cell["n"]
                fill = "#ffffff" if assayed else NOT_ASSAYED
                out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{fill}"/>')
                seen = (result["presence"].get((t.id, m.id, ph)) or {}).get(direction)
                if seen:
                    out.append(f'<circle cx="{x + CELL_W / 2:.1f}" cy="{y + CELL_H / 2:.1f}" r="4.2" fill="none" '
                               f'stroke="{brand.MUTED}" stroke-width="1.2"/>')
        # the white grid between cells, as in the figure
        for j in range(len(cols) + 1):
            x = x0 + j * CELL_W
            out.append(f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + grid_h}" stroke="#ffffff" stroke-width="1.5"/>')
        for i in range(len(taxa) + 1):
            y = top + i * CELL_H
            out.append(f'<line x1="{x0}" y1="{y}" x2="{x0 + panel_w}" y2="{y}" stroke="#ffffff" stroke-width="1.5"/>')
        out.append(f'<rect x="{x0}" y="{top}" width="{panel_w}" height="{grid_h}" fill="none" '
                   f'stroke="{brand.LINE}" stroke-width="1"/>')
    top = origins[1][1] + 34 + col_label_h if stacked else 34 + col_label_h
    # the key
    ky = top + grid_h + 26
    rule = result["value_rule"]
    medium = " / ".join(rule["media"]) if rule["rule"] != "all" else "every medium"
    items = [(_gray(0.5), "", f"measured in {medium}" + ("" if booleans else ", mM (square-root shading)")),
             ("#ffffff", "", f"measured, no change beyond {limit:g} mM"),
             (NOT_ASSAYED, "", "not assayed"),
             ("#ffffff", "circle", "seen only in another medium")]
    x = label_w
    for fill, mark, text in items:
        item_w = 22 + _text_width(text, 11) + 22
        if x > label_w and x + item_w > width:
            x, ky = label_w, ky + 20
        out.append(f'<rect x="{x:.0f}" y="{ky - 10}" width="16" height="12" fill="{fill}" stroke="{brand.LINE}"/>')
        if mark:
            out.append(f'<circle cx="{x + 8:.0f}" cy="{ky - 4}" r="3.6" fill="none" stroke="{brand.MUTED}" '
                       'stroke-width="1.2"/>')
        out.append(f'<text x="{x + 22:.0f}" y="{ky}" font-size="11" fill="{brand.MUTED}">{html.escape(text)}</text>')
        x += item_w
    phase = net.meta.get("phase")
    window = net.meta.get("window")
    what = (f"window {window[0]:g} to {window[1]:g} h" if window else
            {"both": "both growth phases", "exponential": "exponential phase",
             "stationary": "stationary phase"}.get(phase, phase or ""))
    second_info = result.get("second_window") or {}
    if second_info.get("metabolites"):
        what += f" (* {second_info['label']})"
    caption = (f"Net change over the {what}, mean over replicates. foodnet {net.meta.get('tool_version', '')}, "
               f"{net.meta.get('derived_on', '')}; mGrowthDB studies {', '.join(sorted(net.studies))}.")
    for line in _wrap(caption, width - label_w - 10, 11):
        ky += 18
        out.append(f'<text x="{label_w:.0f}" y="{ky}" font-size="11" fill="{brand.MUTED}">{html.escape(line)}</text>')
    height = ky + 20
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" '
            f'width="{width:.0f}" height="{height:.0f}" font-family="{FONT}" role="img" '
            f'aria-label="Consumed and produced matrices">',
            f'<rect width="{width:.0f}" height="{height:.0f}" fill="#ffffff"/>']
    return "\n".join(head + out + ["</svg>"]) + "\n"
