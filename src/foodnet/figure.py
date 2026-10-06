"""The consumed and produced matrices as one image, in the style of the hand-checked reference matrices.

Karoline, 2026-10-04: "as the first result, generate and show a downloadable image of the 2 matrices".
An SVG drawn with the standard library, so the tool keeps no runtime dependencies; it opens in a browser
and imports into Inkscape, Illustrator and PowerPoint.

The cell states are the figure's, and there are exactly these:

  * **a number on a gray**: a change beyond the detection limit in the value medium, in mM, two decimals
    below 1 and one above; the gray darkens with the amount on a square-root scale, since the largest
    turnover exceeds the smallest by two orders of magnitude;
  * **white**: measured in the value medium, and no change beyond the limit in this direction;
  * **pale orange**: never assayed for this taxon in the value medium, a color outside the gray scale, so the
    absence of a measurement never reads as a measured zero;
  * **a question mark on pale yellow**: assayed, but inconclusive: its replicates, or its experiments, do not
    agree on what happened (Karoline, 2026-10-06);
  * **a dash on pale green**: assayed, but its cultures gave no phase (no end of exponential growth, or no
    stationary phase reached), so there is no value in the phase asked for;
  * **a dot in the corner**: the cell rests on one replicate (single_replicate);
  * **an open circle**: the change was seen in another medium, drawn on white (the value medium measured no
    change: seen_elsewhere) or on orange (presence only).

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
# a value taken from another medium (Advanced settings: entries seen only in another medium, "value"): a
# background of its own, outside the gray scale and apart from the not-assayed orange
OTHER_MEDIUM = "#E4DCF1"
INCONCLUSIVE = "#FBEFC5"
NO_PHASE = "#DCEDE3"
MARKS = {"inconclusive": (INCONCLUSIVE, "?"), "no_phase": (NO_PHASE, "\u2013")}
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


TIP_FONT, TIP_LINE, TIP_PAD, TIP_MAX_W = 11, 15, 7, 460
BOLD = ' font-weight="600"'


def _hover_tips(tips, width: float, height: float) -> list:
    """A box per cell with its tooltip text, hidden until the mouse is over the cell.

    A native SVG title is shown only by a browser that chooses to (an embedded browser often does not; Karoline,
    2026-10-06: "the hover didn't land"), so the text is drawn as well: each box follows all the cells, so it
    lies above them, and a rule in the image's own style sheet shows it while its cell is hovered. No
    JavaScript, on the page and in the downloaded file alike. The boxes carry visibility="hidden" as an
    attribute, so a program that ignores the style sheet (a vector editor, a printer) leaves them hidden."""
    rules, boxes = [], []
    for k, (x, y, text) in enumerate(tips, start=1):
        lines = [w for line in text.split("\n") for w in _wrap(line, TIP_MAX_W - 2 * TIP_PAD, TIP_FONT)]
        w = min(TIP_MAX_W, max(_text_width(line, TIP_FONT) for line in lines) + 2 * TIP_PAD)
        h = len(lines) * TIP_LINE + 2 * TIP_PAD - 4
        bx = min(max(4.0, x + CELL_W / 2 - w / 2), max(4.0, width - w - 4))
        below = y + CELL_H + 4
        by = below if below + h <= height - 4 else max(4.0, y - h - 4)
        rows = "".join(f'<text x="{bx + TIP_PAD:.1f}" y="{by + TIP_PAD + 9 + i * TIP_LINE:.1f}" font-size="{TIP_FONT}" '
                       f'fill="{brand.INK}"{BOLD if i == 0 else ""}>{html.escape(line)}</text>'
                       for i, line in enumerate(lines))
        boxes.append(f'<g id="fn-t{k}" class="fn-tip" visibility="hidden" pointer-events="none">'
                     f'<rect x="{bx:.1f}" y="{by:.1f}" width="{w:.1f}" height="{h:.1f}" rx="4" fill="#ffffff" '
                     f'stroke="{brand.MUTED}" stroke-width="0.8"/>{rows}</g>')
        rules.append(f"#fn-c{k}:hover ~ #fn-t{k}")
    if not boxes:
        return []
    style = (f"<style>{', '.join(rules)} {{ visibility: visible; }} "
             ".fn-cell:hover > rect { stroke: #1F1D1A; stroke-width: 1.2; }</style>")
    return [style, *boxes]


def _gray(share: float) -> str:
    level = round(255 - share * (255 - 45))
    return f"#{level:02x}{level:02x}{level:02x}"


def _number(v: float) -> str:
    return f"{v:.1f}" if v >= 1 else f"{v:.2f}"


def _studies(result: dict, ids) -> str:
    """Study ids with their titles, as mGrowthDB names them."""
    studies = result["network"].studies
    return "; ".join(f"{sid} ({studies[sid].citation})" if sid in studies and studies[sid].citation not in ("", sid)
                     else sid for sid in ids)


def tooltip(result: dict, taxon, met, phase: str, direction: str, evidence: str, value) -> str:
    """The text a mouseover shows for one cell: what it is, its value and the studies behind it."""
    from .matrix import interval
    head = f"{taxon.name}, {met.name}, {direction} ({interval(result, met.id, phase)})"
    cell = result["cells"].get((taxon.id, met.id, phase))
    lines = [head]
    if evidence == "single_replicate":
        evidence = "measured" if value else "below_limit"
    if cell is not None and cell["n"]:
        limit = cell.get("limit") or result["settings"].get("detection_limit", 0.2)
        experiments = (f", {cell['n_experiments']} experiments" if (cell.get("n_experiments") or 0) > 1 else "")
        if evidence == "inconclusive":
            lines.append(f"inconclusive: mean {cell['mean']:+.3g} mM"
                         + ("" if cell["sd"] is None else f" \u00b1 {cell['sd']:.2g}")
                         + f", {cell['n']} replicate(s){experiments}; they do not agree on a change beyond "
                         f"{limit:g} mM, or on none")
        elif evidence == "measured":
            if result["settings"].get("booleans"):
                amount = "yes"
            else:
                spread = "" if cell["sd"] is None else f" \u00b1 {cell['sd']:.2g}"
                amount = f"{value:.3g}{spread} mM"
            lines.append(f"{amount}, {cell['n']} replicate(s){experiments}")
        else:
            lines.append(f"measured, no change beyond {limit:g} mM in this direction "
                         f"(mean {cell['mean']:+.3g} mM, {cell['n']} replicate(s))")
        lines.append("Studies: " + _studies(result, cell["studies"]))
        lines.append("Experiments: " + ", ".join(cell["experiments"]))
        lines.append("Medium: " + "; ".join(cell["media"]))
        if cell["cautions"]:
            lines.append("Cautions: " + ", ".join(cell["cautions"]))
        lines += cell["notes"]
    elif cell is not None:
        lines.append("assayed, but its cultures gave no phase (no end of exponential growth, or no stationary "
                     "phase reached), so no value in this phase; a time window gives one")
        lines.append("Studies: " + _studies(result, cell["studies"]))
        lines.append("Experiments: " + ", ".join(cell["experiments"]))
    else:
        lines.append("not assayed in the medium the values come from")
    for entry in (result["presence"].get((taxon.id, met.id, phase)) or {}).get(direction, []):
        lines.append(f"Seen in {entry['medium']} ({entry['mean']:+.3g} mM, {entry['n']} replicate(s)): "
                     + _studies(result, entry["studies"]))
    return "\n".join(lines)


def matrices_svg(result: dict) -> str:
    """Both matrices, side by side or stacked when that would be too wide, with their key, as an SVG."""
    net = result["network"]
    pair = pair_rows(net, result)
    taxa = taxa_rows(net, result)
    cols = columns(net, result)
    booleans = bool(result["settings"].get("booleans"))
    limit = result["settings"].get("detection_limit", 0.2)
    # the gray scale is set by the value medium's measurements only; a value from another medium is not
    # comparable with them and is drawn on its own background
    values = [v for d in ("consumed", "produced") for row, ev in zip(pair[d], pair[f"evidence_{d}"], strict=True)
              for v, e in zip(row, ev, strict=True)
              if e in ("measured", "single_replicate") and isinstance(v, (int, float))]
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
    tips = []          # (x, y, text) per cell, drawn last so a tip lies above every cell
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
                # each cell is a group with a title, which a browser shows on mouseover: what the cell is and
                # the studies behind it (Karoline, 2026-10-06: "a mouseover will show the source studies")
                text = tooltip(result, t, m, ph, direction, evidence, v)
                tips.append((x, y, text))
                out.append(f'<g id="fn-c{len(tips)}" class="fn-cell"><title>{html.escape(text)}</title>')
                lone = evidence == "single_replicate"
                if lone:
                    # one replicate: drawn as what it says, with a dot in its corner
                    evidence = "measured" if v else "below_limit"
                    out.append(f'<circle cx="{x + CELL_W - 4:.1f}" cy="{y + 4:.1f}" r="1.8" fill="{brand.MUTED}"/>')
                if evidence == "measured":
                    share = 1.0 if booleans else math.sqrt(min(1.0, v / vmax))
                    fill = _gray(0.75 if booleans else share)
                    out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{fill}"/>')
                    if not booleans:
                        ink = "#ffffff" if share > 0.55 else brand.INK
                        out.append(f'<text x="{x + CELL_W / 2:.1f}" y="{y + CELL_H / 2 + 4:.1f}" font-size="10.5" '
                                   f'text-anchor="middle" fill="{ink}">{_number(v)}</text>')
                    out.append("</g>")
                    continue
                if evidence == "presence_only" and isinstance(v, float) and not booleans:
                    out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{OTHER_MEDIUM}"/>')
                    out.append(f'<text x="{x + CELL_W / 2:.1f}" y="{y + CELL_H / 2 + 4:.1f}" font-size="10.5" '
                               f'text-anchor="middle" fill="{brand.INK}">{_number(v)}</text>')
                    out.append("</g>")
                    continue
                if evidence in MARKS:
                    fill, mark = MARKS[evidence]
                    out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{fill}"/>')
                    out.append(f'<text x="{x + CELL_W / 2:.1f}" y="{y + CELL_H / 2 + 4:.1f}" font-size="11" '
                               f'text-anchor="middle" fill="{brand.MUTED}">{mark}</text>')
                    out.append("</g>")
                    continue
                cell = result["cells"].get((t.id, m.id, ph))
                assayed = cell is not None and cell["n"]
                fill = "#ffffff" if assayed else NOT_ASSAYED
                out.append(f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" fill="{fill}"/>')
                seen = (result["presence"].get((t.id, m.id, ph)) or {}).get(direction)
                if seen:
                    out.append(f'<circle cx="{x + CELL_W / 2:.1f}" cy="{y + CELL_H / 2:.1f}" r="4.2" fill="none" '
                               f'stroke="{brand.MUTED}" stroke-width="1.2"/>')
                out.append("</g>")
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
    shown = {e for d in ("consumed", "produced") for row in pair[f"evidence_{d}"] for e in row}
    if any(isinstance(v, float) and e == "presence_only" for d in ("consumed", "produced")
           for row, ev in zip(pair[d], pair[f"evidence_{d}"], strict=True) for v, e in zip(row, ev, strict=True)):
        items[-1] = (OTHER_MEDIUM, "", "value from another medium, not comparable with the gray scale")
    if "single_replicate" in shown:
        items.append(("#ffffff", "dot", "one replicate"))
    if "inconclusive" in shown:
        items.append((INCONCLUSIVE, "?", "inconclusive: replicates or experiments disagree"))
    if "no_phase" in shown:
        items.append((NO_PHASE, "\u2013", "assayed, no phase in its cultures"))
    x = label_w
    for fill, mark, text in items:
        item_w = 22 + _text_width(text, 11) + 22
        if x > label_w and x + item_w > width:
            x, ky = label_w, ky + 20
        out.append(f'<rect x="{x:.0f}" y="{ky - 10}" width="16" height="12" fill="{fill}" stroke="{brand.LINE}"/>')
        if mark == "circle":
            out.append(f'<circle cx="{x + 8:.0f}" cy="{ky - 4}" r="3.6" fill="none" stroke="{brand.MUTED}" '
                       'stroke-width="1.2"/>')
        elif mark == "dot":
            out.append(f'<circle cx="{x + 13:.0f}" cy="{ky - 7}" r="1.8" fill="{brand.MUTED}"/>')
        elif mark:
            out.append(f'<text x="{x + 8:.0f}" y="{ky}" font-size="10" text-anchor="middle" '
                       f'fill="{brand.MUTED}">{mark}</text>')
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
    caption = (f"Net change over the {what}, mean over experiments (over replicates within one). "
               f"foodnet {net.meta.get('tool_version', '')}, "
               f"{net.meta.get('derived_on', '')}; mGrowthDB studies {', '.join(sorted(net.studies))}.")
    for line in _wrap(caption, width - label_w - 10, 11):
        ky += 18
        out.append(f'<text x="{label_w:.0f}" y="{ky}" font-size="11" fill="{brand.MUTED}">{html.escape(line)}</text>')
    height = ky + 20
    out += _hover_tips(tips, width, height)
    head = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" '
            f'width="{width:.0f}" height="{height:.0f}" font-family="{FONT}" role="img" '
            f'aria-label="Consumed and produced matrices">',
            f'<rect width="{width:.0f}" height="{height:.0f}" fill="#ffffff"/>']
    return "\n".join(head + out + ["</svg>"]) + "\n"
