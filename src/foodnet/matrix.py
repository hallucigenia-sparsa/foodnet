"""The matrices, and the parameters a consumer-resource model takes.

Karoline, 2026-10-04: "Instead of an adjacency matrix, there will be 2 matrix formats: one matrix with taxa
as rows and metabolites as columns and another with 2 matrices: 1 for consumption and the other for
production. Matrix entries and arc widths are given by the amount of metabolites produced/removed".

  * **The signed matrix** (`signed_csv`): taxa as rows, metabolites as columns, each cell the mean net
    change in mM, positive when the taxon produced the compound and negative when it consumed it.
  * **The pair** (`pair_package`): a consumed matrix and a produced matrix, both of non-negative amounts,
    so a taxon that does both in different phases shows both.

The cell vocabulary, the same in every matrix:

  * a number: a change beyond the detection limit, measured in the medium the values come from;
  * 0: measured there, and no change beyond the limit in this direction (in the pair, a compound that moved
    the other way is 0 here and a number in the other matrix);
  * NA: no value. The compound was never assayed for this taxon in the value medium; or its change was
    seen only in another medium (presence only); or it was assayed but is inconclusive (its replicates, or
    its experiments, do not agree on what happened); or it was assayed but the culture gave
    no phase (no end of exponential growth, or no stationary phase reached).

The evidence matrices of the pair say which, cell by cell (`EVIDENCE_WORDS`): measured, below_limit,
single_replicate (a number or 0 that rests on one replicate), seen_elsewhere (0 in the value medium, but another
medium showed this direction), inconclusive, no_phase, presence_only, not_assayed.

With "Report everything as booleans" a number becomes 1 (or -1 in the signed matrix) and a presence-only
cell counts as 1 too, since a boolean asks only whether it happened.

With the phase choice "Both" each metabolite has one column per phase, named "acetate (exponential)" and
"acetate (stationary)".
"""
from __future__ import annotations

import csv
import io
import json
import zipfile

from .model import FoodNetwork

NA = "NA"
EVIDENCE_WORDS = ("measured", "below_limit", "single_replicate", "seen_elsewhere", "inconclusive", "no_phase",
                  "presence_only", "not_assayed")
# v1 (0.2.0): adds each taxon's biomass change over the phase, which miaSim's yields need, and the caveats
# that make a CRM refuse to build without being told (mixed media, the stationary phase)
CRM_FORMAT = "foodnet.crm/v1"


def taxa_rows(net: FoodNetwork, result: dict | None = None) -> list:
    """The taxon nodes, by name: every taxon with data (`result["taxa_nodes"]`), not only those with arcs."""
    nodes = (result or {}).get("taxa_nodes") or net.taxa()
    return sorted(nodes, key=lambda n: (n.name.lower(), n.id))


def metabolite_nodes(net: FoodNetwork, result: dict) -> list:
    """Every metabolite measured, by name, not only those with arcs: a measured no change is information."""
    return sorted(result.get("metabolite_nodes") or net.metabolites(), key=lambda n: (n.name.lower(), n.id))


def _labels(nodes) -> list:
    """Display names, with the id added where two nodes share a name."""
    names = [n.name or n.id for n in nodes]
    return [f"{name} [{n.id}]" if names.count(name) > 1 else name for name, n in zip(names, nodes, strict=True)]


def phases_of(result: dict) -> list:
    """The phases this search's cells hold, in the order exponential, stationary, window."""
    found = {ph for (_, _, ph) in result["cells"]} | {ph for (_, _, ph) in result["presence"]}
    return [p for p in ("exponential", "stationary", "window") if p in found]


def second_window_metabolites(result: dict) -> set:
    return set((result.get("second_window") or {}).get("metabolites") or ())


def interval(result: dict, met_id: str, phase: str) -> str:
    """What a column's values were measured over, in words: the phase, the main window or the second one."""
    second = result.get("second_window") or {}
    if met_id in second_window_metabolites(result):
        return second.get("label", "second window")
    if phase == "window":
        start, end = (result["network"].meta.get("window") or [None, None])
        return f"{start:g} to {end:g} h" if start is not None else "window"
    return phase


def columns(net: FoodNetwork, result: dict, phases=None) -> list:
    """(metabolite node, phase, label) for each column: metabolites by name, one column per phase.

    A metabolite in the second time window has one column, its window (Karoline, 2026-10-05). When the
    columns are measured over more than one interval, each label says its own: "acetate (exponential)",
    "trehalose (0 h to the last sample)"."""
    phases = phases or phases_of(result) or ["exponential"]
    second = second_window_metabolites(result)
    mets = metabolite_nodes(net, result)
    labels = _labels(mets)
    cols = []
    for m, label in zip(mets, labels, strict=True):
        for ph in (["window"] if m.id in second else [p for p in phases if p != "window" or not second]
                   or ["window"]):
            cols.append((m, ph, label))
    intervals = {interval(result, m.id, ph) for m, ph, _ in cols}
    if len(intervals) == 1:
        return cols
    return [(m, ph, f"{label} ({interval(result, m.id, ph)})") for m, ph, label in cols]


def _booleans(result: dict) -> bool:
    return bool(result["settings"].get("booleans"))


PRESENCE_ENTRIES = ("na", "true", "value")
TRUE = "TRUE"


def _own_limits(s: dict) -> dict:
    from .search import compound_limits_of
    try:
        return compound_limits_of(s)
    except ValueError:
        return {}


def presence_entries(result: dict) -> str:
    """How a cell seen only in another medium is written (Karoline, 2026-10-06): "na" (the cautious default),
    "true", or "value", the change measured there."""
    choice = result["settings"].get("presence_entries", "na")
    return choice if choice in PRESENCE_ENTRIES else "na"


def presence_value(entries) -> float:
    """The change seen in other media, in mM: the mean of the media's own values, each medium counted once."""
    return sum(e["mean"] for e in entries) / len(entries)


def entry(result: dict, taxon: str, met: str, ph: str, direction: str) -> tuple:
    """(value, evidence) of one cell of the consumed or the produced matrix."""
    cell = result["cells"].get((taxon, met, ph))
    seen = (result["presence"].get((taxon, met, ph)) or {}).get(direction)
    booleans = _booleans(result)
    if cell is not None and cell["n"]:
        single = cell["n"] == 1        # decided on one replicate (Karoline, 2026-10-06): its own evidence
        if cell["direction"] == direction:
            return (1 if booleans else abs(cell["mean"])), ("single_replicate" if single else "measured")
        if cell.get("state") == "inconclusive":
            return None, "inconclusive"
        if seen and booleans:
            return 1, "presence_only"
        # measured in the value medium without a change this way: 0 there, whatever another medium showed,
        # and the evidence says when one did
        return 0, ("seen_elsewhere" if seen else "single_replicate" if single else "below_limit")
    if seen:
        if booleans:
            return 1, "presence_only"
        choice = presence_entries(result)
        if choice == "true":
            return TRUE, "presence_only"
        if choice == "value":
            return abs(presence_value(seen)), "presence_only"
        return None, "presence_only"
    if cell is not None:
        return None, "no_phase"            # assayed, but no phase (or no stationary phase) in its cultures
    return None, "not_assayed"


def signed_entry(result: dict, taxon: str, met: str, ph: str):
    """The signed matrix's cell: the mean change, 0, or None (NA)."""
    cell = result["cells"].get((taxon, met, ph))
    booleans = _booleans(result)
    if cell is not None and cell["n"]:
        if cell.get("state") == "inconclusive":
            return None
        if cell["direction"] is None:
            seen = result["presence"].get((taxon, met, ph)) or {}
            if booleans and len(seen) == 1:
                return 1 if "produced" in seen else -1
            return 0
        if booleans:
            return 1 if cell["direction"] == "produced" else -1
        return cell["mean"]
    seen = result["presence"].get((taxon, met, ph)) or {}
    if booleans and len(seen) == 1:
        return 1 if "produced" in seen else -1
    if len(seen) == 1 and not booleans:
        # one direction seen in other media: written as the setting says. TRUE carries no sign, so the
        # direction is read from the consumed and produced matrices; both directions seen stays NA
        choice = presence_entries(result)
        if choice == "true":
            return TRUE
        if choice == "value":
            return presence_value(next(iter(seen.values())))
    return None


def _number(value) -> str:
    if value is None:
        return NA
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return f"{value:.6g}"


_FORMULA = ("=", "+", "-", "@", "\t", "\r")


def _inert(cell):
    """A text cell a spreadsheet would run as a formula (names come from the database as deposited), with a
    leading apostrophe; numbers, negative ones included, stay as they are (OWASP, CSV injection)."""
    if not isinstance(cell, str) or not cell.startswith(_FORMULA):
        return cell
    try:
        float(cell)
        return cell
    except ValueError:
        return "'" + cell


def _csv(header, rows) -> str:
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow([_inert(c) for c in header])
    writer.writerows([_inert(c) for c in row] for row in rows)
    return out.getvalue()


def signed_rows(net: FoodNetwork, result: dict, phases=None) -> tuple:
    taxa = taxa_rows(net, result)
    cols = columns(net, result, phases)
    rows = [[signed_entry(result, t.id, m.id, ph) for m, ph, _ in cols] for t in taxa]
    return _labels(taxa), [label for _, _, label in cols], rows


def signed_csv(net: FoodNetwork, result: dict) -> str:
    """The taxa by metabolites matrix as CSV: produced positive, consumed negative, NA for no value."""
    names, header, rows = signed_rows(net, result)
    return _csv([corner(result), *header], [[n, *map(_number, r)] for n, r in zip(names, rows, strict=True)])


def pair_rows(net: FoodNetwork, result: dict, phases=None) -> dict:
    """{"taxa", "columns", "consumed", "produced", "evidence_consumed", "evidence_produced"}."""
    taxa = taxa_rows(net, result)
    cols = columns(net, result, phases)
    out = {"taxa": _labels(taxa), "columns": [label for _, _, label in cols]}
    for direction in ("consumed", "produced"):
        values, evidence = [], []
        for t in taxa:
            row = [entry(result, t.id, m.id, ph, direction) for m, ph, _ in cols]
            values.append([v for v, _ in row])
            evidence.append([e for _, e in row])
        out[direction], out[f"evidence_{direction}"] = values, evidence
    return out


def corner(result: dict) -> str:
    """The first header cell of every matrix: "taxon", with where the values come from and when, so a CSV
    passed on alone still says it (a reader that takes the first column as row names drops the cell)."""
    rule = result.get("value_rule") or {}
    media = " / ".join(rule.get("media") or []) or "no medium"
    which = "every medium" if rule.get("rule") == "all" else media
    meta = result["network"].meta
    return (f"taxon [{'INCOMPLETE; ' if result.get('errors') else ''}values from {which}; "
            f"foodnet {meta.get('tool_version', '')} on {str(meta.get('derived_at', ''))[:10]}]")


def _matrix_csv(taxa, header, rows, fmt=_number, first: str = "taxon") -> str:
    return _csv([first, *header], [[n, *map(fmt, r)] for n, r in zip(taxa, rows, strict=True)])


def counts(result: dict) -> dict:
    """How the cells of the pair split by evidence, for the page and the README."""
    pair = pair_rows(result["network"], result)
    tally = {w: 0 for w in EVIDENCE_WORDS}
    for direction in ("consumed", "produced"):
        for row in pair[f"evidence_{direction}"]:
            for e in row:
                tally[e] += 1
    return tally


def conflicts(result: dict) -> list:
    """One line per value that pools experiments that disagree."""
    net = result["network"]
    out = []
    for (taxon, met, ph), cell in sorted(result["cells"].items()):
        if "conflict" in cell["cautions"]:
            name = net.nodes[taxon].name if taxon in net.nodes else taxon
            out.append(f"{name}, {net.nodes[met].name if met in net.nodes else met}, {ph} phase: "
                       + "; ".join(cell["notes"]))
    return out


def presence_lines(result: dict) -> list:
    net = result["network"]
    out = []
    for (taxon, met, ph), by in sorted(result["presence"].items()):
        value = result["cells"].get((taxon, met, ph))
        for direction, entries in sorted(by.items()):
            if value is not None and value["direction"] == direction:
                continue                      # the value medium shows it already, with a number
            out.append(f"{net.nodes[taxon].name if taxon in net.nodes else taxon} {direction} "
                       f"{net.nodes[met].name if met in net.nodes else met} ({ph} phase) in "
                       + ", ".join(sorted({e['medium'] for e in entries})))
    return out


def rates_csv(result: dict) -> str:
    net = result["network"]
    rows = []
    for t in taxa_rows(net, result):
        r = result["rates"].get(t.id)
        if r:
            rows.append([t.name, f"{r['rate']:.6g}", r["unit"], r["n"], r["source"], " ".join(r["studies"]),
                         r["method"]])
        else:
            rows.append([t.name, NA, "1/h", 0, result["without_a_rate"].get(t.id, ""), "", ""])
    first = "taxon [INCOMPLETE]" if result.get("errors") else "taxon"
    return _csv([first, "growth_rate", "unit", "replicates", "source", "studies", "method"], rows)


def initial_csv(result: dict) -> str:
    net = result["network"]
    rows = []
    for m in metabolite_nodes(net, result):
        v = result["initial"].get(m.id)
        rows.append([m.name, m.chebi_id, *(["NA", "NA", "NA", 0] if v is None else
                                           [f"{v['mean']:.6g}", f"{v['min']:.6g}", f"{v['max']:.6g}", v["n"]])])
    return _csv(["metabolite", "chebi_id", "mean_mM", "min_mM", "max_mM", "replicates"], rows)


def cautions_csv(result: dict, phases=None) -> str:
    """One row per value with a caution or a note: what the evidence matrices cannot hold (one replicate,
    boundaries far apart, coarse sampling, amounts that differ, a conflict and the experiments behind it)."""
    names = result.get("names") or {}
    rows = []
    for (taxon, met, ph), cell in sorted(result["cells"].items()):
        if phases and ph not in phases:
            continue
        if cell["cautions"] or cell["notes"]:
            rows.append([names.get(taxon, taxon), names.get(met, met), ph, cell.get("state") or "",
                         " ".join(cell["cautions"]), "; ".join(cell["notes"])])
    return _csv(["taxon", "metabolite", "phase", "state", "cautions", "notes"], rows)


def biomass_csv(result: dict) -> str:
    net = result["network"]
    rows = []
    for t in taxa_rows(net, result):
        b = (result.get("biomass") or {}).get(t.id)
        rows.append([t.name, *(["NA", "NA", "", "NA", 0] if b is None else
                               [f"{b['start']:.6g}", f"{b['change']:.6g}", b["unit"], f"{b['hours']:.6g}", b["n"]])])
    return _csv(["taxon", "biomass_start", "biomass_change", "unit", "phase_hours", "replicates"], rows)


def crm_phase(result: dict) -> str:
    """The one phase a CRM is parameterized from: the window when one is set, the exponential phase with
    "Both" (a consumer-resource model describes growth), else the phase chosen. The metabolites of the second
    time window keep their own window."""
    phases = [p for p in phases_of(result) if p != "window" or not second_window_metabolites(result)
              or result["network"].meta.get("window")] or ["exponential"]
    return phases[0] if len(phases) == 1 else "exponential"


def readme(result: dict, which: str = "matrices") -> str:
    net = result["network"]
    s = result["settings"]
    rule = result["value_rule"]
    phase = net.meta.get("phase")
    window = net.meta.get("window")
    tally = counts(result)
    lines = [f"{'Consumer-resource model parameters' if which == 'crm' else 'Consumption and production matrices'}"
             f" from foodnet {net.meta.get('tool_version', '')}",
             f"Derived {net.meta.get('derived_at', '')} from {net.meta.get('source_db', 'mGrowthDB')}.", ""]
    if result.get("errors"):
        lines += ["INCOMPLETE: records could not be read from mGrowthDB, so these numbers may lack data. Run the "
                  "search again.", *(f"  * {e}" for e in result["errors"]), ""]
    if result.get("warnings"):
        lines += ["Read first (the page showed these above the result):",
                  *(f"  * {w}" for w in result["warnings"]), ""]
    if which == "crm":
        lines += [f"Phase: {crm_phase(result)}" + (" (the search asked for both phases; a consumer-resource model "
                                                    "describes growth, so the exponential phase is used)"
                                                    if phase == "both" else ""), ""]
    elif window:
        lines += [f"Window: {window[0]:g} h to {window[1]:g} h (set in Advanced settings; it replaces the phases).", ""]
    else:
        lines += [f"Phase: {phase}. Exponential growth ends at the first sample where the culture has risen "
                  f"{s['fraction']:.0%} of the way from its start to its maximum, or earlier, where its growth rate "
                  "has fallen below a tenth of its maximum over two consecutive intervals; the stationary phase runs "
                  "from there to the last metabolite sample.", ""]
    second = result.get("second_window") or {}
    if second.get("metabolites"):
        named = ", ".join(n.name for n in result.get("metabolite_nodes", []) if n.id in second["metabolites"])
        lines += [f"Second time window: {named} measured from {second['label']}, not over the phase or main "
                  "window; their columns say so.", ""]
    values = "booleans (1 = it happened, 0 = measured and it did not, NA = no evidence either way)" \
        if s["booleans"] else ("mM, the mean net change over the phase: the mean of the experiments' means when "
                               "several experiments give a value, else the mean of the replicates")
    spread = ("; with three or more replicates a change needs a one-sided 90% confidence interval on the mean "
              "beyond it (no change: inside it), with two both replicates beyond it, and one replicate decides "
              "nothing; experiments must not contradict each other, or the value is inconclusive"
              if s.get("judge_confidence", True) else "")
    own = ("; " + ", ".join(f"{k} {v:g} mM" for k, v in sorted(_own_limits(s).items())) + " have limits of their own"
           if _own_limits(s) else "")
    lines += [f"Values: {values}.",
              f"Detection limit: a mean change below {s['detection_limit']:g} mM counts as no change{spread}{own}.",
              "Rows are taxa, columns metabolites. In consumed.csv and produced.csv every number is a magnitude, "
              "never negative; in the signed matrix a produced compound is positive and a consumed one negative. "
              "The first header cell names the medium the values come from.",
              "",
              "NA is never zero. A cell is NA when the compound was not assayed for that taxon in the value "
              "medium; when it was assayed but is inconclusive (its replicates, or its experiments, do not agree "
              "on what happened); when it was assayed but its cultures gave no phase (no end of "
              "exponential growth, or no stationary phase)"
              + {"na": "; or when its change was seen only in another medium (presence only)",
                 "true": "; a change seen only in another medium is written TRUE (Advanced settings), "
                         "and its direction is the matrix it stands in",
                 "value": "; a change seen only in another medium is written as the amount measured "
                          "there (Advanced settings), not comparable with the value medium's amounts"
                 }[presence_entries(result)]
              + ". A 0 is a measurement in the value medium; when another medium showed a change that way, its "
              "evidence is seen_elsewhere. The evidence_*.csv files say which, cell by cell: "
              + ", ".join(EVIDENCE_WORDS) + ".",
              "Cells: " + ", ".join(f"{tally[w]} {w.replace('_', ' ')}" for w in EVIDENCE_WORDS) + ".", ""]
    if rule["rule"] == "all":
        lines.append("Media: every medium gives values (Ignore media differences was set): " + ", ".join(rule["media"]))
    elif rule["rule"] == "selected":
        lines.append("Media: values come from the media, experiments or studies named in the second box: "
                     + ", ".join(rule["media"]) + ". Every other medium gives presence only.")
    elif rule["rule"] == "majority_in_scope":
        lines.append("Media: values come from the medium that holds data for the most taxa in the studies or "
                     "experiments named in the second box, " + " / ".join(rule["media"])
                     + ". Their other media give presence only.")
    else:
        lines.append("Media: values come from the medium that holds data for the most taxa, "
                     + " / ".join(rule["media"]) + ". Every other medium gives presence only.")
    lines.append("")
    for title, items in (("Values whose experiments disagree (inconclusive: NA, and no arc):", conflicts(result)),
                         ("Duplicate deposits counted once:", result["duplicates"]),
                         ("Seen in another medium and not in the value medium (presence_only where the value medium "
                          "has no value, seen_elsewhere where it measured no change):", presence_lines(result))):
        if items:
            lines += [title, *(f"  * {x}" for x in items), ""]
    if which == "crm":
        lines += ["Growth rates: each taxon's maximum specific growth rate in monoculture (1/h), the median over "
                  "the replicates whose metabolites gave the values, or else over another monoculture in the same "
                  "medium. growth_rates.csv says which, per taxon. A taxon without a rate needs one from elsewhere.",
                  *(["Entries seen only in another medium were set to TRUE in the matrices; a CRM needs amounts, "
                     "so they are missing here (their evidence says presence_only)."]
                    if presence_entries(result) == "true" else []),
                  "Initial concentrations: each metabolite's concentration at the first sample of the value-medium "
                  "cultures, averaged (initial_concentrations.csv), in mM.",
                  "Biomass: each taxon's growth over the same phase, from its growth curve, in that curve's unit "
                  "(biomass.csv), with the phase's length in hours. A simulation's starting abundance for a taxon "
                  "must be in the same unit.",
                  *(["Values come from every medium pooled (Ignore media differences): the initial concentrations "
                     "mix media and describe none of them, so the R package refuses to build a CRM from these "
                     "unless told to."] if result["value_rule"]["rule"] == "all" else []),
                  *(["These are stationary-phase amounts: uptake without growth, which a consumer-resource model "
                     "reads as growth. The R package refuses to build a CRM from them unless told to."]
                    if crm_phase(result) == "stationary" else []),
                  "",
                  "These are measured amounts. miaSim's simulateConsumerResource takes an efficiency matrix E and "
                  "has no uptake rate: a taxon takes up each resource at up to 1 mM per unit of abundance per hour, "
                  "so the unit of abundance decides how fast it eats. The foodnet R package's crm_efficiency() and "
                  "as_miasim() give each taxon a unit of its own (crm_scale), from these amounts, the biomass "
                  "changes and the growth rates, so that a taxon alone grows at its measured rate, gains its "
                  "measured biomass and makes its measured by-products; crm_backcheck() simulates each taxon alone "
                  "and says how close it comes. foodnet measures no Monod constants, and miaSim's uptake of each "
                  "resource follows them: choose them, and check them with crm_backcheck().", ""]
    return "\n".join(lines).rstrip() + "\n"


def pair_package(result: dict) -> bytes:
    """consumed.csv, produced.csv, their evidence matrices and a README, in one zip."""
    net = result["network"]
    pair = pair_rows(net, result)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        for direction in ("consumed", "produced"):
            z.writestr(f"{direction}.csv", _matrix_csv(pair["taxa"], pair["columns"], pair[direction],
                                                       first=corner(result)))
            z.writestr(f"evidence_{direction}.csv",
                       _matrix_csv(pair["taxa"], pair["columns"], pair[f"evidence_{direction}"], str,
                                   first=corner(result)))
        z.writestr("signed.csv", signed_csv(net, result))
        z.writestr("cautions.csv", cautions_csv(result))
        from .figure import matrices_svg
        z.writestr("matrices.svg", matrices_svg(result))
        z.writestr("README.txt", readme(result))
    return buffer.getvalue()


def _numbers(rows) -> list:
    return [[None if isinstance(v, str) else v for v in row] for row in rows]


def crm_payload(result: dict) -> dict:
    """What Send to R posts and the R package reads: the CRM's matrices with their caveats as data."""
    net = result["network"]
    ph = crm_phase(result)
    pair = pair_rows(net, result, [ph])
    taxa = taxa_rows(net, result)
    mets = [m for m, _, _ in columns(net, result, [ph])]
    initial = [None if result["initial"].get(m.id) is None else result["initial"][m.id]["mean"] for m in mets]
    rates = [result["rates"].get(t.id, {}).get("rate") for t in taxa]
    biomass = [(result.get("biomass") or {}).get(t.id) for t in taxa]
    presence = []
    for i, t in enumerate(taxa):
        for j, m in enumerate(mets):
            for direction in ("consumed", "produced"):
                if pair[f"evidence_{direction}"][i][j] == "presence_only":
                    media = sorted({e["medium"] for e in (result["presence"].get((t.id, m.id, ph)) or {})
                                    .get(direction, [])})
                    presence.append({"taxon": pair["taxa"][i], "resource": m.name, "direction": direction,
                                     "media": media})
    return {
        "format": CRM_FORMAT, "tool": net.meta.get("tool", "foodnet"), "tool_version": net.meta.get("tool_version", ""),
        "derived_at": net.meta.get("derived_at", ""), "source_db": net.meta.get("source_db", ""),
        "phase": ph, "window": net.meta.get("window"),
        "values": "booleans" if _booleans(result) else "mM",
        "detection_limit_mM": result["settings"]["detection_limit"],
        "taxa": pair["taxa"], "taxon_ids": [t.id for t in taxa],
        "resources": [m.name for m in mets], "resource_ids": [m.id for m in mets],
        # a CRM takes numbers: a TRUE entry is no amount, so it is sent as missing (the evidence matrices
        # still say presence_only)
        "consumed": _numbers(pair["consumed"]), "produced": _numbers(pair["produced"]),
        "evidence_consumed": pair["evidence_consumed"], "evidence_produced": pair["evidence_produced"],
        "growth_rates": rates, "growth_rate_unit": "1/h",
        "growth_rate_detail": {pair["taxa"][i]: result["rates"][t.id] for i, t in enumerate(taxa)
                               if t.id in result["rates"]},
        "initial_concentrations": initial, "initial_unit": "mM",
        # each taxon's growth over the same phase, in its growth curve's unit: what turns the amounts into
        # miaSim's yields (crm_efficiency in the R package); a simulation's starting abundance is in that unit
        "biomass_change": [None if b is None else b["change"] for b in biomass],
        "biomass_start": [None if b is None else b["start"] for b in biomass],
        "biomass_unit": [None if b is None else b["unit"] for b in biomass],
        "phase_hours": [None if b is None else b["hours"] for b in biomass],
        "caveats": {"presence_only": presence, "conflicts": conflicts(result),
                    "duplicates": list(result["duplicates"]),
                    "without_a_rate": [pair["taxa"][i] for i, t in enumerate(taxa) if t.id not in result["rates"]],
                    "media": list(result["value_rule"]["media"]), "value_rule": result["value_rule"]["rule"],
                    "searched_both_phases": net.meta.get("phase") == "both",
                    # a medium that does not exist, and uptake without growth: the R package refuses to build a
                    # CRM from either unless told to
                    "mixed_media": result["value_rule"]["rule"] == "all",
                    "stationary_phase": ph == "stationary",
                    "incomplete": bool(result.get("errors")), "errors": list(result.get("errors") or []),
                    "warnings": list(result.get("warnings") or []),
                    "inconclusive": [{"taxon": pair["taxa"][i], "resource": m.name, "direction": d}
                                     for i, _ in enumerate(taxa) for j, m in enumerate(mets)
                                     for d in ("consumed", "produced")
                                     if pair[f"evidence_{d}"][i][j] == "inconclusive"]},
        "readme": readme(result, "crm"),
        "studies": sorted(net.studies), "settings": result["settings"],
    }


def crm_package(result: dict) -> bytes:
    """The CRM parameters as files: the pair for the CRM's phase, rates, initial concentrations, README."""
    net = result["network"]
    pair = pair_rows(net, result, [crm_phase(result)])
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        for direction in ("consumed", "produced"):
            z.writestr(f"{direction}.csv", _matrix_csv(pair["taxa"], pair["columns"], pair[direction],
                                                       first=corner(result)))
            z.writestr(f"evidence_{direction}.csv",
                       _matrix_csv(pair["taxa"], pair["columns"], pair[f"evidence_{direction}"], str,
                                   first=corner(result)))
        z.writestr("growth_rates.csv", rates_csv(result))
        z.writestr("initial_concentrations.csv", initial_csv(result))
        z.writestr("biomass.csv", biomass_csv(result))
        z.writestr("cautions.csv", cautions_csv(result, [crm_phase(result)]))
        z.writestr("crm.json", json.dumps(crm_payload(result), indent=1))
        z.writestr("README.txt", readme(result, "crm"))
    return buffer.getvalue()
