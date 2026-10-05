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
  * NA: no value. Either the compound was never assayed for this taxon in the value medium, or its change
    was seen only in another medium (presence only). The evidence matrices of the pair say which.

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
EVIDENCE_WORDS = ("measured", "below_limit", "presence_only", "not_assayed")
CRM_FORMAT = "foodnet.crm/v0"


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


def columns(net: FoodNetwork, result: dict, phases=None) -> list:
    """(metabolite node, phase, label) for each column: metabolites by name, one column per phase."""
    phases = phases or phases_of(result) or ["exponential"]
    mets = metabolite_nodes(net, result)
    labels = _labels(mets)
    if len(phases) == 1:
        return [(m, phases[0], label) for m, label in zip(mets, labels, strict=True)]
    return [(m, ph, f"{label} ({ph})") for m, label in zip(mets, labels, strict=True) for ph in phases]


def _booleans(result: dict) -> bool:
    return bool(result["settings"].get("booleans"))


def entry(result: dict, taxon: str, met: str, ph: str, direction: str) -> tuple:
    """(value, evidence) of one cell of the consumed or the produced matrix."""
    cell = result["cells"].get((taxon, met, ph))
    seen = (result["presence"].get((taxon, met, ph)) or {}).get(direction)
    booleans = _booleans(result)
    if cell is not None and cell["n"]:
        if cell["direction"] == direction:
            return (1 if booleans else abs(cell["mean"])), "measured"
        if seen and booleans:
            return 1, "presence_only"
        return 0, "below_limit"
    if seen:
        return (1 if booleans else None), "presence_only"
    return None, "not_assayed"


def signed_entry(result: dict, taxon: str, met: str, ph: str):
    """The signed matrix's cell: the mean change, 0, or None (NA)."""
    cell = result["cells"].get((taxon, met, ph))
    booleans = _booleans(result)
    if cell is not None and cell["n"]:
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
    return None


def _number(value) -> str:
    if value is None:
        return NA
    if isinstance(value, int):
        return str(value)
    return f"{value:.6g}"


def _csv(header, rows) -> str:
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return out.getvalue()


def signed_rows(net: FoodNetwork, result: dict, phases=None) -> tuple:
    taxa = taxa_rows(net, result)
    cols = columns(net, result, phases)
    rows = [[signed_entry(result, t.id, m.id, ph) for m, ph, _ in cols] for t in taxa]
    return _labels(taxa), [label for _, _, label in cols], rows


def signed_csv(net: FoodNetwork, result: dict) -> str:
    """The taxa by metabolites matrix as CSV: produced positive, consumed negative, NA for no value."""
    names, header, rows = signed_rows(net, result)
    return _csv(["taxon", *header], [[n, *map(_number, r)] for n, r in zip(names, rows, strict=True)])


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


def _matrix_csv(taxa, header, rows, fmt=_number) -> str:
    return _csv(["taxon", *header], [[n, *map(fmt, r)] for n, r in zip(taxa, rows, strict=True)])


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
    return _csv(["taxon", "growth_rate", "unit", "replicates", "source", "studies", "method"], rows)


def initial_csv(result: dict) -> str:
    net = result["network"]
    rows = []
    for m in metabolite_nodes(net, result):
        v = result["initial"].get(m.id)
        rows.append([m.name, m.chebi_id, *(["NA", "NA", "NA", 0] if v is None else
                                           [f"{v['mean']:.6g}", f"{v['min']:.6g}", f"{v['max']:.6g}", v["n"]])])
    return _csv(["metabolite", "chebi_id", "mean_mM", "min_mM", "max_mM", "replicates"], rows)


def crm_phase(result: dict) -> str:
    """The one phase a CRM is parameterized from: the window when one is set, the exponential phase with
    "Both" (a consumer-resource model describes growth), else the phase chosen."""
    phases = phases_of(result) or ["exponential"]
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
    if which == "crm":
        lines += [f"Phase: {crm_phase(result)}" + (" (the search asked for both phases; a consumer-resource model "
                                                    "describes growth, so the exponential phase is used)"
                                                    if phase == "both" else ""), ""]
    elif window:
        lines += [f"Window: {window[0]:g} h to {window[1]:g} h (set in Advanced settings; it replaces the phases).", ""]
    else:
        lines += [f"Phase: {phase}. Exponential growth ends at the first sample where the culture reaches "
                  f"{s['fraction']:.0%} of its maximal abundance; the stationary phase runs from there to the last "
                  "metabolite sample.", ""]
    values = "booleans (1 = it happened, 0 = measured and it did not, NA = no evidence either way)" \
        if s["booleans"] else "mM, the mean net change over the phase across replicates"
    lines += [f"Values: {values}.",
              f"Detection limit: a mean change below {s['detection_limit']:g} mM counts as no change.",
              "Rows are taxa, columns metabolites. In consumed.csv and produced.csv every number is a magnitude, "
              "never negative; in the signed matrix a produced compound is positive and a consumed one negative.",
              "",
              "NA is never zero. A cell is NA when the compound was not assayed for that taxon in the value "
              "medium, or when its change was seen only in another medium (presence only). The evidence_*.csv "
              "files say which, cell by cell: measured, below_limit, presence_only, not_assayed.",
              f"Cells: {tally['measured']} measured, {tally['below_limit']} measured below the limit or the other "
              f"way, {tally['presence_only']} presence only, {tally['not_assayed']} not assayed.", ""]
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
    for title, items in (("Values that pool experiments that disagree (they are pooled anyway; read them first):",
                          conflicts(result)),
                         ("Duplicate deposits counted once:", result["duplicates"]),
                         ("Seen only in another medium (NA in the value matrices):", presence_lines(result))):
        if items:
            lines += [title, *(f"  * {x}" for x in items), ""]
    if which == "crm":
        lines += ["Growth rates: each taxon's maximum specific growth rate in monoculture (1/h), the median over "
                  "the replicates whose metabolites gave the values, or else over another monoculture in the same "
                  "medium. growth_rates.csv says which, per taxon. A taxon without a rate needs one from elsewhere.",
                  "Initial concentrations: each metabolite's concentration at the first sample of the value-medium "
                  "cultures, averaged (initial_concentrations.csv), in mM.",
                  "",
                  "These are measured amounts, not model parameters. A consumer-resource model such as miaSim's "
                  "simulateConsumerResource takes an efficiency matrix E (positive for consumption, negative for "
                  "production). The foodnet R package's crm_efficiency() builds one from these matrices; how to scale "
                  "it is a modeling choice, so nothing here does it for you.", ""]
    return "\n".join(lines).rstrip() + "\n"


def pair_package(result: dict) -> bytes:
    """consumed.csv, produced.csv, their evidence matrices and a README, in one zip."""
    net = result["network"]
    pair = pair_rows(net, result)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
        for direction in ("consumed", "produced"):
            z.writestr(f"{direction}.csv", _matrix_csv(pair["taxa"], pair["columns"], pair[direction]))
            z.writestr(f"evidence_{direction}.csv",
                       _matrix_csv(pair["taxa"], pair["columns"], pair[f"evidence_{direction}"], str))
        z.writestr("signed.csv", signed_csv(net, result))
        from .figure import matrices_svg
        z.writestr("matrices.svg", matrices_svg(result))
        z.writestr("README.txt", readme(result))
    return buffer.getvalue()


def crm_payload(result: dict) -> dict:
    """What Send to R posts and the R package reads: the CRM's matrices with their caveats as data."""
    net = result["network"]
    ph = crm_phase(result)
    pair = pair_rows(net, result, [ph])
    taxa = taxa_rows(net, result)
    mets = [m for m, _, _ in columns(net, result, [ph])]
    initial = [None if result["initial"].get(m.id) is None else result["initial"][m.id]["mean"] for m in mets]
    rates = [result["rates"].get(t.id, {}).get("rate") for t in taxa]
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
        "consumed": pair["consumed"], "produced": pair["produced"],
        "evidence_consumed": pair["evidence_consumed"], "evidence_produced": pair["evidence_produced"],
        "growth_rates": rates, "growth_rate_unit": "1/h",
        "growth_rate_detail": {pair["taxa"][i]: result["rates"][t.id] for i, t in enumerate(taxa)
                               if t.id in result["rates"]},
        "initial_concentrations": initial, "initial_unit": "mM",
        "caveats": {"presence_only": presence, "conflicts": conflicts(result),
                    "duplicates": list(result["duplicates"]),
                    "without_a_rate": [pair["taxa"][i] for i, t in enumerate(taxa) if t.id not in result["rates"]],
                    "media": list(result["value_rule"]["media"]), "value_rule": result["value_rule"]["rule"],
                    "searched_both_phases": net.meta.get("phase") == "both"},
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
            z.writestr(f"{direction}.csv", _matrix_csv(pair["taxa"], pair["columns"], pair[direction]))
            z.writestr(f"evidence_{direction}.csv",
                       _matrix_csv(pair["taxa"], pair["columns"], pair[f"evidence_{direction}"], str))
        z.writestr("growth_rates.csv", rates_csv(result))
        z.writestr("initial_concentrations.csv", initial_csv(result))
        z.writestr("crm.json", json.dumps(crm_payload(result), indent=1))
        z.writestr("README.txt", readme(result, "crm"))
    return buffer.getvalue()
