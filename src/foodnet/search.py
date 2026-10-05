"""One search: taxa (and optionally media) in, a taxon and metabolite network out.

`run_query` is what the page and the command line both call, so the two cannot disagree. It is pure apart
from the client it is given: names to taxon ids (foodnet.taxonomy), taxon ids to studies (mGrowthDB search),
the batch monocultures of those taxa (foodnet.reading), and the derivation (foodnet.derive).
"""
from __future__ import annotations

from collections import defaultdict

from . import crm, rates
from . import derive as d
from . import phase as phases
from . import selection as selecting
from .mgrowthdb import MGrowthDBError, data_versions, provenance
from .model import Edge, FoodNetwork, Node, Study, genus_name, genus_species
from .reading import UNREAD, read_cultures
from .taxonomy import resolve_species, species_index, split_entries

PHASE_LABELS = {"exponential": "Exponential phase", "stationary": "Stationary phase", "both": "Both"}
DEFAULTS = {
    "phase": "exponential",
    # None: no time window; a window (start and end, in hours) overrides the phases
    "window_start": None, "window_end": None,
    "fraction": phases.FRACTION, "no_growth_factor": phases.NO_GROWTH_FACTOR,
    "detection_limit": d.DETECTION_LIMIT,
    "ignore_media": False, "booleans": False,
    # off by default: a rate costs a fit per growth curve; CRM mode turns it on
    "report_rates": False, "rate_method": rates.DEFAULT_METHOD, "rate_window": rates.DEFAULT_WINDOW,
    "merge_arcs": False, "min_studies": 1, "merge_genera": False,
    "conditions": "", "exclude_studies": "",
    # off: a filled second box limits the search to what matches it; on: everything else the taxa were grown
    # in adds presence-only evidence (Karoline, 2026-10-05)
    "outside_evidence": False,
    # nothing is left out by default: a metabolite with a known measurement problem is handled in mGrowthDB
    # (Karoline, 2026-10-04: "we don't want to skip any metabolites by default")
    "exclude_metabolites": "",
    "include_non_batch": False, "spike_factor": phases.SPIKE_FACTOR, "correction": "bh",
}
EXAMPLE = ("Escherichia coli LF82", "Bacteroides fragilis", "Roseburia intestinalis")


def window_of(s: dict) -> tuple | None:
    if s.get("window_start") is None or s.get("window_end") is None:
        return None
    return float(s["window_start"]), float(s["window_end"])


def _node_for_taxon(taxon: dict) -> Node:
    return Node(id=taxon["id"], kind="taxon", name=taxon["name"], identity=taxon["identity"],
                taxon_id=taxon.get("taxon_id", ""), species=taxon.get("species", ""))


def _metabolite_node(mid: str, name: str, chebi: str) -> Node:
    return Node(id=mid, kind="metabolite", name=name, identity="chebi" if chebi else "metabolite_name",
                chebi_id=chebi)


def _arc_edge(arc: dict) -> Edge:
    produced = arc["direction"] == d.PRODUCED
    return Edge(source=arc["taxon"] if produced else arc["metabolite"],
                target=arc["metabolite"] if produced else arc["taxon"],
                direction=arc["direction"], phase=arc["phase"], evidence=arc["evidence"],
                amount=d.finite(arc["amount"]), change=d.finite(arc["change"]), sd=d.finite(arc["sd"]),
                n=arc["n"], p_value=arc.get("p_value"), q_value=arc.get("q_value"),
                window_start=arc.get("window_start"), window_end=arc.get("window_end"),
                exponential_h=d.finite(arc.get("exponential_h")), medium=arc["medium"],
                study_ids=tuple(arc["study_ids"]), experiments=tuple(arc["experiments"]),
                cautions=tuple(arc["cautions"]), notes=tuple(arc["notes"]), merged_arcs=arc.get("merged_arcs"),
                merged_taxa=tuple(arc.get("merged_taxa", ())))


def run_query(client, entries, settings: dict | None = None, index=None, progress=None,
              all_studies: bool = False) -> dict:
    """Taxa (names or NCBI taxon ids) to a network, with everything the page, the downloads and the R side
    need. Failures that concern one study are collected in "errors", so one bad study does not lose the rest.
    `all_studies` ignores the entries and reads every batch monoculture with metabolites in mGrowthDB."""
    def say(done, total, message):
        if progress:
            progress(done, total, message)

    s = {**DEFAULTS, **(settings or {})}
    names = [] if all_studies else split_entries(entries)
    say(0, None, "Looking up the taxa in mGrowthDB")
    index = species_index(client) if index is None else index
    resolved = resolve_species(names, index)
    current = getattr(index, "current", {})
    resolved["resolved"] = [(entry, {t: current.get(t, n) for t, n in matches.items()})
                            for entry, matches in resolved["resolved"]]
    errors = []
    selection = selecting.parse(s["conditions"])
    studies = list(selection["studies"])
    if all_studies:
        studies = list(dict.fromkeys(studies + list(getattr(index, "studies", []))))
    elif resolved["taxon_ids"]:
        if hasattr(index, "studies_of"):
            found = index.studies_of(resolved["taxon_ids"])
        else:
            try:
                found = client.search(strain_ncbi_ids=",".join(str(t) for t in resolved["taxon_ids"])).get(
                    "studies", [])
            except MGrowthDBError as e:
                errors.append(f"search failed: {e}")
                found = []
        studies = list(dict.fromkeys(studies + list(found)))
    # What is entered in the second box is a limit: only data matching it are used (Karoline, 2026-10-05,
    # after "when I gave a list of studies, the results also included studies that were not in my list":
    # "by default, when something is entered in the 2nd field, only data matching what was entered are
    # shown ... but in advanced settings, we can switch on showing supporting evidence from other studies").
    # With the box empty, all data are considered. With ids only, only their studies need to be read.
    limited = not selecting.empty(selection) and not s["outside_evidence"]
    if limited and not selection["media"]:
        listed = set(selection["studies"])
        for eid in selection["experiments"]:
            try:
                listed.add(str(client.get_experiment(eid).get("studyId", "")))
            except MGrowthDBError as e:
                errors.append(f"{eid}: {e}")
        studies = [sid for sid in dict.fromkeys(list(selection["studies"]) + studies) if sid in listed]
    excluded = {sid.strip().upper() for sid in s["exclude_studies"].split(",") if sid.strip()}
    left_out = [sid for sid in studies if sid.upper() in excluded]
    studies = [sid for sid in studies if sid.upper() not in excluded]

    wanted = {genus_species(name) for _, matches in resolved["resolved"] for name in matches.values()}
    wanted_ids = {str(t) for t in resolved["taxon_ids"]}

    def keep(name, taxon):
        # a strain the search asked for, by taxon id or by genus and species (a renamed strain keeps its id)
        return all_studies or str(taxon) in wanted_ids or genus_species(name) in wanted

    from .fetch import prefetch_studies
    prefetch_studies(client, studies, keep, s["include_non_batch"], progress=say)
    excluded_metabolites = [m.strip() for m in s["exclude_metabolites"].split(",") if m.strip()]
    read = read_cultures(client, studies, keep, excluded_metabolites, s["include_non_batch"], progress=say)
    cultures, skipped = read["cultures"], list(read["skipped"])
    grown_in = sorted({c.medium or "unnamed medium" for c in cultures})   # for the "nothing matched" note
    if limited:
        def named(c):
            return selecting.matches(read["experiments"].get(c.experiment, {"id": c.experiment, "studyId": c.study}),
                                     selection)
        left = [c for c in cultures if not named(c)]
        if left:
            skipped.append(("experiments outside the second box", ", ".join(sorted({c.experiment for c in left}))
                            + ": left out, since the second box limits the search to what matches it (Advanced "
                            "settings: Include supporting evidence outside the second box)"))
        cultures = [c for c in cultures if named(c)]
        read["growth_only"] = [c for c in read["growth_only"] if named(c)]
    failed = [x for x in skipped if x[1].startswith(UNREAD)]
    if failed:
        errors.append(f"{len(failed)} record(s) could not be read from mGrowthDB (for example {failed[0][0]}: "
                      f"{failed[0][1]}); the result is incomplete, so run the search again")

    say(len(studies), len(studies), "Deriving production and consumption")
    window = window_of(s)
    rows, skips = d.changes(cultures, s["phase"], window, s["fraction"], s["no_growth_factor"], s["spike_factor"])
    skipped += skips
    dropped, duplicate_lines = d.duplicates(cultures)
    rows = [r for r in rows if (cultures[r["culture"]].experiment, r["metabolite"]) not in dropped]
    rule = d.value_set(cultures, read["experiments"], selection if not all_studies or s["conditions"] else None,
                       s["ignore_media"])
    chosen = set(rule["chosen"])
    value_rows = [r for r in rows if r["culture"] in chosen]
    other_rows = [r for r in rows if r["culture"] not in chosen]
    limit = s["detection_limit"]
    value_cells = d.pool(value_rows, cultures, limit)
    d.adjust(value_cells, s["correction"])
    presence_cells = {} if s["ignore_media"] else d.presence(other_rows, cultures, limit)
    per_study = None
    if not s["merge_arcs"]:
        by_study = defaultdict(list)
        for r in value_rows:
            by_study[cultures[r["culture"]].study].append(r)
        per_study = {sid: d.pool(group, cultures, limit) for sid, group in by_study.items()}
        # one correction runs over one family of tests, the pooled cells; a per-study arc takes its cell's
        # q-value when that study is the only one behind the cell (the two values are then the same), and
        # has none when the cell pools several studies
        for sid, cells in per_study.items():
            for key, cell in cells.items():
                pooled = value_cells.get(key) or {}
                cell["q_value"] = pooled.get("q_value") if pooled.get("studies") == [sid] else None

    taxa = {c.taxon["id"]: c.taxon for c in cultures}
    metabolites = {}
    for c in cultures:
        for mid, met in c.metabolites.items():
            metabolites.setdefault(mid, (met["name"], met["chebi_id"]))
    matrix_cells, matrix_presence = value_cells, presence_cells
    node_taxa = {tid: _node_for_taxon(t) for tid, t in taxa.items()}
    if s["merge_genera"]:
        matrix_cells = d.merge_genus_cells(value_cells, taxa)
        matrix_presence = d.merge_genus_presence(presence_cells, taxa)
        if per_study is not None:
            per_study = {sid: d.merge_genus_cells(cells, taxa) for sid, cells in per_study.items()}
        node_taxa = {}
        for tid, t in taxa.items():
            gid = d.genus_of(taxa, tid)
            node_taxa.setdefault(gid, Node(id=gid, kind="taxon", name=genus_name(t["name"]), identity="genus"))
    arc_list = d.arcs(matrix_cells, matrix_presence, per_study, s["booleans"])
    if s["merge_genera"]:
        for arc in arc_list:
            arc["merged_taxa"] = (matrix_cells.get((arc["taxon"], arc["metabolite"], arc["phase"])) or {}).get(
                "merged_taxa", [])
    arc_list, below_min = d.min_studies_filter(arc_list, s["min_studies"])

    used_studies = sorted({sid for c in cultures for sid in [c.study]})
    meta = {**provenance(), "source_db": "mGrowthDB (live)", "query": "all" if all_studies else "taxa",
            "taxa": names, "studies": studies, "settings": dict(s), "selection": selection,
            "phase": "window" if window else s["phase"], "window": list(window) if window else None,
            "detection_limit_mM": limit, "value_rule": {k: v for k, v in rule.items() if k != "chosen"},
            "duplicates": duplicate_lines, "values": "booleans" if s["booleans"] else "mM",
            "hidden": {"below_min_studies": below_min}}
    net = FoodNetwork(meta=meta)
    # the network holds the nodes its arcs touch; the matrices and the CRM hold every taxon and metabolite
    # measured, since a measured no change is information there and clutter in a picture
    taxa_nodes = list(node_taxa.values())
    metabolite_nodes = [_metabolite_node(mid, name, chebi)
                        for mid, (name, chebi) in sorted(metabolites.items(), key=lambda kv: kv[1][0])]
    touched = {end for arc in arc_list for end in (arc["taxon"], arc["metabolite"])}
    for node in taxa_nodes + metabolite_nodes:
        if node.id in touched:
            net.add_node(node)
    for sid in used_studies:
        try:
            st = client.get_study(sid)
        except MGrowthDBError:
            st = {}
        net.add_study(Study(id=sid, citation=st.get("name", sid), url=st.get("url", "")))
    for arc in arc_list:
        net.add_edge(_arc_edge(arc))
    net.meta["data"] = data_versions(client, studies, net.meta["derived_at"])

    organism_rates, without_rate = {}, {}
    if s["report_rates"]:
        say(len(studies), len(studies), "Reading the growth rates")
        value_cultures = [cultures[i] for i in sorted(chosen)]
        others = [c for i, c in enumerate(cultures) if i not in chosen] + read["growth_only"]
        organism_rates, without_rate, rate_skips = crm.growth_rates(
            list(taxa), value_cultures, others, set(rule["keys"]), s["ignore_media"], s["rate_method"],
            s["rate_window"], s["spike_factor"])
        skipped += rate_skips
        if s["merge_genera"]:
            organism_rates, without_rate = _genus_rates(organism_rates, without_rate, taxa)
        net.meta["growth_rates"] = {"method": rates.method_name(s["rate_method"], s["rate_window"]),
                                    "rates": organism_rates, "without_a_rate": without_rate}
    initial = crm.initial_concentrations([cultures[i] for i in sorted(chosen)])
    net.meta["initial_concentrations_mM"] = {k: round(v["mean"], 6) for k, v in initial.items()}

    warnings = _warnings(value_cells, presence_cells, rule, window, cultures, chosen, grown_in,
                         not selecting.empty(selection))
    say(len(studies), len(studies), "Preparing the result")
    return {"entries": names, "settings": dict(s), "resolved": resolved["resolved"],
            "reasons": resolved["reasons"], "suggestions": resolved["suggestions"], "genera": resolved["genera"],
            "unresolved": resolved["unresolved"], "taxon_ids": resolved["taxon_ids"], "all": all_studies,
            "excluded": left_out, "studies": studies, "network": net, "cells": matrix_cells,
            "taxa_nodes": taxa_nodes, "metabolite_nodes": metabolite_nodes,
            "presence": matrix_presence, "value_rule": rule, "duplicates": duplicate_lines,
            "rates": organism_rates, "without_a_rate": without_rate, "initial": initial,
            "cultures": len(cultures), "value_cultures": len(chosen), "warnings": warnings,
            "skipped": skipped, "errors": errors}


def _genus_rates(found, missing, taxa):
    """Rates per genus: the median of its strains' rates."""
    import statistics
    by_genus = defaultdict(list)
    for tid, r in found.items():
        by_genus[d.genus_of(taxa, tid)].append(r)
    out = {g: {**rs[0], "rate": statistics.median(r["rate"] for r in rs), "n": sum(r["n"] for r in rs)}
           for g, rs in by_genus.items()}
    gone = {d.genus_of(taxa, tid): why for tid, why in missing.items() if d.genus_of(taxa, tid) not in out}
    return out, gone


def _warnings(value_cells, presence_cells, rule, window, cultures, chosen, grown_in=(), boxed=False) -> list:
    """What a reader must see above the result, not only in the report."""
    out = []
    if boxed and not chosen and grown_in:
        # the second box matched none of these taxa's monocultures: say so, and what it could have matched,
        # rather than showing a result that looks like the previous one (Karoline, 2026-10-05: "the search
        # is not updated when I relaunch the same species but with another input in the medium field")
        out.append("Nothing in the second box matches a monoculture of these taxa, so nothing gives a value. "
                   "Their monocultures were grown in: " + " · ".join(grown_in) + ".")
    short = sorted({cultures[i].experiment_name or cultures[i].experiment for i in chosen
                    if any(phases.short_record(m["series"]) for m in cultures[i].metabolites.values())})
    if short:
        out.append(f"Metabolites were recorded for less than {phases.SHORT_RECORD_H:g} h in {len(short)} "
                   f"experiment(s) ({', '.join(short[:5])}{' and more' if len(short) > 5 else ''}); their values "
                   "are used, and the arcs carry the caution short_record.")
    beyond = sum(1 for c in value_cells.values() if "window_beyond_data" in c["cautions"])
    if beyond:
        out.append(f"{beyond} value(s) end after the last metabolite sample, so the last sample stands in for "
                   "the end of the phase or window (caution window_beyond_data).")
    conflicts = sum(1 for c in value_cells.values() if "conflict" in c["cautions"])
    if conflicts:
        out.append(f"{conflicts} value(s) pool experiments that disagree on the direction or on whether the "
                   "compound changed beyond the detection limit (caution conflict); the report names them.")
    if rule["tie"]:
        out.append("Media tied for the most taxa: " + ", ".join(rule["tie"]) + f". Values come from "
                   f"{rule['keys'][0]}; name a medium in the second box to choose.")
    return out
