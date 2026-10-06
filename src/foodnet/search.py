"""One search: taxa (and optionally media) in, a taxon and metabolite network out.

`run_query` is what the page and the command line both call, so the two cannot disagree. It is pure apart
from the client it is given: names to taxon ids (foodnet.taxonomy), taxon ids to studies (mGrowthDB search),
the batch monocultures of those taxa (foodnet.reading), and the derivation (foodnet.derive).
"""
from __future__ import annotations

from collections import defaultdict

from . import crm, rates
from . import derive as d
from . import media as media_rules
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
    # a second window for the metabolites named here (comma separated), from its start to its end in hours,
    # an end of None meaning until the last sample (Karoline, 2026-10-05, for trehalose)
    "second_window_metabolites": "", "second_window_start": 0.0, "second_window_end": None,
    "detection_limit": d.DETECTION_LIMIT,
    # a change needs a one-sided 90% confidence interval on the mean beyond the limit, no change the interval
    # inside it, and experiments must not contradict each other; else inconclusive (Karoline, 2026-10-06).
    # Off: the mean alone decides, as in 0.1.0
    "judge_confidence": True,
    # detection limits of their own, for compounds measured at another scale: "thiamine=0.01, riboflavin=0.005"
    "compound_limits": "",
    "ignore_media": False, "booleans": False,
    # off by default: a rate costs a fit per growth curve; CRM mode turns it on
    "report_rates": False, "rate_method": rates.DEFAULT_METHOD, "rate_window": rates.DEFAULT_WINDOW,
    "merge_arcs": False, "min_studies": 1, "merge_genera": False,
    "conditions": "", "exclude_studies": "", "exclude_experiments": "",
    # tell media apart by the alterations their descriptions state and by a recorded atmosphere (Karoline,
    # 2026-10-06), since mGrowthDB does not report a medium's composition systematically
    "strict_media": True,
    # off: a filled second box limits the search to what matches it; on: everything else the taxa were grown
    # in adds presence-only evidence (Karoline, 2026-10-05)
    "outside_evidence": False,
    # how a matrix cell seen only in another medium is written: "na" (cautious), "true", or "value" (the change
    # measured there, drawn on its own background in the image) (Karoline, 2026-10-06)
    "presence_entries": "na",
    # nothing is left out by default: a metabolite with a known measurement problem is handled in mGrowthDB
    # (Karoline, 2026-10-04: "we don't want to skip any metabolites by default")
    "exclude_metabolites": "",
    "include_non_batch": False, "spike_factor": phases.SPIKE_FACTOR, "correction": "bh",
}
EXAMPLE = ("Escherichia coli LF82", "Bacteroides fragilis", "Roseburia intestinalis")


def second_window_of(s: dict) -> dict | None:
    names = [n.strip() for n in (s.get("second_window_metabolites") or "").split(",") if n.strip()]
    if not names:
        return None
    return {"names": names, "start": float(s.get("second_window_start") or 0.0),
            "end": None if s.get("second_window_end") is None else float(s["second_window_end"])}


def second_window_label(second: dict | None) -> str:
    if not second:
        return ""
    end = "the last sample" if second["end"] is None else f"{second['end']:g} h"
    return f"{second['start']:g} h to {end}"


def compound_limits_of(s: dict) -> dict:
    """{name, lowercased: limit in mM} from the setting "thiamine=0.01, riboflavin: 0.005"; a malformed entry
    raises ValueError with what is wrong, so the page and the command line can say so."""
    out = {}
    for part in (s.get("compound_limits") or "").replace("\n", ",").replace(";", ",").split(","):
        if not part.strip():
            continue
        name, sep, value = part.replace(":", "=").partition("=")
        try:
            limit = float(value)
        except ValueError:
            limit = None
        if not sep or not name.strip() or limit is None or limit < 0:
            raise ValueError(f"a compound's detection limit is written name=mM, as thiamine=0.01; {part.strip()!r} "
                             "is not")
        out[" ".join(name.casefold().split())] = limit
    return out


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
                n=arc["n"], n_experiments=arc.get("n_experiments"), p_value=arc.get("p_value"),
                q_value=arc.get("q_value"),
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
    # experiments the user excludes by id (Advanced settings), like whole studies
    excluded_exps = {e.strip().upper() for e in s["exclude_experiments"].replace("\n", ",").split(",") if e.strip()}
    if excluded_exps:
        gone = sorted({c.experiment for c in cultures + read["growth_only"] if c.experiment.upper() in excluded_exps})
        if gone:
            skipped.append(("excluded experiments", ", ".join(gone) + ": left out, as Exclude these experiments asks"))
        cultures = [c for c in cultures if c.experiment.upper() not in excluded_exps]
        read["growth_only"] = [c for c in read["growth_only"] if c.experiment.upper() not in excluded_exps]
    # what makes two experiments' media one medium: the names, and with the strict rule (default) the
    # alterations their descriptions state and the recorded atmosphere (foodnet.media)
    everyone = cultures + read["growth_only"]
    idents = media_rules.assign_atmospheres([
        media_rules.identity(read["experiments"].get(c.experiment, {}), s["strict_media"]) for c in everyone])
    for c, ident in zip(everyone, idents, strict=True):
        c.medium_key, c.medium = ident["key"], ident["label"]
    grown_in = sorted({c.medium or "unnamed medium" for c in cultures})   # for the "nothing matched" note
    # Two questions, kept apart (Karoline, 2026-10-05, reproducing a hand-checked reference from four ids and
    # Wilkins-Chalgren): ids set the SCOPE, the data looked at; a medium name chooses the VALUE MEDIUM
    # within it. Without a medium name the majority rule runs inside the scope, so the scope's other media
    # give presence, as in the reference. A medium name alone is a scope as well.
    ids = {"studies": selection["studies"], "experiments": selection["experiments"], "media": []}
    media = {"studies": [], "experiments": [], "media": selection["media"]}
    scope = ids if not selecting.empty(ids) else media

    def record(c):
        return read["experiments"].get(c.experiment, {"id": c.experiment, "studyId": c.study})

    if limited:
        def named(c):
            return selecting.matches(record(c), scope)
        left = [c for c in cultures if not named(c)]
        if left:
            skipped.append(("experiments outside the second box", ", ".join(sorted({c.experiment for c in left}))
                            + ": left out, since the second box limits the search to what matches it (Advanced "
                            "settings: Include supporting evidence outside the second box)"))
        cultures = [c for c in cultures if named(c)]
        read["growth_only"] = [c for c in read["growth_only"] if named(c)]
    lost = getattr(index, "failed", [])
    if lost:
        for entry in resolved["unresolved"]:
            resolved["reasons"][entry] = (resolved["reasons"].get(entry, "") + " (the species list could not be read "
                                          "whole, so it may be held in a study that was not read)").strip()
    if lost and not all_studies:
        errors.append(f"{len(lost)} record(s) of the species list could not be read from mGrowthDB (for example "
                      f"{lost[0][0]}: {lost[0][1]}), so taxa held only there are missing from this search; the result "
                      "is incomplete, so run the search again")
    elif lost:
        errors.append(f"{len(lost)} record(s) could not be read while listing mGrowthDB's studies (for example "
                      f"{lost[0][0]}: {lost[0][1]}); the result is incomplete, so run the search again")
    failed = [x for x in skipped if x[1].startswith(UNREAD)]
    if failed:
        errors.append(f"{len(failed)} record(s) could not be read from mGrowthDB (for example {failed[0][0]}: "
                      f"{failed[0][1]}); the result is incomplete, so run the search again")

    say(len(studies), len(studies), "Deriving production and consumption")
    window = window_of(s)
    second = second_window_of(s)
    rows, skips = d.changes(cultures, s["phase"], window, s["fraction"], s["no_growth_factor"], s["spike_factor"],
                            second, s["detection_limit"])
    second_mids = sorted({r["metabolite"] for r in rows if r.get("second")})
    unmatched = [n for n in (second or {}).get("names", [])
                 if not any(d.second_window_matches(met, [n]) for c in cultures for met in c.metabolites.values())]
    skipped += skips
    dropped, duplicate_lines = d.duplicates(cultures)
    rows = [r for r in rows if (cultures[r["culture"]].experiment, r["metabolite"]) not in dropped]
    valued = {r["culture"] for r in rows if r["change"] is not None}
    # the value medium is chosen among the cultures inside the scope; with outside evidence on, the cultures
    # outside it stay in the list and give presence only
    inside = [i for i, c in enumerate(cultures) if selecting.empty(selection) or selecting.matches(record(c), scope)]
    rule = d.value_set([cultures[i] for i in inside], read["experiments"],
                       media if selection["media"] else None, s["ignore_media"],
                       valued=[j for j, i in enumerate(inside) if i in valued])
    rule["chosen"] = [inside[i] for i in rule["chosen"]]
    if not selecting.empty(selection) and not selection["media"] and rule["rule"] == "majority":
        rule["rule"] = "majority_in_scope"
    chosen = set(rule["chosen"])
    value_rows = [r for r in rows if r["culture"] in chosen]
    other_rows = [r for r in rows if r["culture"] not in chosen]
    limit, spread = s["detection_limit"], s["judge_confidence"]
    by_name = compound_limits_of(s)
    limits = {mid: by_name[n] for c in cultures for mid, met in c.metabolites.items()
              for n in by_name if d.second_window_matches(met, [n])}
    unknown_limits = [n for n in by_name
                      if not any(d.second_window_matches(met, [n]) for c in cultures for met in c.metabolites.values())]
    value_cells = d.pool(value_rows, cultures, limit, spread, limits)
    d.adjust(value_cells, s["correction"])
    presence_cells = {} if s["ignore_media"] else d.presence(other_rows, cultures, limit, spread, limits)
    per_study = None
    if not s["merge_arcs"]:
        by_study = defaultdict(list)
        for r in value_rows:
            by_study[cultures[r["culture"]].study].append(r)
        per_study = {sid: d.pool(group, cultures, limit, spread, limits) for sid, group in by_study.items()}
        # the arcs are per study, so their tests are the family corrected together: every study's cells
        # (the pooled cells, corrected above, are the matrices' values)
        family = {(sid, *key): cell for sid, cells in per_study.items() for key, cell in cells.items()}
        d.adjust(family, s["correction"])

    # a strain is shown by its current name, the one its most recent study uses (grownet #24): study
    # SMGDB00000004 still calls taxon 411483 Faecalibacterium prausnitzii A2-165
    for c in cultures + read["growth_only"]:
        if c.taxon.get("identity") == "ncbi" and str(c.taxon.get("taxon_id", "")).isdigit():
            name = current.get(int(c.taxon["taxon_id"]))
            if name and name != c.taxon["name"]:
                c.taxon = {**c.taxon, "name": name, "species": genus_species(name)}
    taxa = {c.taxon["id"]: c.taxon for c in cultures}
    metabolites = {}
    for c in cultures:
        for mid, met in c.metabolites.items():
            metabolites.setdefault(mid, (met["name"], met["chebi_id"]))
    matrix_cells, matrix_presence = value_cells, presence_cells
    node_taxa = {tid: _node_for_taxon(t) for tid, t in taxa.items()}
    if s["merge_genera"]:
        matrix_cells = d.merge_genus_cells(value_cells, taxa, limit, limits)
        matrix_presence = d.merge_genus_presence(presence_cells, taxa)
        if per_study is not None:
            per_study = {sid: d.merge_genus_cells(cells, taxa, limit, limits) for sid, cells in per_study.items()}
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
    # the growth over the phase a CRM is parameterized from (foodnet.matrix.crm_phase)
    crm_ph = "window" if window else ("exponential" if s["phase"] == "both" else s["phase"])
    biomass = crm.biomass_changes([(i, cultures[i]) for i in sorted(chosen)], value_rows, crm_ph)
    if s["merge_genera"]:
        biomass = {}            # strains of one genus grow in units and to densities that do not average
    net.meta["initial_concentrations_mM"] = {k: round(v["mean"], 6) for k, v in initial.items()}

    second_info = None if not second else {**second, "metabolites": second_mids, "label": second_window_label(second),
                                           "unmatched": unmatched}
    net.meta["second_window"] = None if not second_info else {
        "metabolites": [metabolites[m][0] for m in second_mids if m in metabolites], "start": second["start"],
        "end": second["end"]}
    warnings = _warnings(value_cells, presence_cells, rule, window, cultures, chosen, grown_in,
                         not selecting.empty(selection))
    warnings += _value_medium_warnings(rule, cultures, chosen, valued, taxa, rows, window)
    if unknown_limits:
        warnings.append("Detection limits of their own name " + ", ".join(unknown_limits) + ", which no culture of "
                        "these taxa measured; check the spelling (the report lists every metabolite read).")
    if unmatched:
        warnings.append("The second time window names " + ", ".join(unmatched) + ", which no culture of these taxa "
                        "measured; check the spelling (the report lists every metabolite read).")
    if second_info and second_mids:
        warnings.append(f"{', '.join(metabolites[m][0] for m in second_mids if m in metabolites)}: measured over the "
                        f"second time window, {second_info['label']}, not over the phase or main window.")
    # a result built on what could not all be read says so in its files, not only on the screen (the command
    # line also exits with an error, unless told to accept it)
    net.meta["incomplete"] = bool(errors)
    net.meta["errors"] = list(errors)
    net.meta["warnings"] = list(warnings)       # what the page shows above the result, kept with the files
    say(len(studies), len(studies), "Preparing the result")
    return {"entries": names, "settings": dict(s), "resolved": resolved["resolved"],
            "reasons": resolved["reasons"], "suggestions": resolved["suggestions"], "genera": resolved["genera"],
            "unresolved": resolved["unresolved"], "taxon_ids": resolved["taxon_ids"], "all": all_studies,
            "excluded": left_out, "studies": studies, "network": net, "cells": matrix_cells,
            "taxa_nodes": taxa_nodes, "metabolite_nodes": metabolite_nodes, "second_window": second_info,
            # every taxon and metabolite with data by name, whether or not it has an arc in the network
            "names": {n.id: n.name for n in taxa_nodes + metabolite_nodes},
            "presence": matrix_presence, "value_rule": rule, "duplicates": duplicate_lines,
            "rates": organism_rates, "without_a_rate": without_rate, "initial": initial, "biomass": biomass,
            "cultures": len(cultures), "value_cultures": len(chosen), "warnings": warnings,
            "skipped": skipped, "errors": errors}


def _value_medium_warnings(rule, cultures, chosen, valued, taxa, rows, window) -> list:
    """What a reader must know about where the values come from (Karoline, 2026-10-06, after a review showed
    a taxon's values changing with the other taxa searched): the taxa whose own best medium is another one,
    a close choice, and the taxa whose cultures in the value medium give no phase."""
    out = []
    if rule["rule"] in ("majority", "majority_in_scope") and rule["keys"]:
        value_key = rule["keys"][0]
        elsewhere = sorted((taxa[t]["name"], label) for t, label in d.elsewhere(cultures, valued, value_key).items()
                           if t in taxa)
        if elsewhere:
            out.append("Values come from " + " / ".join(rule["media"]) + ", the medium with values for the most "
                       "taxa. " + "; ".join(f"{name} has more data in {label}" for name, label in elsewhere)
                       + ": its values here are from the value medium or missing, and would change if it were "
                       "searched with other taxa. Name a medium in the second box to fix the choice.")
        runner = rule.get("runner_up")
        if runner and not rule["tie"] and runner[1] >= rule["taxa_per_medium"].get(value_key, 0) - 1:
            out.append(f"The value medium won narrowly: {rule['taxa_per_medium'].get(value_key, 0)} taxa against "
                       f"{runner[1]} for {rule.get('labels', {}).get(runner[0], runner[0])}.")
    if window is None:
        lost = defaultdict(set)
        for r in rows:
            if r["culture"] in chosen and "no_phase" in r["cautions"]:
                lost[cultures[r["culture"]].taxon["id"]].add(r["culture"])
        whole = sorted(taxa[t]["name"] for t, idx in lost.items()
                       if t in taxa and not any(i in valued for i in chosen if cultures[i].taxon["id"] == t))
        if whole:
            out.append(", ".join(whole) + ": no end of exponential growth was found on their growth curves (they "
                       "did not grow by the no-growth factor, or the curve is too short), so their metabolites give "
                       "no value in the phases (no_phase in the matrices). A time window in Advanced settings needs "
                       "no phase and gives them values.")
    return out


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
        out.append(f"{beyond} value(s) start before the first or end after the last metabolite sample, so the "
                   "nearest sample stands in for that end of the phase or window (caution window_beyond_data).")
    conflicts = sum(1 for c in value_cells.values() if "conflict" in c["cautions"])
    if conflicts:
        out.append(f"{conflicts} value(s) pool experiments that disagree on the direction or on whether the "
                   "compound changed beyond the detection limit (caution conflict): they are inconclusive, NA in "
                   "the matrices and no arc; the report names the experiments.")
    unsure = sum(1 for c in value_cells.values()
                 if c.get("state") == "inconclusive" and "conflict" not in c["cautions"])
    if unsure:
        out.append(f"{unsure} value(s) are inconclusive: their replicates do not agree on a change beyond the "
                   "detection limit, or on none, so they are neither an arc nor a measured zero (NA, evidence "
                   "inconclusive).")
    coarse = sum(1 for c in value_cells.values() if "coarse_sampling" in c["cautions"])
    if coarse:
        out.append(f"{coarse} value(s) rest on a phase boundary placed on fewer than three growth samples (caution "
                   "coarse_sampling): the phase may end anywhere between two samples.")
    n = sum(1 for k, c in value_cells.items() if k[2] == "exponential" and "still_changing" in c["cautions"])
    if n:
        out.append(f"{n} exponential-phase value(s) miss uptake or production that went on after growth slowed "
                   "(caution still_changing; the growth rate ended exponential growth well before the culture "
                   "stopped): their compounds' use is split between the phases, so a model built on the "
                   "exponential phase alone, as the CRM parameters are, underrates it. A time window over both "
                   "gives the whole change.")
    apart = sum(1 for c in value_cells.values() if "boundaries_differ" in c["cautions"])
    if apart:
        out.append(f"{apart} value(s) average replicates whose phase boundaries lie further apart than a sampling "
                   "interval (caution boundaries_differ); consider a time window.")
    if rule.get("also_matched"):
        out.append("The second box also matched " + "; ".join(rule["also_matched"]) + ", which "
                   + ("differ from the medium giving the values by what their descriptions say was added or taken "
                      "away, or by their atmosphere; they give presence only. Name them in the second box to use "
                      "them for values, or exclude them in Advanced settings."))
    if rule["tie"]:
        labels = rule.get("labels", {})
        out.append("Media tied for the most taxa: " + "; ".join(labels.get(k, k) for k in rule["tie"])
                   + f". Values come from {labels.get(rule['keys'][0], rule['keys'][0])} (more replicates with "
                   "values, then the name that sorts first); name a medium in the second box to choose.")
    return out
