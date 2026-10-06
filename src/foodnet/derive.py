"""From monoculture metabolite series to produced and consumed arcs.

The steps, each a decision recorded in docs/METHOD_NOTES.md:

1. **A change per replicate and phase** (`changes`): the metabolite's net change over the growth phase
   asked for (`foodnet.phase`), in mM.
2. **Which data give values** (`value_set`, Karoline, 2026-10-04): the medium that supplies data for most
   of the taxa; or, when the second box is filled, the media (or experiments, or studies) named there; or
   every medium, with "Ignore media differences". Data from any other medium say only whether the compound
   was produced or consumed, never how much ("the majority medium was used for values and the other media
   were only used for presence/absence").
3. **Duplicate deposits** (`duplicates`): the same experiment deposited under two studies is counted once.
   R. intestinalis in Wilkins-Chalgren is study SMGDB00000007 `ri2` and study SMGDB00000002 `RI_WC`, the
   second rounded to two decimals; pooling both counts three replicates twice.
4. **Pooling** (`pool`): replicates of one taxon, metabolite and phase are pooled across the experiments and
   studies of the value set ("can be pooled with the same medium, but when there's a conflict, this needs to
   be reported"). A conflict is experiments that disagree on what happened: one produced the compound, one
   consumed it, or one changed it beyond the detection limit and another did not.
5. **The detection limit** (0.2 mM by default, an advanced setting): a mean change smaller than this, in
   either direction, is a measurement of nothing, not an arc.
6. **Arcs** (`arcs`): a measured arc for a pooled change beyond the limit; a presence_only arc for a change
   seen beyond the limit only in another medium. One arc per study by default, as in grownet; "Merge arcs
   across studies" gives one per taxon, metabolite, phase and direction.

Never assayed is never zero: a compound not on a taxon's panel has no value anywhere, and the matrices write
it as NA.
"""
from __future__ import annotations

import math
import statistics
from collections import defaultdict

from . import phase as phases
from . import selection as selecting
from .model import genus_name
from .stats import CORRECTIONS, paired

DETECTION_LIMIT = 0.2          # mM; where the data of the hand-checked reference matrices break
PHASE_CHOICES = ("exponential", "stationary", "both")
PRODUCED, CONSUMED = "produced", "consumed"
MEASURED, PRESENCE_ONLY = "measured", "presence_only"
# how two series count as the same deposit: every value within this many mM plus this share of itself
DUPLICATE_ABS, DUPLICATE_REL, DUPLICATE_MIN_POINTS = 0.011, 0.01, 3


# ---- 1. a change per replicate and phase ----------------------------------------------------------

def boundaries(cultures, fraction: float = phases.FRACTION, factor: float = phases.NO_GROWTH_FACTOR,
               spike_factor: float = phases.SPIKE_FACTOR) -> tuple:
    """({culture index: boundary or None}, skipped): where exponential growth ends in each culture.

    A culture without a usable growth curve takes the median boundary of the other replicates of its
    experiment, and says so (`phase_from_other_replicates`); with none to take, it has no boundary and is
    reported."""
    found, skipped, own = {}, [], {}
    for i, c in enumerate(cultures):
        if c.growth is None:
            continue
        ratio = phases.spike(c.growth["values"], spike_factor)
        if ratio:
            skipped.append((c.label, f"growth curve spikes ({ratio:.0f} times its neighbors); its phase boundary "
                            "is taken from the other replicates"))
            continue
        try:
            own[i] = phases.exponential_end(c.growth["times"], c.growth["values"], fraction, factor)
        except phases.NoBoundary as e:
            skipped.append((c.label, f"no end of exponential growth: {e}"))
    by_experiment = defaultdict(list)
    for i, b in own.items():
        by_experiment[cultures[i].experiment].append(b)
    for i, c in enumerate(cultures):
        if i in own:
            found[i] = own[i]
            continue
        siblings = by_experiment.get(c.experiment)
        if siblings:
            found[i] = {"end": statistics.median(b["end"] for b in siblings),
                        "last": all(b["last"] for b in siblings), "borrowed": True}
            if c.growth is None:
                skipped.append((c.label, "no growth curve in this replicate; the median phase boundary of its "
                                "experiment's other replicates is used"))
        else:
            found[i] = None
            if c.growth is None:
                skipped.append((c.label, "no growth curve in this replicate or its experiment, so its growth "
                                "phases are unknown; its metabolites are left out (a time window in Advanced "
                                "settings does not need one)"))
    return found, skipped


def second_window_matches(met: dict, names) -> bool:
    """Whether a metabolite is one the second window names: by its name, or the name its study recorded it
    under (so "acetic acid" and "acetate" both name the joined acetate), ignoring case and spacing."""
    wanted = {" ".join(n.casefold().split()) for n in names}
    return any(" ".join((x or "").casefold().split()) in wanted
               for x in (met.get("name"), met.get("recorded_name")))


def changes(cultures, phase: str = "exponential", window: tuple | None = None,
            fraction: float = phases.FRACTION, factor: float = phases.NO_GROWTH_FACTOR,
            spike_factor: float = phases.SPIKE_FACTOR, second: dict | None = None) -> tuple:
    """(rows, skipped): one row per culture, metabolite and phase.

    A row: {"culture" (index), "taxon", "metabolite", "metabolite_name", "chebi_id", "phase", "change",
    "start", "end", "initial", "exponential_h", "cautions"}. `window` (start, end) in hours replaces the
    phases. `exponential_h` is how long the culture grew exponentially, from its first growth sample to the
    end of exponential growth (Karoline, 2026-10-04: "report the duration of the exponential phase"); it is
    reported with a time window too, where the boundary only describes the culture.

    `second` is the second window, {"names", "start", "end"} with `end` None for "until the last sample":
    the metabolites it names take it instead of the phases or the main window (Karoline, 2026-10-05, for
    trehalose, which the hand-checked reference measures over the whole run). Their rows carry phase "window" and
    "second": True."""
    rows, skipped = [], []
    bounds, skips = boundaries(cultures, fraction, factor, spike_factor)
    if window is None:
        skipped += skips
    for i, c in enumerate(cultures):
        boundary = bounds.get(i)
        duration = exponential_hours(c, boundary)
        if window is not None:
            boundary = None
        for mid, met in c.metabolites.items():
            series = met["series"]
            short = phases.short_record(series)
            in_second = bool(second and second.get("names") and second_window_matches(met, second["names"]))
            if in_second:
                end = second["end"] if second.get("end") is not None else series[-1][0]
                windows = phases.phase_windows(series, None, phase, (second.get("start") or 0.0, end))
            else:
                windows = phases.phase_windows(series, boundary, phase, window)
            for name, (start, end, cautions) in windows.items():
                cautions = list(cautions)
                if short:
                    cautions.append("short_record")
                if boundary and boundary.get("borrowed"):
                    cautions.append("phase_from_other_replicates")
                if start is None:
                    rows.append({"culture": i, "taxon": c.taxon["id"], "metabolite": mid,
                                 "metabolite_name": met["name"], "chebi_id": met["chebi_id"], "phase": name,
                                 "change": None, "start": None, "end": None, "initial": series[0][1],
                                 "exponential_h": duration, "cautions": cautions})
                    continue
                d = phases.change(series, start, end)
                if d["beyond"] and "window_beyond_data" not in cautions:
                    cautions.append("window_beyond_data")
                rows.append({"culture": i, "taxon": c.taxon["id"], "metabolite": mid,
                             "metabolite_name": met["name"], "chebi_id": met["chebi_id"], "phase": name,
                             "change": d["change"], "start": start, "end": end, "initial": d["initial"],
                             "exponential_h": duration, "second": in_second, "cautions": cautions})
    return rows, skipped


def exponential_hours(culture, boundary) -> float | None:
    """Hours from the culture's first growth sample (its first metabolite sample without a curve) to the end
    of exponential growth, or None when there is no boundary."""
    if not boundary:
        return None
    if culture.growth:
        first = culture.growth["times"][0]
    else:
        first = min((m["series"][0][0] for m in culture.metabolites.values()), default=0.0)
    return boundary["end"] - first


# ---- 2. which data give values -------------------------------------------------------------------

def value_set(cultures, experiments: dict, selection: dict | None = None, ignore_media: bool = False) -> dict:
    """The rule that decides which cultures give values, and its outcome.

    Returns {"rule": "all" | "selected" | "majority", "media": [names], "keys": [medium keys],
    "taxa_per_medium": {key: n}, "tie": [keys tied for most taxa], "chosen": [culture indices]}.

    The majority medium is the one whose monocultures cover the most taxa (strains) with metabolite data.
    A tie goes to the medium with more replicates, then to the name that sorts first, and is reported.
    """
    taxa = defaultdict(set)
    count = defaultdict(int)
    names = defaultdict(set)
    for c in cultures:
        taxa[c.medium_key].add(c.taxon["id"])
        count[c.medium_key] += 1
        names[c.medium_key].add(c.medium or "unnamed medium")
    per_medium = {k: len(v) for k, v in taxa.items()}
    if ignore_media:
        chosen = list(range(len(cultures)))
        return {"rule": "all", "media": sorted({n for v in names.values() for n in v}), "keys": sorted(taxa),
                "taxa_per_medium": per_medium, "tie": [], "chosen": chosen}
    if not selecting.empty(selection):
        # each entry gives the values from the medium it matches for the most taxa: "Wilkins-Chalgren" also
        # matches Wilkins-Chalgren with mucin added, which the strict medium rule tells apart (Karoline,
        # 2026-10-06), and a reader typing a medium means one medium. The others it matches are named, and
        # give presence only. Ids are matched exactly and choose their experiments as they are.
        chosen, also = set(), set()
        for kind in ("media", "experiments", "studies"):
            for entry in selection.get(kind, ()):
                one = {"media": [], "experiments": [], "studies": [], kind: [entry]}
                matched = [i for i, c in enumerate(cultures)
                           if selecting.matches(experiments.get(c.experiment, {"id": c.experiment, "studyId": c.study}),
                                                one)]
                if kind != "media":
                    chosen.update(matched)
                    continue
                by_key = defaultdict(list)
                for i in matched:
                    by_key[cultures[i].medium_key].append(i)
                if not by_key:
                    continue
                best = sorted(by_key, key=lambda k: (-len({cultures[i].taxon["id"] for i in by_key[k]}),
                                                     -len(by_key[k]), k))[0]
                chosen.update(by_key[best])
                also.update(cultures[i].medium for k, idx in by_key.items() if k != best for i in idx)
        chosen = sorted(chosen)
        keys = sorted({cultures[i].medium_key for i in chosen})
        media = sorted({cultures[i].medium for i in chosen})
        return {"rule": "selected", "media": media, "keys": keys, "taxa_per_medium": per_medium, "tie": [],
                "chosen": chosen, "also_matched": sorted(also - set(media))}
    if not taxa:
        return {"rule": "majority", "media": [], "keys": [], "taxa_per_medium": {}, "tie": [], "chosen": []}
    ranked = sorted(taxa, key=lambda k: (-per_medium[k], -count[k], k))
    best = ranked[0]
    tie = [k for k in ranked if per_medium[k] == per_medium[best]]
    chosen = [i for i, c in enumerate(cultures) if c.medium_key == best]
    return {"rule": "majority", "media": sorted(names[best]), "keys": [best], "taxa_per_medium": per_medium,
            "tie": tie if len(tie) > 1 else [], "chosen": chosen}


# ---- 3. duplicate deposits ------------------------------------------------------------------------

def _same_series(a, b) -> bool:
    times = {t for t, _ in a} & {t for t, _ in b}
    if len(times) < DUPLICATE_MIN_POINTS:
        return False
    da, db = dict(a), dict(b)
    return all(abs(da[t] - db[t]) <= DUPLICATE_ABS + DUPLICATE_REL * abs(da[t]) for t in times)


def _contained(inner, outer, shared) -> bool:
    """Whether every replicate series of `inner` matches one of `outer`, for every shared compound."""
    return all(all(any(_same_series(ci.metabolites[m]["series"], co.metabolites[m]["series"])
                       for co in outer if m in co.metabolites)
                   for ci in inner if m in ci.metabolites)
               for m in shared)


def duplicates(cultures) -> tuple:
    """(dropped, found): the (experiment, metabolite) pairs to leave out, and one line per duplicate found.

    Two experiments of one taxon in one medium are the same deposit when they share at least two compounds
    and, for every shared compound, every replicate series of one matches a replicate series of the other at
    their common time points (within 0.011 mM plus 1%, which covers rounding to two decimals). One direction
    is enough: study SMGDB00000002's RI_WC holds one acetate series twice, so it is contained in study
    SMGDB00000007's ri2 while ri2 is not contained in it. The deposit that holds the other is kept; when each
    holds the other, the one with more replicates (the earlier study on a tie). The other gives only the
    compounds the kept one did not measure."""
    by_exp = defaultdict(list)
    for c in cultures:
        by_exp[c.experiment].append(c)
    keys = sorted(by_exp, key=lambda e: (by_exp[e][0].study, e))
    dropped, found = set(), []
    for x, ea in enumerate(keys):
        for eb in keys[x + 1:]:
            a, b = by_exp[ea], by_exp[eb]
            if a[0].taxon["id"] != b[0].taxon["id"] or a[0].medium_key != b[0].medium_key:
                continue
            shared = set().union(*(c.metabolites for c in a)) & set().union(*(c.metabolites for c in b))
            if len(shared) < 2:
                continue
            b_in_a, a_in_b = _contained(b, a, shared), _contained(a, b, shared)
            if not (b_in_a or a_in_b):
                continue
            if b_in_a != a_in_b:                 # one holds the other: keep the one that holds
                keep, drop = (a, b) if b_in_a else (b, a)
            else:
                keep, drop = (a, b) if len(a) >= len(b) else (b, a)
            for m in shared:
                dropped.add((drop[0].experiment, m))
            found.append(f"{drop[0].experiment_name or drop[0].experiment} ({drop[0].study}) repeats "
                         f"{keep[0].experiment_name or keep[0].experiment} ({keep[0].study}) on "
                         f"{len(shared)} compound(s): counted once, from {keep[0].study}")
    return dropped, found


# ---- 4. pooling ------------------------------------------------------------------------------------

def _class(mean: float | None, limit: float) -> int:
    if mean is None or abs(mean) < limit:
        return 0
    return 1 if mean > 0 else -1


WORD = {1: "produced", -1: "consumed", 0: "below the limit"}


def pool(rows, cultures, limit: float = DETECTION_LIMIT) -> dict:
    """{(taxon, metabolite, phase): cell} from the rows of one group of cultures.

    A cell: {"mean", "sd", "n", "values", "direction" (produced, consumed or None), "experiments",
    "studies", "media", "cautions", "notes", "start", "end", "p_value", "initial"}."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r["taxon"], r["metabolite"], r["phase"])].append(r)
    cells = {}
    for key, members in groups.items():
        valued = [r for r in members if r["change"] is not None]
        cautions = sorted({c for r in members for c in r["cautions"]})
        experiments = sorted({cultures[r["culture"]].experiment for r in members})
        studies = sorted({cultures[r["culture"]].study for r in members})
        media = sorted({cultures[r["culture"]].medium for r in members})
        exponential = _mean_of(r.get("exponential_h") for r in members)
        if not valued:
            cells[key] = {"mean": None, "sd": None, "n": 0, "values": [], "direction": None,
                          "experiments": experiments, "studies": studies, "media": media, "cautions": cautions,
                          "notes": [], "start": None, "end": None, "p_value": None, "initial": None,
                          "exponential_h": exponential}
            continue
        values = [r["change"] for r in valued]
        mean = statistics.mean(values)
        sd = statistics.stdev(values) if len(values) > 1 else None
        test = paired(values, [0.0] * len(values))
        notes = []
        per_exp = defaultdict(list)
        for r in valued:
            per_exp[cultures[r["culture"]].experiment].append(r["change"])
        classes = {e: _class(statistics.mean(v), limit) for e, v in per_exp.items()}
        if len(set(classes.values())) > 1:
            cautions = sorted(set(cautions) | {"conflict"})
            names = {c.experiment: (c.experiment_name or c.experiment) for c in cultures}
            notes.append("experiments disagree: " + "; ".join(
                f"{names.get(e, e)} {WORD[k]} ({statistics.mean(per_exp[e]):+.2f} mM)"
                for e, k in sorted(classes.items())))
        if len(values) == 1:
            cautions = sorted(set(cautions) | {"single_replicate"})
        k = _class(mean, limit)
        cells[key] = {"mean": mean, "sd": sd, "n": len(values), "values": values,
                      "direction": PRODUCED if k > 0 else CONSUMED if k < 0 else None,
                      "experiments": experiments, "studies": studies, "media": media, "cautions": cautions,
                      "notes": notes, "start": statistics.mean(r["start"] for r in valued),
                      "end": statistics.mean(r["end"] for r in valued),
                      "p_value": test["p"] if test else None,
                      "initial": statistics.mean(r["initial"] for r in valued), "exponential_h": exponential}
    return cells


def presence(rows, cultures, limit: float = DETECTION_LIMIT) -> dict:
    """{(taxon, metabolite, phase): {direction: [{"medium", "mean", "n", "studies", "experiments"}]}} from
    the rows outside the value set: per medium, the direction of the mean change beyond the limit."""
    by_medium = defaultdict(list)
    for r in rows:
        if r["change"] is None:
            continue
        c = cultures[r["culture"]]
        by_medium[(r["taxon"], r["metabolite"], r["phase"], c.medium_key)].append(r)
    out = defaultdict(lambda: defaultdict(list))
    for (taxon, met, ph, _), members in by_medium.items():
        mean = statistics.mean(r["change"] for r in members)
        k = _class(mean, limit)
        if not k:
            continue
        out[(taxon, met, ph)][PRODUCED if k > 0 else CONSUMED].append({
            "medium": cultures[members[0]["culture"]].medium, "mean": mean, "n": len(members),
            "studies": sorted({cultures[r["culture"]].study for r in members}),
            "experiments": sorted({cultures[r["culture"]].experiment for r in members}),
            "exponential_h": _mean_of(r.get("exponential_h") for r in members),
            "cautions": sorted({x for r in members for x in r["cautions"]})})
    return {k: dict(v) for k, v in out.items()}


def adjust(cells: dict, correction: str = "bh") -> None:
    """Set each cell's q_value: its p_value corrected over every cell of the search (reported, not used)."""
    keys = [k for k, c in cells.items() if c.get("p_value") is not None]
    adjusted = CORRECTIONS[correction][1]([cells[k]["p_value"] for k in keys])
    for k, q in zip(keys, adjusted, strict=True):
        cells[k]["q_value"] = q


# ---- 5 and 6. arcs ---------------------------------------------------------------------------------

def _measured_arc(taxon, met, ph, cell, studies=None, booleans=False) -> dict:
    direction = cell["direction"]
    return {"taxon": taxon, "metabolite": met, "phase": ph, "direction": direction, "evidence": MEASURED,
            "amount": None if booleans else abs(cell["mean"]),
            "change": None if booleans else cell["mean"],
            "sd": None if booleans else cell["sd"], "n": cell["n"],
            "p_value": None if booleans else cell.get("p_value"),
            "q_value": None if booleans else cell.get("q_value"),
            "window_start": cell["start"], "window_end": cell["end"], "exponential_h": cell.get("exponential_h"),
            "medium": "; ".join(cell["media"]),
            "study_ids": studies or cell["studies"], "experiments": cell["experiments"],
            "cautions": list(cell["cautions"]), "notes": list(cell["notes"])}


def _presence_arc(taxon, met, ph, direction, entries, value_cell) -> dict:
    cautions = sorted({x for e in entries for x in e["cautions"]} & {"short_record", "window_beyond_data",
                                                                     "stationary_not_reached"})
    notes = [f"seen in {e['medium']} ({e['mean']:+.2f} mM over {e['n']} replicate(s)); another medium than the "
             "values come from, so only its direction counts" for e in entries]
    if value_cell is not None and value_cell["n"]:
        if value_cell["direction"] is None:
            cautions.append("not_detected_in_value_medium")
        else:
            notes.append(f"in the value medium it was {value_cell['direction']} instead")
    return {"taxon": taxon, "metabolite": met, "phase": ph, "direction": direction, "evidence": PRESENCE_ONLY,
            "amount": None, "change": None, "sd": None, "n": sum(e["n"] for e in entries), "p_value": None,
            "q_value": None, "window_start": None, "window_end": None,
            "exponential_h": _mean_of(e.get("exponential_h") for e in entries),
            "medium": "; ".join(sorted({e["medium"] for e in entries})),
            "study_ids": sorted({s for e in entries for s in e["studies"]}),
            "experiments": sorted({x for e in entries for x in e["experiments"]}),
            "cautions": sorted(set(cautions)), "notes": notes}


def arcs(value_cells: dict, presence_cells: dict, per_study_cells: dict | None = None,
         booleans: bool = False) -> list:
    """The arcs of a search.

    `value_cells` are pooled over the whole value set; `per_study_cells`, when given (merging off), are the
    same pooled per study, {study: cells}, and give one measured arc per study. Presence arcs follow the same
    split: one per study unless merged. A presence arc is added for a direction the value medium did not show
    (absent there, below the limit, or the other way)."""
    out = []
    if per_study_cells is None:
        for (taxon, met, ph), cell in value_cells.items():
            if cell["direction"]:
                out.append(_measured_arc(taxon, met, ph, cell, booleans=booleans))
                out[-1]["merged_arcs"] = len(cell["studies"])
    else:
        for study, cells in sorted(per_study_cells.items()):
            for (taxon, met, ph), cell in cells.items():
                if cell["direction"]:
                    out.append(_measured_arc(taxon, met, ph, cell, [study], booleans))
    for (taxon, met, ph), by_direction in presence_cells.items():
        value_cell = value_cells.get((taxon, met, ph))
        for direction, entries in by_direction.items():
            if value_cell is not None and value_cell["direction"] == direction:
                continue
            if per_study_cells is None:
                out.append(_presence_arc(taxon, met, ph, direction, entries, value_cell))
                out[-1]["merged_arcs"] = len(out[-1]["study_ids"])
            else:
                by_study = defaultdict(list)
                for e in entries:
                    for s in e["studies"]:
                        by_study[s].append(e)
                for study, group in sorted(by_study.items()):
                    arc = _presence_arc(taxon, met, ph, direction, group, value_cell)
                    arc["study_ids"] = [study]
                    out.append(arc)
    return out


def min_studies_filter(arc_list, minimum: int = 1) -> tuple:
    """The arcs resting on at least `minimum` studies, and how many were left out."""
    kept = [a for a in arc_list if len(a["study_ids"]) >= minimum]
    return kept, len(arc_list) - len(kept)


# ---- merging to genus -------------------------------------------------------------------------------

def genus_of(taxa: dict, taxon: str) -> str:
    return f"genus:{genus_name(taxa[taxon]['name'])}"


def merge_genus_cells(cells: dict, taxa: dict) -> dict:
    """Cells keyed by genus: the median of the member taxa's mean changes (assayed taxa only), with the
    direction from that median, and every caution, study and experiment carried along."""
    groups = defaultdict(list)
    for (taxon, met, ph), cell in cells.items():
        groups[(genus_of(taxa, taxon), met, ph)].append((taxon, cell))
    out = {}
    for key, members in groups.items():
        means = [c["mean"] for _, c in members if c["mean"] is not None]
        mean = statistics.median(means) if means else None
        cautions = sorted({x for _, c in members for x in c["cautions"]})
        direction = None
        if mean is not None:
            k = 1 if mean > 0 else -1
            direction = PRODUCED if k > 0 else CONSUMED
            if not any(c["direction"] == direction for _, c in members):
                direction = None
        merged = {"mean": mean, "sd": None, "n": sum(c["n"] for _, c in members), "values": [],
                  "direction": direction, "experiments": sorted({x for _, c in members for x in c["experiments"]}),
                  "studies": sorted({x for _, c in members for x in c["studies"]}),
                  "media": sorted({x for _, c in members for x in c["media"]}), "cautions": cautions,
                  "notes": [f"median over {len(means)} taxon(s): " + ", ".join(sorted(taxa[t]["name"]
                                                                                    for t, _ in members))],
                  "start": _mean_of(c["start"] for _, c in members), "end": _mean_of(c["end"] for _, c in members),
                  "p_value": None, "initial": _mean_of(c["initial"] for _, c in members),
                  "exponential_h": _mean_of(c.get("exponential_h") for _, c in members),
                  "merged_taxa": sorted(t for t, _ in members)}
        out[key] = merged
    return out


def _mean_of(values):
    values = [v for v in values if v is not None]
    return statistics.mean(values) if values else None


def merge_genus_presence(presence_cells: dict, taxa: dict) -> dict:
    out = defaultdict(lambda: defaultdict(list))
    for (taxon, met, ph), by_direction in presence_cells.items():
        for direction, entries in by_direction.items():
            out[(genus_of(taxa, taxon), met, ph)][direction] += entries
    return {k: dict(v) for k, v in out.items()}


def finite(x):
    return None if x is None or (isinstance(x, float) and not math.isfinite(x)) else x
