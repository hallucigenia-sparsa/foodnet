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
import re
import statistics
from collections import defaultdict

from . import compounds
from . import phase as phases
from . import selection as selecting
from .model import genus_name
from .stats import CORRECTIONS, paired, t_quantile

DETECTION_LIMIT = 0.2          # mM; where the data of the hand-checked reference matrices break
PHASE_CHOICES = ("exponential", "stationary", "both")
PRODUCED, CONSUMED = "produced", "consumed"
MEASURED, PRESENCE_ONLY = "measured", "presence_only"
# how two series count as the same deposit: every value within this many mM plus this share of itself
# rounding to two decimals only: a 1% tolerance merged genuine repeats of one protocol (a fourth review round)
DUPLICATE_ABS, DUPLICATE_REL, DUPLICATE_MIN_POINTS = 0.0051, 0.0, 3
DUPLICATE_TIME = 0.02          # hours: sample times that differ by rounding only
# a series identifies a deposit only when it moves by more than this many mM and this share of its level
DUPLICATE_MIN_RANGE, DUPLICATE_REL_RANGE = 1.0, 0.05


# ---- 1. a change per replicate and phase ----------------------------------------------------------

def _own_boundary(c, fraction, factor, spike_factor) -> tuple:
    """(boundary or None, reason): the first of the culture's growth curves, in the order of preference, that
    gives an end of exponential growth. The curve that gave it becomes the culture's growth curve, so its
    growth rate is read from the same curve. A curve too short, spiked, or without growth gives way to the
    next one (a two-point flow cytometry trace must not hide a full OD curve)."""
    reasons = []
    c.not_grown = False
    flat = []
    for curve in c.curves or ([c.growth] if c.growth else []):
        ratio = phases.spike(curve["values"], spike_factor)
        if ratio:
            reasons.append(f"its {curve['technique'] or 'growth'} curve spikes ({ratio:.0f} times its neighbors)")
            continue
        try:
            od = (curve.get("technique") or "").casefold() == "od"
            found = phases.exponential_end(curve["times"], curve["values"], fraction, factor, smoothed=od)
        except phases.NoBoundary as e:
            reasons.append(f"{curve['technique'] or 'growth'} curve: {e}")
            flat.append(str(e).startswith("did not grow") and not phases.grown_by_od(curve))
            continue
        c.growth = curve
        return found, ""
    # every curve said "did not grow", and no OD curve rose by OD_RISE: the culture did not grow
    c.not_grown = bool(flat) and all(flat) and len(flat) == len(reasons)
    return None, "; ".join(reasons)


def boundaries(cultures, fraction: float = phases.FRACTION, factor: float = phases.NO_GROWTH_FACTOR,
               spike_factor: float = phases.SPIKE_FACTOR) -> tuple:
    """({culture index: boundary or None}, skipped): where exponential growth ends in each culture.

    A culture without a usable growth curve takes the median boundary of the other replicates of its
    experiment, and says so (`phase_from_other_replicates`); with none to take, it has no boundary and is
    reported. Replicates of one experiment whose boundaries lie further apart than the experiment's typical
    sampling interval are marked (`differ`), since their values then average different stretches of time."""
    found, skipped, own = {}, [], {}
    for i, c in enumerate(cultures):
        if c.growth is None and not c.curves:
            continue
        boundary, why = _own_boundary(c, fraction, factor, spike_factor)
        if boundary is None:
            skipped.append((c.label, f"no end of exponential growth: {why}; the phase boundary is taken from the "
                            "other replicates of its experiment, if any has one"))
        else:
            own[i] = boundary
    by_experiment = defaultdict(list)
    for i, b in own.items():
        by_experiment[cultures[i].experiment].append((i, b))
    differ = set()
    for exp, members in by_experiment.items():
        ends = [b["end"] for _, b in members]
        steps = [t1 - t0 for i, _ in members for t0, t1 in zip(cultures[i].growth["times"],
                                                              cultures[i].growth["times"][1:], strict=False)]
        if len(ends) > 1 and steps and max(ends) - min(ends) > statistics.median(steps):
            differ.add(exp)
    for i, c in enumerate(cultures):
        if i in own:
            found[i] = {**own[i], "differ": c.experiment in differ}
            continue
        siblings = [b for _, b in by_experiment.get(c.experiment, [])]
        if siblings:
            found[i] = {"end": statistics.median(b["end"] for b in siblings),
                        "last": all(b["last"] for b in siblings), "borrowed": True,
                        "coarse": any(b.get("coarse") for b in siblings), "differ": c.experiment in differ,
                        "by": "rate" if any(b.get("by") == "rate" for b in siblings) else "90%",
                        "moved": any(b.get("moved") for b in siblings),
                        "late": statistics.median(b["late"] for b in siblings if b.get("late") is not None)
                        if any(b.get("late") is not None for b in siblings) else None}
            if c.growth is None:
                skipped.append((c.label, "no growth curve in this replicate; the median phase boundary of its "
                                "experiment's other replicates is used"))
        else:
            found[i] = None
            if c.growth is None:
                skipped.append((c.label, "no growth curve in this replicate or its experiment, so its growth "
                                "phases are unknown; its metabolites give no value in the phases (no_phase; a "
                                "time window in Advanced settings does not need one)"))
    return found, skipped


def second_window_matches(met: dict, names) -> bool:
    """Whether a metabolite is one the second window names: by its name, or the name its study recorded it
    under (so "acetic acid" and "acetate" both name the joined acetate), ignoring case and spacing."""
    def plain(text):
        # "thiamine(1+)", as mGrowthDB records a cation, is thiamine to whoever types it
        return re.sub(r"\(\d*[+\-\u2212]\)$", "", " ".join((text or "").casefold().split())).strip()
    wanted = {plain(n) for n in names}
    return any(plain(x) in wanted for x in (met.get("name"), met.get("recorded_name")))


def _next_change(series, boundary: float, until: float | None = None):
    """The metabolite's change from the end of exponential growth to where the 90% rule alone would have ended
    it (`until`), or to its next sample when that is not later, or None. Where the growth-rate rule ended
    growth, a culture can still be using or making a compound while it slows (E. coli LF82 ferments its glucose
    between 8 and 12 h, after pyruvate runs out, while its cells grow by 15%; Karoline, 2026-10-06: flag it
    rather than move the boundary). The whole stretch counts, so a slow change spread over several samples is
    caught too (Karoline, 2026-10-07: "Check the whole stretch"). `pool` judges it per cell (`still_changing`)."""
    after = [t for t, _ in series if t > boundary]
    if not after:
        return None
    first, _ = phases.value_at(series, boundary)
    target = until if until is not None and until > boundary else after[0]
    target = min(target, series[-1][0])
    if target <= boundary:
        return None
    last, _ = phases.value_at(series, target)
    return last - first


def changes(cultures, phase: str = "exponential", window: tuple | None = None,
            fraction: float = phases.FRACTION, factor: float = phases.NO_GROWTH_FACTOR,
            spike_factor: float = phases.SPIKE_FACTOR, second: dict | None = None,
            limit: float = DETECTION_LIMIT, evaporation: float = phases.EVAPORATION) -> tuple:
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
            elif window is None and boundary is None and c.not_grown and not active(c, limit, evaporation):
                # it did not grow by either rule, and its compounds did not move as metabolism moves them:
                # no value, so drift or a dead inoculum is never an arc
                windows = {ph: (None, None, ["not_grown"]) for ph in
                           (("exponential", "stationary") if phase == "both" else (phase,))}
            elif window is None and boundary is None:
                # no end of exponential growth (a culture read as not grown, as unblanked OD can make it, or no
                # growth curve at all): the change over the whole run, in the column of the phase asked for and
                # marked whole_run, rather than nothing (Karoline, 2026-10-06). With Both it fills the
                # exponential column and the stationary one says no_phase
                column = "stationary" if phase == "stationary" else "exponential"
                why = ("growth_unclear" if c.not_grown else "growth_unknown" if not (c.curves or c.growth)
                       else None)
                windows = {column: (series[0][0], series[-1][0], ["whole_run"] + ([why] if why else []))}
                if phase == "both":
                    windows["stationary"] = (None, None, ["no_phase"])
            else:
                windows = phases.phase_windows(series, boundary, phase, window)
            for name, (start, end, cautions) in windows.items():
                cautions = list(cautions)
                after = None
                if short:
                    cautions.append("short_record")
                if boundary and not in_second and window is None:
                    if boundary.get("borrowed"):
                        cautions.append("phase_from_other_replicates")
                    if boundary.get("coarse"):
                        cautions.append("coarse_sampling")
                    if boundary.get("differ"):
                        cautions.append("boundaries_differ")
                    if boundary.get("moved"):
                        cautions.append("growth_rate_boundary")
                    # the change right after growth slowed, on both phases' rows: the stationary value then
                    # holds it, and the exponential value lacks it (and the CRM takes the exponential value)
                    cut = end if name == "exponential" else start
                    after = (_next_change(series, cut, boundary.get("late"))
                             if (cut is not None and boundary.get("moved")) else None)
                if start is None:
                    rows.append({"culture": i, "taxon": c.taxon["id"], "metabolite": mid,
                                 "metabolite_name": met["name"], "chebi_id": met["chebi_id"], "phase": name,
                                 "change": None, "start": None, "end": None, "initial": series[0][1],
                                 "exponential_h": duration, "cautions": cautions})
                    continue
                d = phases.change(series, start, end)
                # a change beyond the limit but within what the series' own scatter, or (in a culture that did not
                # grow) evaporation, could account for: no measured change, so the value becomes inconclusive. Only
                # changes beyond the limit are judged, so a flat series stays a measured 0, and evaporation only
                # explains rises (a thirteenth review round)
                # the scatter within the window the change spans: across phases a compound made and then used
                # is a trend, not scatter
                inside = [(t, v) for t, v in series if start <= t <= end]
                if abs(d["change"]) >= limit and abs(d["change"]) < 2 * _scatter(inside):
                    cautions.append("within_scatter")
                elif ("growth_unclear" in cautions and d["change"] >= limit
                      and d["change"] <= beyond_evaporation(series, limit, evaporation)):
                    cautions.append("within_evaporation")
                if d["beyond"] and "window_beyond_data" not in cautions:
                    cautions.append("window_beyond_data")
                rows.append({"culture": i, "taxon": c.taxon["id"], "metabolite": mid,
                             "metabolite_name": met["name"], "chebi_id": met["chebi_id"], "phase": name,
                             "change": d["change"], "start": start, "end": end, "initial": d["initial"],
                             "exponential_h": duration, "second": in_second, "cautions": cautions,
                             "after": after if not in_second and window is None else None})
    return rows, skipped


def _ranks(values) -> list:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2
        i = j + 1
    return ranks


TREND = 0.8                    # Spearman's rank correlation with time a trend must reach


def _steady(series, limit: float, sign: int) -> bool:
    """Whether a series moves one way over time, as a trend and not as scatter: Spearman's rank correlation
    with time at least TREND in that direction (a single dip does not undo a rise: butyrate 2.60, 2.23, 2.72,
    3.29, 3.50, 3.75), and more than twice the limit overall."""
    values = [v for _, v in series]
    if len(values) < 4 or sign * (values[-1] - values[0]) <= 2 * limit:
        return False
    a, b = _ranks(list(range(len(values)))), _ranks(values)
    ma, mb = statistics.mean(a), statistics.mean(b)
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return den > 0 and sign * sum((x - ma) * (y - mb) for x, y in zip(a, b, strict=True)) / den >= TREND


def beyond_evaporation(series, limit: float, evaporation: float = phases.EVAPORATION) -> float:
    """The change a culture that did not grow must exceed to count: the detection limit, or the share of the
    compound's level that evaporation could account for over the run (phase.EVAPORATION), whichever is larger."""
    level = max(abs(series[0][1]), abs(series[-1][1]))
    return max(limit, evaporation * level)


def active(culture, limit: float = DETECTION_LIMIT, evaporation: float = phases.EVAPORATION) -> bool:
    """Whether a culture's compounds moved as metabolism moves them over its run: one used up (to below the limit
    from above twice it, or falling steadily) and another made (rising steadily), each by more than evaporation
    could account for (`beyond_evaporation`). A culture whose growth curve shows no growth can still do this (A.
    soehngenii in study SMGDB00000010, on an optical density read without its blank, turns glucose and lactate
    into butyrate; Karoline, 2026-10-06: "Activity counts"), and drift, evaporation or a dead inoculum cannot.
    Endpoints alone are not enough: a twelfth review round found scatter of +/-0.5 mM passing an endpoint test
    in 13 of 20 cultures, and evaporation concentrating one compound while a volatile one fell."""
    series = [m["series"] for m in culture.metabolites.values() if len(m["series"]) > 2]
    # a volatile compound's fall is no uptake: it can leave as vapor
    lasting = [m["series"] for m in culture.metabolites.values()
               if len(m["series"]) > 2 and not compounds.volatile(m.get("name"), m.get("chebi_id"))]

    def moved(s, sign):
        return sign * (s[-1][1] - s[0][1]) > beyond_evaporation(s, limit, evaporation)

    used = any(moved(s, -1) and ((s[0][1] >= 2 * limit and s[-1][1] < limit) or _steady(s, limit, -1))
               for s in lasting)
    made = any(moved(s, 1) and _steady(s, limit, 1) for s in series)
    return used and made


def _scatter(series) -> float:
    """How much a series jumps between samples beyond its trend: the spread (median absolute deviation, scaled to
    a standard deviation) of its sample-to-sample steps, over the square root of two. A steady trend or a single
    sharp step (glucose 10, 10, 0, 0) gives little; scatter around a level (acetate 24.8, 30.9, 25.2, 31.8) gives
    much. 0 for a series of fewer than four samples."""
    v = [x for _, x in series]
    if len(v) < 4:
        return 0.0
    steps = [b - a for a, b in zip(v, v[1:], strict=False)]
    middle = statistics.median(steps)
    return 1.4826 * statistics.median(abs(d - middle) for d in steps) / math.sqrt(2)


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

def value_set(cultures, experiments: dict, selection: dict | None = None, ignore_media: bool = False,
              valued=None) -> dict:
    """The rule that decides which cultures give values, and its outcome.

    Returns {"rule": "all" | "selected" | "majority", "media": [names], "keys": [medium keys],
    "taxa_per_medium": {key: n}, "tie": [keys tied for most taxa], "chosen": [culture indices],
    "runner_up": [key, n] or None}.

    The majority medium is the one whose monocultures give values (a change in the phase or window asked
    for) for the most taxa (strains); `valued` holds the indices of the cultures that give one, and a
    culture without (no phase boundary, say) does not vote (Karoline, 2026-10-06). A tie goes to the medium
    with more replicates giving values, then to the name that sorts first, and is reported.
    """
    taxa = defaultdict(set)
    count = defaultdict(int)
    names = defaultdict(set)
    voters = set(range(len(cultures))) if valued is None else set(valued)
    if not voters:
        voters = set(range(len(cultures)))
    for i, c in enumerate(cultures):
        names[c.medium_key].add(c.medium or "unnamed medium")
        if i in voters:
            taxa[c.medium_key].add(c.taxon["id"])
            count[c.medium_key] += 1
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
        return {"rule": "majority", "media": [], "keys": [], "taxa_per_medium": {}, "tie": [], "chosen": [],
                "runner_up": None}
    ranked = sorted(taxa, key=lambda k: (-per_medium[k], -count[k], k))
    best = ranked[0]
    tie = [k for k in ranked if per_medium[k] == per_medium[best]]
    chosen = [i for i, c in enumerate(cultures) if c.medium_key == best]
    return {"rule": "majority", "media": sorted(names[best]), "keys": [best], "taxa_per_medium": per_medium,
            "tie": tie if len(tie) > 1 else [], "chosen": chosen,
            "runner_up": [ranked[1], per_medium[ranked[1]]] if len(ranked) > 1 else None,
            "labels": {k: " / ".join(sorted(v)) for k, v in names.items()}}


def elsewhere(cultures, valued, value_key: str) -> dict:
    """{taxon: label of another medium}: the taxa with more cultures giving values in another medium than in
    the value medium. They are named in a warning, since their values depend on which other taxa were
    searched with them."""
    count = defaultdict(lambda: defaultdict(int))
    label = {}
    for i in valued:
        c = cultures[i]
        count[c.taxon["id"]][c.medium_key] += 1
        label[c.medium_key] = c.medium or "unnamed medium"
    out = {}
    for taxon, by in count.items():
        best = sorted(by, key=lambda m: (-by[m], m))[0]
        if best != value_key and by[best] > by.get(value_key, 0):
            out[taxon] = label[best]
    return out


# ---- 3. duplicate deposits ------------------------------------------------------------------------

def _common_times(a, b) -> list:
    """[(time in a, time in b)] for the samples of two series taken at the same time, within DUPLICATE_TIME
    hours (8.333 h in one deposit is 8.33 h in another, rounded)."""
    pairs, used = [], set()
    for ta, _ in a:
        near = min(((abs(ta - tb), tb) for tb, _ in b if tb not in used), default=None)
        if near is not None and near[0] <= DUPLICATE_TIME:
            pairs.append((ta, near[1]))
            used.add(near[1])
    return pairs


def _same_series(a, b) -> bool:
    times = _common_times(a, b)
    if len(times) < DUPLICATE_MIN_POINTS:
        return False
    da, db = dict(a), dict(b)
    return all(abs(da[ta] - db[tb]) <= DUPLICATE_ABS + DUPLICATE_REL * abs(da[ta]) for ta, tb in times)


def _in_transit(series) -> int:
    """Whether a series moves enough that matching it is evidence of one deposit: a compound that stays at
    its medium level matches any other culture in the same medium (a review showed two experiments merged on
    two compounds that only wobbled within 1% around 30 and 12 mM), and so does one that only goes from its
    start to exhaustion or to a plateau (glucose 10, 0, 0). So the series must move by more than 1 mM and 5%
    of its level, and its samples in transit, strictly between its lowest and highest values, are what
    count: this returns how many it has (0 for a series that does not move)."""
    values = [v for _, v in series]
    low, high = min(values), max(values)
    margin = max(DUPLICATE_MIN_RANGE, DUPLICATE_REL_RANGE * max(abs(v) for v in values))
    if high - low <= margin:
        return 0
    tolerance = DUPLICATE_ABS + DUPLICATE_REL * high
    return sum(low + tolerance < v < high - tolerance for v in values)


def _contained(inner, outer, shared) -> bool:
    """Whether every replicate series of `inner` matches one of `outer`, for every shared compound."""
    return all(all(any(_same_series(ci.metabolites[m]["series"], co.metabolites[m]["series"])
                       for co in outer if m in co.metabolites)
                   for ci in inner if m in ci.metabolites)
               for m in shared)


def duplicates(cultures) -> tuple:
    """(dropped, found): the (experiment, metabolite) pairs to leave out, and one line per duplicate found.

    Two experiments of one taxon are the same deposit when they share at least two compounds and, for every
    shared compound, every replicate series of one matches a replicate series of the other at their common
    time points (within 0.0051 mM, which covers rounding to two decimals and nothing more, and times within
    DUPLICATE_TIME hours), and the shared compounds move (`_in_transit`): series that sit at the
    medium's level match in any two cultures of one medium, so at least two of the shared samples must be in
    transit (`_in_transit`). Only experiments in one medium are compared: a
    second round of the review showed two genuine experiments in Wilkins-Chalgren with and without mucin
    matching on an exhausted sugar and a plateau, and the value medium's copy dropped. One direction is enough:
    study SMGDB00000002's RI_WC holds one acetate series twice, so it is contained in study SMGDB00000007's
    ri2 while ri2 is not contained in it. The deposit that holds the other is kept; when each holds the
    other, the one with more replicates (the earlier study on a tie). The other gives only the compounds the
    kept one did not measure."""
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
            # at least two samples in transit over the shared compounds: one can match by chance between two
            # cultures of one protocol (a third round of the review merged two genuine repeats on one)
            first = {m: next(c for c in a if m in c.metabolites) for m in shared}
            if sum(_in_transit(first[m].metabolites[m]["series"]) for m in shared) < 2:
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


INCONCLUSIVE = None


CONFIDENCE = 0.9               # one-sided: the interval's lower (or upper) bound must clear the limit


def classify(values, limit: float, agree: bool = True):
    """1 (produced), -1 (consumed), 0 (no change) or None (inconclusive) for the replicates of one experiment,
    changes in mM.

    With `agree` (the default), how many replicates there are decides how (Karoline, 2026-10-06):

      * three or more: a one-sided 90% t-interval on the mean ("Confidence interval"). A change needs its
        near bound beyond the limit, no change needs the whole interval inside it, anything else is
        inconclusive. More replicates narrow the interval, so they make a call easier, never harder, and one
        failed sample widens it without erasing a change four replicates agree on;
      * two (13 of the 32 experiments with metabolites in mGrowthDB, 2026-10-06): both beyond the limit on
        the same side for a change, both inside it for none ("Pairs agree"). The interval on one degree of
        freedom is so wide that clear pairs (-3.2 and -1.5 mM) came out inconclusive;
      * one: inconclusive ("1 is inconclusive"), so the least data never makes the boldest call. Identical
        replicates near the limit count as one: they are one series deposited twice, not certainty (a
        substrate exhausted from one shared start, identical far beyond the limit, stays a change).

    At three or more, replicates that all lie beyond the limit on one side (or all inside it) decide even
    when the interval does not, so a confirming third replicate never undoes what a pair would decide.

    Without `agree` the mean alone decides, as in 0.1.0. Two rules came first and were dropped on review:
    the mean plus and minus one standard deviation (a spread, not an inference), and every replicate beyond
    the limit at any n (one failed sample erased a change)."""
    values = [v for v in values if v is not None]
    if not values:
        return INCONCLUSIVE
    mean = statistics.mean(values)
    if not agree:
        return _class(mean, limit)
    if len(values) > 1 and 0 < abs(mean) < 2 * limit and statistics.stdev(values) == 0:
        values = values[:1]        # near the limit, identical replicates would pass as certainty
    if len(values) == 1:
        return INCONCLUSIVE
    classes = {_class(v, limit) for v in values}
    agreed = classes.pop() if len(classes) == 1 else INCONCLUSIVE
    if len(values) == 2:
        return agreed
    half = t_quantile(CONFIDENCE, len(values) - 1) * statistics.stdev(values) / math.sqrt(len(values))
    low, high = mean - half, mean + half
    if low >= limit:
        return 1
    if high <= -limit:
        return -1
    if -limit < low and high < limit:
        return 0
    # replicates that all agree, as a pair must, are never inconclusive for there being three of them
    return agreed


WORD = {1: "produced", -1: "consumed", 0: "below the limit", None: "inconclusive"}
STATE = {1: "produced", -1: "consumed", 0: "no_change", None: "inconclusive"}


def _units(valued, cultures) -> tuple:
    """({experiment: [replicate changes]}, the values a cell's mean and test are over): the experiments' means
    when the cell pools several experiments, each counted once whatever its replicates (Karoline, 2026-10-06:
    the experiment, not the replicate, is the unit, so a study with ten replicates does not outvote one with
    two), or the replicates of its one experiment."""
    per_exp = defaultdict(list)
    for r in valued:
        per_exp[cultures[r["culture"]].experiment].append(r["change"])
    if len(per_exp) > 1:
        return per_exp, [statistics.mean(v) for v in per_exp.values()]
    return per_exp, [v for vs in per_exp.values() for v in vs]


AMOUNTS_DIFFER = 2.0           # experiments agreeing in direction whose means differ more than this factor
STILL_SHARE = 0.25             # still_changing: the change right after growth slowed is this share of the phase's
START_SHARE = 0.25             # start_differs: replicates' starts differ by more than this share of the phase


def pool(rows, cultures, limit: float = DETECTION_LIMIT, agree: bool = True, limits: dict | None = None) -> dict:
    """{(taxon, metabolite, phase): cell} from the rows of one group of cultures.

    A cell: {"mean", "sd", "n", "n_experiments", "values", "direction" (produced, consumed or None),
    "state" (produced, consumed, no_change or inconclusive), "experiments", "studies", "media", "cautions",
    "notes", "start", "end", "p_value", "low", "high", "initial"}: `low` and `high` are the lowest and highest
    replicate change behind the value (None on one replicate); `start` and `end` the interval it covers, the
    mean over its experiments. `limits` maps a metabolite to its own detection limit.

    Each experiment is judged on its replicates (`classify`). With one experiment, that is the cell. With
    several, every experiment that is not inconclusive must say the same: then that is the cell, its mean is
    the mean of those experiments' means, and experiments whose amounts differ more than twofold add the
    caution amounts_differ (a difference in size is not a difference in what happened). Experiments that say
    different things make the cell inconclusive, with the caution conflict and a note naming them; so does an
    inconclusive experiment whose mean lies beyond the limit on the other side. An inconclusive experiment
    that does not contradict the others is left out and named (experiment_left_out); when every experiment
    is inconclusive, so is the cell. `sd`, the test, `n` and `n_experiments` are over the experiments the
    value rests on."""
    limits = limits or {}
    groups = defaultdict(list)
    for r in rows:
        groups[(r["taxon"], r["metabolite"], r["phase"])].append(r)
    names = {c.experiment: (c.experiment_name or c.experiment) for c in cultures}
    cells = {}
    for key, members in groups.items():
        lim = limits.get(key[1], limit)
        # a replicate within its series' scatter (or evaporation) keeps its vote: dropping it would select
        # replicates by their values and hand the call to the ones left (a fourteenth review round). Only when
        # every value is so does the cell become inconclusive
        noisy = [r for r in members if r["change"] is not None
                 and {"within_scatter", "within_evaporation"} & set(r["cautions"])]
        valued = [r for r in members if r["change"] is not None]
        if noisy and len(noisy) < len(valued):
            noisy = []
        whole = [r for r in valued if "whole_run" in r["cautions"]]
        whole_note = None
        if whole and len(whole) < len(valued):
            # a phase value where any culture has one; the whole-run changes of the others do not mix in
            valued = [r for r in valued if "whole_run" not in r["cautions"]]
            members = [r for r in members if "whole_run" not in r["cautions"]]
            whole_note = ("left out, as their change spans the whole run: " + ", ".join(sorted(
                {cultures[r["culture"]].experiment_name or cultures[r["culture"]].experiment for r in whole})))
        cautions = set(c for r in members for c in r["cautions"])
        experiments = sorted({cultures[r["culture"]].experiment for r in members})
        studies = sorted({cultures[r["culture"]].study for r in members})
        media = sorted({cultures[r["culture"]].medium for r in members})
        exponential = _mean_of(r.get("exponential_h") for r in members)
        if not valued:
            cells[key] = {"mean": None, "sd": None, "n": 0, "n_experiments": 0, "values": [], "direction": None,
                          "state": None, "experiments": experiments, "studies": studies, "media": media,
                          "cautions": sorted(cautions), "notes": [], "start": None, "end": None, "p_value": None,
                          "initial": None, "exponential_h": exponential, "limit": lim}
            continue
        values = [r["change"] for r in valued]
        per_exp, units = _units(valued, cultures)
        classes = {e: classify(v, lim, agree) for e, v in per_exp.items()}
        if agree:
            # a replicate within its series' scatter (or evaporation) keeps its vote, so it can stop a call, but
            # cannot make one: a change needs two replicates outside it that show it themselves (a fifteenth
            # review round: a Bacteroides pair at +0.20 and a flagged +0.21 mM read as butyrate production)
            clean = defaultdict(list)
            for r in valued:
                if not {"within_scatter", "within_evaporation"} & set(r["cautions"]):
                    clean[cultures[r["culture"]].experiment].append(r["change"])
            for e, c in classes.items():
                if c in (1, -1) and sum(_class(v, lim) == c for v in clean.get(e, [])) < 2:
                    classes[e] = INCONCLUSIVE
        said = {c for c in classes.values() if c is not INCONCLUSIVE}
        # an experiment too noisy to decide still disagrees when its mean lies beyond the limit on the other
        # side of what the others say (a review: one replicate made it inconclusive and its -3 mM vanished)
        if len(said) == 1 and next(iter(said)) in (1, -1):
            side = next(iter(said))
            for e, c in classes.items():
                if c is INCONCLUSIVE and _class(statistics.mean(per_exp[e]), lim) == -side:
                    classes[e] = -side
            said = {c for c in classes.values() if c is not INCONCLUSIVE}
        notes = []
        used = list(per_exp)
        if len(said) > 1:
            k = INCONCLUSIVE
            cautions.add("conflict")
            notes.append("experiments disagree: " + "; ".join(
                f"{names.get(e, e)} {WORD[c]} ({statistics.mean(per_exp[e]):+.2f} mM)"
                for e, c in sorted(classes.items(), key=lambda ec: ec[0])))
        elif not said and len({_class(statistics.mean(v), lim) for v in per_exp.values()} - {0}) > 1:
            k = INCONCLUSIVE           # none decides, but their means lie beyond the limit on both sides
            cautions.add("conflict")
            notes.append("experiments disagree: " + "; ".join(
                f"{names.get(e, e)} {WORD[_class(statistics.mean(v), lim)]} ({statistics.mean(v):+.2f} mM)"
                for e, v in sorted(per_exp.items())))
        elif not said:
            k = INCONCLUSIVE
            cautions.add("inconclusive")
            notes.append(f"inconclusive: the replicates do not pin the change down beyond the {lim:g} mM limit, "
                         "or inside it (" + "; ".join(f"{v:+.2f}" for v in values[:8])
                         + (" ..." if len(values) > 8 else "") + " mM)")
        else:
            k = said.pop()
            used = [e for e, c in classes.items() if c == k]
            if len(per_exp) > 1:
                units = [statistics.mean(per_exp[e]) for e in used]
                unsure = sorted(names.get(e, e) for e, c in classes.items() if c is INCONCLUSIVE)
                if unsure:
                    cautions.add("experiment_left_out")
                    notes.append("left out as inconclusive: " + ", ".join(unsure))
                sizes = [abs(u) for u in units]
                if k and len(sizes) > 1 and min(sizes) > 0 and max(sizes) / min(sizes) > AMOUNTS_DIFFER:
                    cautions.add("amounts_differ")
                    notes.append("experiments agree, in amounts from " + ", ".join(
                        f"{names.get(e, e)} {statistics.mean(per_exp[e]):+.2f}" for e in sorted(used)) + " mM")
        values = [v for e in used for v in per_exp[e]] if len(used) < len(per_exp) else values
        if noisy and k is not INCONCLUSIVE:
            k = INCONCLUSIVE
            cautions.add("inconclusive")
            notes.append("inconclusive: the changes lie within their series' own scatter, or within what "
                         "evaporation could account for")
        if k in (1, -1) and len(used) < len(per_exp):
            # the amount is over every experiment that does not contradict, not only the deciding ones: the
            # experiments left out are the ones with the smaller effects, so leaving them out inflates it
            units = [statistics.mean(v) for v in per_exp.values()]
            used = list(per_exp)
            values = [r["change"] for r in valued]
            if _class(statistics.mean(units), lim) != k:
                k = INCONCLUSIVE
                cautions.add("inconclusive")
                notes.append("inconclusive: with the experiments left out included, the mean change is within the "
                             f"{lim:g} mM limit")
        afters = [r["after"] for r in valued if r.get("after") is not None and cultures[r["culture"]].experiment
                  in used]
        if k is not INCONCLUSIVE and afters and "growth_rate_boundary" in cautions:
            moved = statistics.mean(afters)
            # judged as any change, on the replicates' mean, and only when it is a real part of the phase's change
            if _class(moved, lim) and (key[2] == "exponential"
                                       or abs(moved) >= STILL_SHARE * abs(statistics.mean(units))):
                cautions.add("still_changing")
                notes.append(f"still changing by {moved:+.2f} mM between the end of growth by the growth rate and "
                             "where the 90% rule would have ended it"
                             + (" (that change is in the stationary value, not this one)"
                                if key[2] == "exponential" else ""))
        mean = statistics.mean(units)
        sd = statistics.stdev(units) if len(units) > 1 else None
        # bounds on the change: the range of the replicates behind the value (Karoline, 2026-10-07: "bounds
        # resulting from their entries", then "Range of replicates" over an interval on the mean, which was
        # narrower than the replicates and on two experiments absurdly wide); none on one replicate
        spread = (min(values), max(values)) if len(values) > 1 else (None, None)
        # the interval each experiment covered, counted once per experiment as the value is
        by_exp = defaultdict(list)
        for r in valued:
            if cultures[r["culture"]].experiment in used:
                by_exp[cultures[r["culture"]].experiment].append(r)
        start = statistics.mean(statistics.mean(r["start"] for r in rs) for rs in by_exp.values())
        end = statistics.mean(statistics.mean(r["end"] for r in rs) for rs in by_exp.values())
        test = paired(units, [0.0] * len(units))
        if len(values) == 1:
            cautions.add("single_replicate")
        starts = [r["start"] for r in valued if r.get("start") is not None]
        lengths = [r["end"] - r["start"] for r in valued if r.get("start") is not None and r.get("end") is not None]
        if starts and lengths and max(starts) - min(starts) > START_SHARE * statistics.mean(lengths):
            # replicates whose metabolite series start later cover less of the phase (E. coli LF82: two of five
            # replicates have no 0 h sample, so their exponential values run from 4 h; a final review round)
            cautions.add("start_differs")
        if k in (1, -1, 0) and agree and any(len(per_exp[e]) == 2 for e in used if classes.get(e) == k):
            # decided on a pair, which errs more readily than three or more replicates, a change or a zero
            # alike (a pair inside the limit is a false zero more often than a triplicate)
            cautions.add("pair_decided")
        if test and test.get("no_variance"):
            cautions.add("no_variance")
        if whole_note:
            notes.append(whole_note)
        cells[key] = {"mean": mean, "sd": sd, "n": len(values), "n_experiments": len(used), "values": values,
                      "whole_run": "whole_run" in cautions,
                      "direction": PRODUCED if k == 1 else CONSUMED if k == -1 else None, "state": STATE[k],
                      "experiments": experiments, "studies": studies, "media": media, "cautions": sorted(cautions),
                      "notes": notes, "start": start, "end": end,
                      "p_value": test["p"] if test else None,
                      "low": spread[0], "high": spread[1],
                      "initial": statistics.mean(r["initial"] for r in valued), "exponential_h": exponential,
                      "limit": lim}
    return cells


def presence(rows, cultures, limit: float = DETECTION_LIMIT, agree: bool = True, limits: dict | None = None) -> dict:
    """{(taxon, metabolite, phase): {direction: [{"medium", "mean", "n", "studies", "experiments"}]}} from
    the rows outside the value set: per medium, the direction its pooled change shows, judged as a value
    cell is (`pool`); an inconclusive or conflicting medium shows none."""
    by_medium = defaultdict(list)
    for r in rows:
        if r["change"] is None:
            continue
        c = cultures[r["culture"]]
        by_medium[(r["taxon"], r["metabolite"], r["phase"], c.medium_key)].append(r)
    out = defaultdict(lambda: defaultdict(list))
    for (taxon, met, ph, _), members in by_medium.items():
        cell = pool(members, cultures, limit, agree, limits)[(taxon, met, ph)]
        if not cell["direction"]:
            continue
        out[(taxon, met, ph)][cell["direction"]].append({
            "medium": cultures[members[0]["culture"]].medium, "mean": cell["mean"], "n": cell["n"],
            "studies": cell["studies"], "experiments": cell["experiments"],
            "exponential_h": cell["exponential_h"], "cautions": cell["cautions"],
            "whole_run": cell.get("whole_run", False)})
    return {k: dict(v) for k, v in out.items()}


def adjust(cells: dict, correction: str = "bh") -> None:
    """Set each cell's q_value: its p_value corrected over every cell given, one family (reported, not
    used to decide). The search corrects the cells that make its arcs: the pooled cells when arcs are
    merged across studies, else every study's cells together."""
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
            "sd": None if booleans else cell["sd"], "n": cell["n"], "n_experiments": cell.get("n_experiments"),
            "p_value": None if booleans else cell.get("p_value"),
            "q_value": None if booleans else cell.get("q_value"),
            "window_start": cell["start"], "window_end": cell["end"], "exponential_h": cell.get("exponential_h"),
            "medium": "; ".join(cell["media"]),
            "study_ids": studies or cell["studies"], "experiments": cell["experiments"],
            "cautions": list(cell["cautions"]), "notes": list(cell["notes"])}


def _presence_arc(taxon, met, ph, direction, entries, value_cell) -> dict:
    cautions = sorted({x for e in entries for x in e["cautions"]} & {
        "short_record", "window_beyond_data", "stationary_not_reached", "single_replicate", "coarse_sampling",
        "boundaries_differ", "phase_from_other_replicates", "no_variance", "amounts_differ", "still_changing",
        "experiment_left_out", "whole_run"})
    notes = [f"seen in {e['medium']} ({e['mean']:+.2f} mM over {e['n']} replicate(s)); another medium than the "
             "values come from, so only its direction counts" for e in entries]
    if value_cell is not None and value_cell["n"]:
        if value_cell.get("state") == "inconclusive":
            notes.append("in the value medium it was inconclusive")
        elif value_cell["direction"] is None:
            cautions.append("not_detected_in_value_medium")
        else:
            notes.append(f"in the value medium it was {value_cell['direction']} instead")
    if entries and all(e.get("whole_run") for e in entries):
        ph = "whole_run"
    return {"taxon": taxon, "metabolite": met, "phase": ph, "direction": direction, "evidence": PRESENCE_ONLY,
            "amount": None, "change": None, "sd": None, "n": sum(e["n"] for e in entries),
            "n_experiments": len({x for e in entries for x in e["experiments"]}), "p_value": None,
            "q_value": None, "window_start": None, "window_end": None,
            "exponential_h": _mean_of(e.get("exponential_h") for e in entries),
            "medium": "; ".join(sorted({e["medium"] for e in entries})),
            "study_ids": sorted({s for e in entries for s in e["studies"]}),
            "experiments": sorted({x for e in entries for x in e["experiments"]}),
            "cautions": sorted(set(cautions)), "notes": notes}


def _arc_phase(ph: str, cell: dict) -> str:
    """An arc's phase: the cell's, or whole_run when its change spans the whole run."""
    return "whole_run" if cell.get("whole_run") else ph


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
                out.append(_measured_arc(taxon, met, _arc_phase(ph, cell), cell, booleans=booleans))
                out[-1]["merged_arcs"] = len(cell["studies"])
    else:
        for study, cells in sorted(per_study_cells.items()):
            for (taxon, met, ph), cell in cells.items():
                if cell["direction"]:
                    out.append(_measured_arc(taxon, met, _arc_phase(ph, cell), cell, [study], booleans))
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


def merge_genus_cells(cells: dict, taxa: dict, limit: float = DETECTION_LIMIT, limits: dict | None = None) -> dict:
    """Cells keyed by genus: the median of the member taxa's mean changes (taxa with a value only), judged
    against the detection limit like any value. Member taxa that disagree (one produced, one consumed, or one
    changed and another did not) make the genus cell inconclusive, with the caution conflict; a member that
    was inconclusive itself is named and does not vote. `n` is the replicates behind the cell, and
    `n_experiments` its experiments; the note names the taxa."""
    limits = limits or {}
    groups = defaultdict(list)
    for (taxon, met, ph), cell in cells.items():
        groups[(genus_of(taxa, taxon), met, ph)].append((taxon, cell))
    out = {}
    for key, members in groups.items():
        lim = limits.get(key[1], limit)
        voting = [(t, c) for t, c in members if c["mean"] is not None and c.get("state") != "inconclusive"]
        means = [c["mean"] for _, c in voting]
        mean = statistics.median(means) if means else None
        cautions = set().union(*(c["cautions"] for _, c in members)) - {"conflict", "inconclusive"}
        notes = [f"median over {len(means)} taxon(s): " + ", ".join(sorted(taxa[t]["name"] for t, _ in voting))]
        states = {c["state"] for _, c in voting}
        k = _class(mean, lim) if means else INCONCLUSIVE
        if len(states) > 1:
            k = INCONCLUSIVE
            cautions.add("conflict")
            notes.append("taxa disagree: " + "; ".join(f"{taxa[t]['name']} {WORD[_class(c['mean'], lim)]} "
                                                       f"({c['mean']:+.2f} mM)" for t, c in sorted(voting)))
        unsure = sorted(taxa[t]["name"] for t, c in members if c.get("state") == "inconclusive")
        if unsure:
            notes.append("inconclusive, so not counted: " + ", ".join(unsure))
            if not voting:
                cautions.add("inconclusive")
        merged = {"mean": mean, "sd": None, "n": sum(c["n"] for _, c in members),
                  "n_experiments": sum(c.get("n_experiments") or 0 for _, c in members), "values": [],
                  "direction": PRODUCED if k == 1 else CONSUMED if k == -1 else None,
                  "state": STATE[k] if means or unsure else None,
                  "experiments": sorted({x for _, c in members for x in c["experiments"]}),
                  "studies": sorted({x for _, c in members for x in c["studies"]}),
                  "media": sorted({x for _, c in members for x in c["media"]}), "cautions": sorted(cautions),
                  "notes": notes,
                  "start": _mean_of(c["start"] for _, c in members), "end": _mean_of(c["end"] for _, c in members),
                  "p_value": None, "initial": _mean_of(c["initial"] for _, c in members),
                  "exponential_h": _mean_of(c.get("exponential_h") for _, c in members),
                  "merged_taxa": sorted(t for t, _ in members), "limit": lim}
        if not any(c["n"] for _, c in members):
            merged["state"] = None
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
