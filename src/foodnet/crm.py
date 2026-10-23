"""What a consumer-resource model takes besides the matrices: growth rates and initial concentrations.

A consumer-resource model (CRM) needs each taxon's maximum growth rate (Karoline, 2026-10-04: "CRMs also
take growth rates, so these have to be collected and sent as well"). They are taken, in this order:

  1. from the growth curves of the replicates whose metabolites gave the values (the same cultures);
  2. failing that, from another batch monoculture of the taxon in the same medium (one without metabolite
     data, or one outside the replicates used).

A taxon with neither has no rate, and says why; nothing is borrowed from another medium, since a rate is
specific to its medium. The rate is grownet's: easylinear by default, as mGrowthDB computes the rates it
reports, the median over the replicates.

The initial concentrations are the medium's: each metabolite's concentration at the first sample of the
value-set cultures, averaged over them (Karoline, 2026-10-04: "yes, send initial medium concentrations").

The biomass change (`biomass_changes`) is what a consumer-resource model needs to turn the amounts into
yields (Karoline, 2026-10-06, after a review showed that miaSim reads the size of a positive entry of its
efficiency matrix as a yield, so uptake shares could not reproduce a monoculture): each taxon's growth over
the same phase its metabolites were measured over, in the unit of its growth curve, with the phase's length.
"""
from __future__ import annotations

import math
import statistics
from collections import defaultdict

from . import rates
from .phase import SPIKE_FACTOR, spike, value_at

FROM_METABOLITE_REPLICATES = "the replicates with metabolite data"
FROM_SAME_MEDIUM = "another monoculture in the same medium"


def technique(growth: dict) -> str:
    """The measurement a growth curve is: technique, unit, its flow cytometry gate where it has one, and whether
    it counts the strain or the whole culture ("fc Cells/mL live (strain)", "od (culture)"), so rates and biomass
    from different techniques are not mixed up (Karoline, 2026-10-22: "the measurement technique, to prevent
    mismatches across techniques")."""
    level = "strain" if growth.get("level") == "strain" else "culture"
    return f"{_curve_unit(growth)} ({level})"


def _rate_values(cultures, method: str, spike_factor: float, skipped: list) -> list:
    """[(rate, first hour, last hour the fit used, technique, study)] of the cultures' growth curves."""
    feature = rates.fitted(method)
    values = []
    for c in cultures:
        if c.growth is None:
            continue
        if getattr(c, "not_grown", False):
            # as for declines and maintenance (Karoline, 2026-10-22: "growth rates should follow the same rules")
            skipped.append((c.label, "the culture did not grow; no growth rate from it"))
            continue
        if spike(c.growth["values"], spike_factor):
            skipped.append((c.label, "growth curve spikes; no growth rate from it"))
            continue
        try:
            value, start, end = feature(c.growth["times"], c.growth["values"])
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no growth rate: {e}"))
            continue
        if value > 0:
            values.append((value, start, end, technique(c.growth), c.study))
        else:
            skipped.append((c.label, f"non-positive growth rate ({value:g})"))
    return values


def growth_rates(taxa, value_cultures, other_cultures, value_keys, ignore_media: bool = False,
                 rate_method: str = rates.DEFAULT_METHOD, window: int = rates.DEFAULT_WINDOW,
                 spike_factor: float = SPIKE_FACTOR) -> tuple:
    """({taxon: {"rate", "unit", "n", "source", "studies"}}, {taxon: why there is none}, skipped).

    `value_cultures` gave the values; `other_cultures` are every other batch monoculture read (growth-only
    and metabolite ones outside the value set); `value_keys` are the medium keys of the value set."""
    method = rates.method_name(rate_method, window)
    skipped = []
    by_taxon = defaultdict(list)
    for c in value_cultures:
        by_taxon[c.taxon["id"]].append(c)
    others = defaultdict(list)
    for c in other_cultures:
        if ignore_media or c.medium_key in value_keys:
            others[c.taxon["id"]].append(c)
    found, missing = {}, {}
    for taxon in taxa:
        values = _rate_values(by_taxon.get(taxon, []), method, spike_factor, skipped)
        source = FROM_METABOLITE_REPLICATES
        if not values:
            values = _rate_values(others.get(taxon, []), method, spike_factor, skipped)
            source = FROM_SAME_MEDIUM
        if values:
            # the growth curves of the technique and unit most replicates use, as declines, maintenance and
            # biomass take them (Karoline, 2026-10-22: "growth rates should follow the same rules")
            kind = majority([v[3].rsplit(" (", 1)[0] for v in values])
            left_out = sum(1 for v in values if v[3].rsplit(" (", 1)[0] != kind)
            values = [v for v in values if v[3].rsplit(" (", 1)[0] == kind]
            # the hours the fits used, to set beside the values' intervals (Karoline, 2026-10-22: "add
            # fit_start_h and fit_end_h"): the earliest first and the latest last hour of any replicate's fit,
            # so no hour is named that no fit used (a review: means of 2-6 h and 20-30 h gave 11-18 h)
            found[taxon] = {"rate": statistics.median(v[0] for v in values), "unit": "1/h", "n": len(values),
                            "source": source, "studies": sorted({v[4] for v in values}), "method": method,
                            "fit_start_h": min(v[1] for v in values),
                            "fit_end_h": max(v[2] for v in values),
                            "technique": "; ".join(sorted({v[3] for v in values})),
                            "other_techniques_left_out": left_out}
        elif by_taxon.get(taxon) or others.get(taxon):
            missing[taxon] = "its growth curves in this medium give no rate (see the report)"
        else:
            missing[taxon] = "no growth curve of it in this medium"
    return found, missing, skipped


def _declines(cultures, spike_factor: float, skipped: list, fraction: float = 0.9) -> list:
    found = []
    for c in cultures:
        if c.growth is None or spike(c.growth["values"], spike_factor):
            continue                                    # a spiking curve is already named for its growth rate
        if getattr(c, "not_grown", False):
            skipped.append((c.label, "no decline rate: the culture did not grow"))
            continue
        try:
            found.append({**rates.decline(c.growth["times"], c.growth["values"], fraction), "study": c.study,
                          "experiment": c.experiment, "technique": technique(c.growth)})
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no decline rate: {e}"))
    return found


def majority(units: list) -> str:
    """The unit most entries share; a tie goes to the technique foodnet prefers for a growth curve
    (reading.GROWTH_TECHNIQUES), then to the name, so the same search always picks the same (a review: a tie
    went by Python's string hashing, which changes from run to run)."""
    from .reading import GROWTH_TECHNIQUES

    def rank(unit):
        head = unit.split(" ", 1)[0].casefold()
        return GROWTH_TECHNIQUES.index(head) if head in GROWTH_TECHNIQUES else len(GROWTH_TECHNIQUES)
    return min(set(units), key=lambda u: (-units.count(u), rank(u), u))


def summarize_declines(found: list, source: str) -> dict:
    """One taxon's decline rate from its replicates' (`rates.decline`), over the growth curves of the technique
    most of them use (as biomass and maintenance are; a review: an OD and a cell-count decline were averaged):
    the mean, its scatter between replicates (standard deviation, standard error) and n, with the number of
    experiments they come from (Karoline, 2026-10-22: replicates as the unit, labeled, since replicates of one
    experiment share inoculum and sampling); `phase` is "death" where most replicates leave their plateau band
    before the last sample, else "plateau"; the window is the replicates' mean start and last sample."""
    # the technique and its unit decide; whether a curve counts the strain or the culture is a label, since in a
    # monoculture both count the one strain (a review: B. fragilis' culture-level counts were left out of a genus)
    def measure(f):
        return f.get("technique", "").rsplit(" (", 1)[0]
    kind = majority([measure(f) for f in found])
    left_out = [f for f in found if measure(f) != kind]
    found = [f for f in found if measure(f) == kind]
    values = [f["rate"] for f in found]
    sd = statistics.stdev(values) if len(values) > 1 else None
    deaths = sum(1 for f in found if f["death"])
    return {"rate": statistics.mean(values), "unit": "1/h", "sd": sd,
            "se": None if sd is None else sd / math.sqrt(len(values)), "n": len(values),
            "n_experiments": len({f.get("experiment") for f in found}),
            "phase": "death" if 2 * deaths > len(values) else "plateau", "death_replicates": deaths,
            "replicates_to_zero": sum(1 for f in found if f.get("to_zero")),
            "fit_r2": (statistics.median(r2s) if (r2s := [f["r2"] for f in found if f.get("r2") is not None])
                       else None),
            "start_h": statistics.mean(f["start"] for f in found), "end_h": statistics.mean(f["end"] for f in found),
            "source": source, "studies": sorted({f["study"] for f in found}),
            "technique": "; ".join(sorted({f.get("technique", "") for f in found})),
            "other_techniques_left_out": len(left_out), "replicates": found + left_out}


def decline_rates(taxa, value_cultures, other_cultures, value_keys, ignore_media: bool = False,
                  spike_factor: float = SPIKE_FACTOR, fraction: float = 0.9) -> tuple:
    """({taxon: decline rate (summarize_declines)}, {taxon: why there is none}, skipped): the decline after each
    taxon's growth peak, from cultures in the same order as its growth rate (the replicates with metabolite
    data, else another monoculture in the same medium). A taxon without one is absent, never 0."""
    skipped = []
    by_taxon = defaultdict(list)
    for c in value_cultures:
        by_taxon[c.taxon["id"]].append(c)
    others = defaultdict(list)
    for c in other_cultures:
        if ignore_media or c.medium_key in value_keys:
            others[c.taxon["id"]].append(c)
    found, missing = {}, {}
    for taxon in taxa:
        declines = _declines(by_taxon.get(taxon, []), spike_factor, skipped, fraction)
        source = FROM_METABOLITE_REPLICATES
        # another monoculture only where the value replicates have no growth curve: never in place of curves
        # that gave no usable decline (a review: curves falling to zero gave way to a flat one elsewhere)
        if not declines and not any(c.growth for c in by_taxon.get(taxon, [])):
            declines, source = _declines(others.get(taxon, []), spike_factor, skipped, fraction), FROM_SAME_MEDIUM
        if declines:
            found[taxon] = summarize_declines(declines, source)
        elif by_taxon.get(taxon) or others.get(taxon):
            missing[taxon] = "its growth curves in this medium give no usable decline (see the report)"
        else:
            missing[taxon] = "no growth curve of it in this medium"
    return found, missing, skipped


def _curve_unit(growth: dict) -> str:
    """The technique, its unit and its cell gate where it has one ("fc Cells/mL live"): live, total and ungated
    counts are different measures (Karoline, 2026-10-22: "Yes, gate in the key")."""
    technique, curve_unit = growth.get("technique") or "", growth.get("unit") or ""
    named = technique if curve_unit.casefold() in ("", technique.casefold()) else f"{technique} {curve_unit}"
    return f"{named} {growth['gate']}" if growth.get("gate") else named


def _mean_biomass(times, values, start: float, end: float) -> float:
    """The curve's time average from `start` to `end` (both sample times), by the trapezoid rule."""
    pts = [(t, v) for t, v in zip(times, values, strict=True) if start <= t <= end]
    area = sum((t1 - t0) * (v0 + v1) / 2 for (t0, v0), (t1, v1) in zip(pts, pts[1:], strict=False))
    return area / (end - start)


def maintenance(value_cultures, limit: float, limits: dict | None = None,
                spike_factor: float = SPIKE_FACTOR, dropped=frozenset(), fraction: float = 0.9,
                name_units: bool = True, raw: bool = False) -> tuple:
    """({taxon: {metabolite: summary}}, skipped): the Pirt maintenance coefficient, a separate characteristic
    for a fit to use or not (Karoline, 2026-10-09): what each culture took up of a compound per unit of biomass
    per hour while it sat on its growth plateau (`rates.plateau`). Per culture: the compound's fall over the
    plateau (a change within its detection limit is 0, "zeros are zeros"), over the curve's mean biomass there
    and the plateau's hours; positive is uptake, negative release. Only a compound sampled at or before the
    plateau's start and at or after its end counts, never a culture that did not grow, nor an (experiment,
    compound) pair `dropped` as a duplicate deposit (a review: one deposit counted twice). Per taxon: the
    cultures of the unit most of its cultures use (one unit for all its compounds), and per compound the mean
    with standard deviation, standard error, n, the experiments, the window, the fewest growth samples a
    plateau holds, and the compound's median concentration at the plateau's start: where it was below the
    detection limit (below_limit_replicates; used up or never there), a 0 measures no demand (a
    review). On a plateau, growth and death may balance, so the coefficient bounds Pirt's m_s from above."""
    skipped, by_taxon = [], defaultdict(list)
    for c in value_cultures:
        if c.growth is None or getattr(c, "not_grown", False) or spike(c.growth["values"], spike_factor):
            continue
        try:
            start, end = rates.plateau(c.growth["times"], c.growth["values"], fraction)
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no maintenance coefficient: {e}"))
            continue
        biomass = _mean_biomass(c.growth["times"], c.growth["values"], start, end)
        if biomass > 0:
            by_taxon[c.taxon["id"]].append((c, start, end, biomass))
    found = defaultdict(lambda: defaultdict(list))
    for taxon, plateaus in by_taxon.items():
        # every culture with its own unit: summarize_maintenance keeps the taxon's majority unit and counts the rest
        for c, start, end, biomass in plateaus:
            unit = _curve_unit(c.growth)
            samples = sum(1 for t in c.growth["times"] if start <= t <= end)
            curve = list(zip(c.growth["times"], c.growth["values"], strict=True))
            first_x = value_at(curve, start)[0]
            # how much the biomass changed over the window: a "plateau" in a death phase falls (a review)
            change = value_at(curve, end)[0] / first_x - 1 if first_x > 0 else None
            for mid, met in c.metabolites.items():
                series = met["series"]
                if (c.experiment, mid) in dropped or not series or series[0][0] > start or series[-1][0] < end:
                    continue
                at_start = value_at(series, start)[0]
                at_end = value_at(series, end)[0]
                fall = at_start - at_end
                lim = (limits or {}).get(mid, limit)
                if abs(fall) < lim:
                    fall = 0.0
                found[taxon][mid].append({"value": fall / biomass / (end - start), "unit": unit, "start": start,
                                          "end": end, "study": c.study, "experiment": c.experiment,
                                          "culture": c.label, "samples": samples, "at_start": at_start,
                                          "limit": lim, "gone_by_end": at_start >= lim > at_end,
                                          "biomass_relative_change": change})
    if name_units:
        skipped += unit_skips(found)
    if raw:
        # every culture's own entry, for a merge to genus to pool before choosing its unit (a review: a strain's
        # culture whose compound had no entry in the strain's unit never reached the genus)
        return found, skipped
    return summarize_maintenance(found), skipped


def unit_skips(found: dict) -> list:
    """The cultures each taxon (or genus) leaves out of its maintenance for their unit, as report lines (a review:
    with merged genera, a strain's cultures left out at genus level went unnamed)."""
    lines = []
    for by_met in found.values():
        unit = _taxon_unit(by_met)
        for culture, other in sorted({(e["culture"], e["unit"]) for entries in by_met.values() for e in entries}):
            if other != unit:
                lines.append((culture, f"no maintenance coefficient: its growth curve is in {other}, and most of the "
                                       f"cultures of its taxon (or genus) are in {unit}"))
    return lines


def _taxon_unit(by_met: dict) -> str:
    """The unit most of a taxon's cultures use, each culture counted once whatever it measured (a review: counted
    per compound, two OD cultures measuring twelve compounds outvoted three qPCR cultures)."""
    per_culture = {e.get("culture", id(e)): e["unit"] for entries in by_met.values() for e in entries}
    return majority(list(per_culture.values()))


def summarize_maintenance(found: dict) -> dict:
    """{taxon: {metabolite: summary}} from {taxon: {metabolite: [per-culture entries]}} (see `maintenance`)."""
    out = {}
    for taxon, by_met in found.items():
        unit = _taxon_unit(by_met)
        for mid, entries in by_met.items():
            same = [e for e in entries if e["unit"] == unit]
            if not same:
                continue
            values = [e["value"] for e in same]
            sd = statistics.stdev(values) if len(values) > 1 else None
            at_start = statistics.median(e["at_start"] for e in same)
            gone = sum(1 for e in same if e["at_start"] < e["limit"])
            out.setdefault(taxon, {})[mid] = {
                "coefficient": statistics.mean(values), "unit": f"mM h-1 per ({unit or 'biomass unit'})",
                "biomass_unit": unit, "sd": sd, "se": None if sd is None else sd / math.sqrt(len(values)),
                "n": len(values), "n_experiments": len({e["experiment"] for e in same}),
                "start_h": statistics.mean(e["start"] for e in same),
                "end_h": statistics.mean(e["end"] for e in same), "studies": sorted({e["study"] for e in same}),
                "plateau_samples": min(e["samples"] for e in same), "start_mM": at_start,
                # how many replicates had the compound below the detection limit when the plateau began (used
                # up, or never there): their 0 measures no demand (a review: one flag from the mean hid a
                # culture that took it up)
                "below_limit_replicates": gone,
                "below_limit_by_end": sum(1 for e in same if e["gone_by_end"]),
                # end over start, minus 1: not biomass.csv's biomass_change, an amount in the curve's unit
                "biomass_relative_change": (statistics.mean(ch) if (ch := [e["biomass_relative_change"] for e in same
                                                                           if e["biomass_relative_change"] is not None])
                                            else None),
                # each replicate's own coefficient, so a spread set by the window's placement can be seen
                "replicate_values": [e["value"] for e in same],
                # cultures of this compound in another unit than the taxon's (or genus's) one (a review)
                "other_units_left_out": len(entries) - len(same), "cultures": entries}
    return out


def initial_concentrations(value_cultures, dropped=frozenset()) -> dict:
    """{metabolite: {"name", "mean", "min", "max", "n"}}: the concentration at each value-set culture's first
    metabolite sample, in mM; a duplicate deposit's repeated compounds (`dropped`) count once, as in the values."""
    found = defaultdict(list)
    names = {}
    for c in value_cultures:
        for mid, met in c.metabolites.items():
            if (c.experiment, mid) in dropped:
                continue
            found[mid].append(met["series"][0][1])
            names[mid] = met["name"]
    return {mid: {"name": names[mid], "mean": statistics.mean(v), "min": min(v), "max": max(v), "n": len(v)}
            for mid, v in found.items()}


def biomass_changes(value_cultures, rows, phase: str) -> dict:
    """{taxon: {"change", "start", "start_h", "unit", "hours", "n"}}: the growth of each taxon over `phase`, from the
    value-set cultures' growth curves, interpolated at the start and end of the window their metabolite
    values cover (the window most of the culture's compounds share, since they can end at different
    samples), averaged over the cultures of each experiment and then over the experiments. A taxon whose curves come in
    several techniques or units takes the one most of its cultures use (n says how many); one without a
    curve or a window is absent."""
    windows, whole, seen = {}, set(), {}
    for r in rows:
        if r["phase"] == phase and r["change"] is not None and not r.get("second"):
            # the window most of the culture's compounds share, whatever order they come in (a review: formate
            # ending at 20 h and glucose at 120 h in one culture gave one or the other), ties to the longest
            seen.setdefault(r["culture"], []).append((r["start"], r["end"]))
            if "whole_run" in r["cautions"]:
                whole.add(r["culture"])
    for culture, found in seen.items():
        windows[culture] = max(sorted(set(found)), key=lambda w: (found.count(w), w[1] - w[0]))
    # as in the values (foodnet.derive.pool): where a taxon has cultures with a phase, its whole-run cultures
    # are left out, so the biomass comes from the same cultures and window as the amounts (a tenth review
    # round: study SMGDB00000010 deposits R. intestinalis and a Lachnospiraceae strain under one taxon)
    phased = {c.taxon["id"] for index, c in value_cultures if index in windows and index not in whole}
    for index, c in value_cultures:
        if index in whole and c.taxon["id"] in phased:
            windows.pop(index, None)
    by_taxon = defaultdict(list)
    for index, c in value_cultures:
        if c.growth is None or index not in windows:
            continue
        start, end = windows[index]
        series = list(zip(c.growth["times"], c.growth["values"], strict=True))
        x0, _ = value_at(series, start)
        x1, _ = value_at(series, end)
        # the unit is the technique and its unit together: qPCR gene copies and flow cytometry events can both
        # be "Cells/mL" and still measure different things (DNA from lysed cells counts in one, not the other)
        by_taxon[c.taxon["id"]].append({"unit": _curve_unit(c.growth), "experiment": c.experiment,
                                        "start": x0, "change": x1 - x0, "hours": end - start,
                                        "window_start": start, "series": series})
    out = {}
    for taxon, found in by_taxon.items():
        units = [f["unit"] for f in found]
        unit = majority(units)          # a tie goes to the preferred technique, as for declines and maintenance
        same = [f for f in found if f["unit"] == unit]

        def per_experiment(key, same=same):
            # each experiment counted once, as the values and their intervals are
            by_exp = {}
            for f in same:
                by_exp.setdefault(f["experiment"], []).append(f[key])
            return statistics.mean(statistics.mean(v) for v in by_exp.values())
        start, change, hours = per_experiment("start"), per_experiment("change"), per_experiment("hours")
        # the starting abundance a simulation takes: every culture read at one hour, the median of its window's
        # start, and averaged on the log scale (Karoline, 2026-10-23: "fix biomass_start first"; a review: an
        # arithmetic mean of abundances read at 0 h and at 4 h gave E. coli LF82 1.26e7 where its cultures held
        # 2.9e5 at 0 h, so a simulation started it about 40 times too high)
        # the median over experiments, each its own median, as the logs are averaged (a review)
        by_exp_starts = {}
        for f in same:
            by_exp_starts.setdefault(f["experiment"], []).append(f["window_start"])
        at_h = statistics.median(statistics.median(v) for v in by_exp_starts.values())
        logs, ends, end_logs, exact = {}, {}, {}, set()
        for f in same:
            x = _log_at(f["series"], at_h)
            y = _log_at(f["series"], f["window_start"] + f["hours"])
            if x is not None and y is not None:
                logs.setdefault(f["experiment"], []).append(x)
                ends.setdefault(f["experiment"], []).append(f["window_start"] + f["hours"])
                end_logs.setdefault(f["experiment"], []).append(y)
                exact |= {v for u, v in f["series"] if u == at_h and v > 0}
        flat = [x for v in logs.values() for x in v]
        if not flat:
            common = None
        elif max(flat) == min(flat) and len(exact) == 1:
            common = exact.pop()        # one sampled value: itself, not exp(log(x)) (a review: 999.9999 for 1000)
        else:
            common = math.exp(statistics.mean(statistics.mean(v) for v in logs.values()))
        # the phase's length from that hour to the cultures' mean window end (a review: the start had moved to 0 h
        # while the length still came from starts at 0 h and at 4 h); without a start, no such length
        hours = (statistics.mean(statistics.mean(v) for v in ends.values()) - at_h) if ends else None
        # the rate from biomass_start to the cultures' mean abundance at their windows' ends, both as logs averaged
        # alike, over that length: one trajectory from the measured start to the measured end (reviews: averaging
        # slopes and hours apart overstated B. thetaiotaomicron's growth by 19%; the start plus the per-window
        # change ended stationary phases far from the measured end, once below zero)
        end_log = statistics.mean(statistics.mean(v) for v in end_logs.values()) if end_logs else None
        start_log = statistics.mean(statistics.mean(v) for v in logs.values()) if logs else None
        phase_rate = ((end_log - start_log) / hours
                      if end_log is not None and start_log is not None and hours and hours > 0 else None)
        out[taxon] = {"change": change, "start": common, "start_h": at_h, "unit": unit, "hours": hours,
                      # the abundance at the windows' own starts, which the change runs from
                      "window_start_abundance": start,
                      "n": len(same),
                      # cultures whose curve gives nothing at biomass_start_h (a 0 there, or no sample by then)
                      "start_left_out": len(same) - sum(len(v) for v in logs.values()),
                      # cultures of another unit, left out (a review: B. thetaiotaomicron's ungated copy, whose
                      # values the amounts hold, was left out of its biomass without a word)
                      "other_units_left_out": len(found) - len(same),
                      # the mean specific growth rate over the phase: a maximum rate below it cannot be right
                      "phase_rate": phase_rate}
    return out


def _log_at(series, t: float):
    """The natural log of a growth curve's abundance at hour t, interpolated on the log scale between the samples
    around it (after the last sample, the last); None where it reads 0 at t, or has no positive sample at or before
    t (a measured 0 has no log, and a later sample does not replace it; a review)."""
    if any(u == t and v <= 0 for u, v in series):
        return None
    positive = [(u, math.log(v)) for u, v in series if v > 0]
    if not positive or t < positive[0][0]:
        return None           # nothing at or before t: a later, denser sample does not stand in (a review)
    if t == positive[0][0]:
        return positive[0][1]
    if t >= positive[-1][0]:
        return positive[-1][1]
    for (t0, y0), (t1, y1) in zip(positive, positive[1:], strict=False):
        if t0 <= t <= t1:
            return y0 if t1 == t0 else y0 + (y1 - y0) * (t - t0) / (t1 - t0)
    return positive[-1][1]


def culture_changes(rows, cultures, limit: float, limits: dict | None = None) -> dict:
    """{(taxon, phase): {culture label: {metabolite: change in mM}}}: each value-medium culture's own changes,
    from which the electron balance takes its range (foodnet.chemistry.electron_balance: each culture's share
    over the resources the pooled share counts; Karoline, 2026-10-08, "Per-culture shares", "Pooled
    resources"). As in the values
    (foodnet.derive.pool), a whole-run change is left out where another culture of the taxon has a phase value
    for the compound; second-window changes are left out, since the balance does not count them. A change inside
    the detection limit (`limits` per compound, else `limit`) is 0, as in the matrices ("zeros are zeros")."""
    groups = defaultdict(list)
    for r in rows:
        if r["change"] is not None and not r.get("second"):
            groups[(r["taxon"], r["metabolite"], r["phase"])].append(r)
    out = defaultdict(dict)
    for (taxon, met, phase), members in groups.items():
        if any("whole_run" not in r["cautions"] for r in members):
            members = [r for r in members if "whole_run" not in r["cautions"]]
        for r in members:
            change = r["change"] if abs(r["change"]) >= (limits or {}).get(met, limit) else 0.0
            out[(taxon, phase)].setdefault(cultures[r["culture"]].label, {})[met] = change
    return dict(out)
