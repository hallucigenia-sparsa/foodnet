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


def _rate_values(cultures, method: str, spike_factor: float, skipped: list) -> list:
    """[(rate, first hour, last hour the fit used)] of the cultures' growth curves."""
    feature = rates.fitted(method)
    values = []
    for c in cultures:
        if c.growth is None:
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
            values.append((value, start, end))
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
        source, used = FROM_METABOLITE_REPLICATES, by_taxon.get(taxon, [])
        if not values:
            values = _rate_values(others.get(taxon, []), method, spike_factor, skipped)
            source, used = FROM_SAME_MEDIUM, others.get(taxon, [])
        if values:
            # the hours the fits used, the replicates' mean, to set beside the values' intervals (Karoline,
            # 2026-10-22: "add fit_start_h and fit_end_h")
            found[taxon] = {"rate": statistics.median(v for v, _, _ in values), "unit": "1/h", "n": len(values),
                            "source": source, "studies": sorted({c.study for c in used}), "method": method,
                            "fit_start_h": statistics.mean(a for _, a, _ in values),
                            "fit_end_h": statistics.mean(b for _, _, b in values)}
        elif by_taxon.get(taxon) or others.get(taxon):
            missing[taxon] = "its growth curves in this medium give no rate (see the report)"
        else:
            missing[taxon] = "no growth curve of it in this medium"
    return found, missing, skipped


def _declines(cultures, spike_factor: float, skipped: list) -> list:
    found = []
    for c in cultures:
        if c.growth is None or spike(c.growth["values"], spike_factor):
            continue                                    # a spiking curve is already named for its growth rate
        try:
            found.append({**rates.decline(c.growth["times"], c.growth["values"]), "study": c.study})
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no decline rate: {e}"))
    return found


def summarize_declines(found: list, source: str) -> dict:
    """One taxon's decline rate from its replicates' (`rates.decline`): the mean, its scatter (standard deviation
    and standard error) and n, so a real drift can be told from noise; `phase` is "death" where most replicates
    leave their plateau before the last sample, else "plateau" (a lower bound); the window is the replicates'
    mean peak and last sample."""
    values = [f["rate"] for f in found]
    sd = statistics.stdev(values) if len(values) > 1 else None
    deaths = sum(1 for f in found if f["death"])
    return {"rate": statistics.mean(values), "unit": "1/h", "sd": sd,
            "se": None if sd is None else sd / math.sqrt(len(values)), "n": len(values),
            "phase": "death" if 2 * deaths > len(values) else "plateau", "death_replicates": deaths,
            "start_h": statistics.mean(f["start"] for f in found), "end_h": statistics.mean(f["end"] for f in found),
            "source": source, "studies": sorted({f["study"] for f in found}), "replicates": found}


def decline_rates(taxa, value_cultures, other_cultures, value_keys, ignore_media: bool = False,
                  spike_factor: float = SPIKE_FACTOR) -> tuple:
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
        declines, source = _declines(by_taxon.get(taxon, []), spike_factor, skipped), FROM_METABOLITE_REPLICATES
        if not declines:
            declines, source = _declines(others.get(taxon, []), spike_factor, skipped), FROM_SAME_MEDIUM
        if declines:
            found[taxon] = summarize_declines(declines, source)
        elif by_taxon.get(taxon) or others.get(taxon):
            missing[taxon] = "its growth curves in this medium give no usable decline (see the report)"
        else:
            missing[taxon] = "no growth curve of it in this medium"
    return found, missing, skipped


def _curve_unit(growth: dict) -> str:
    technique, curve_unit = growth.get("technique") or "", growth.get("unit") or ""
    return technique if curve_unit.casefold() in ("", technique.casefold()) else f"{technique} {curve_unit}"


def _mean_biomass(times, values, start: float, end: float) -> float:
    """The curve's time average from `start` to `end` (both sample times), by the trapezoid rule."""
    pts = [(t, v) for t, v in zip(times, values, strict=True) if start <= t <= end]
    area = sum((t1 - t0) * (v0 + v1) / 2 for (t0, v0), (t1, v1) in zip(pts, pts[1:], strict=False))
    return area / (end - start)


def maintenance(value_cultures, limit: float, limits: dict | None = None,
                spike_factor: float = SPIKE_FACTOR) -> tuple:
    """({taxon: {metabolite: summary}}, skipped): the Pirt maintenance coefficient, a separate characteristic
    for a fit to use or not (Karoline, 2026-10-09): what each culture took up of a compound per unit of biomass
    per hour while it sat on its growth plateau (`rates.plateau`), where it no longer grew. Per culture: the
    compound's fall over the plateau (a change within its detection limit is 0, "zeros are zeros"), over the
    curve's mean biomass there and the plateau's hours; positive is uptake, negative release. Only a compound
    sampled at or before the plateau's start and at or after its end counts. Per taxon and compound: the mean
    over its cultures in the unit most of them use, with standard deviation, standard error, n and the mean
    window; the biomass is in that curve unit, so the coefficient is mM per unit per hour."""
    skipped, found = [], defaultdict(lambda: defaultdict(list))
    for c in value_cultures:
        if c.growth is None or spike(c.growth["values"], spike_factor):
            continue
        try:
            start, end = rates.plateau(c.growth["times"], c.growth["values"])
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no maintenance coefficient: {e}"))
            continue
        biomass = _mean_biomass(c.growth["times"], c.growth["values"], start, end)
        if biomass <= 0:
            continue
        for mid, met in c.metabolites.items():
            series = met["series"]
            if not series or series[0][0] > start or series[-1][0] < end:
                continue
            fall = value_at(series, start)[0] - value_at(series, end)[0]
            if abs(fall) < (limits or {}).get(mid, limit):
                fall = 0.0
            found[c.taxon["id"]][mid].append({"value": fall / biomass / (end - start), "unit": _curve_unit(c.growth),
                                              "start": start, "end": end, "study": c.study})
    return summarize_maintenance(found), skipped


def summarize_maintenance(found: dict) -> dict:
    """{taxon: {metabolite: summary}} from {taxon: {metabolite: [per-culture entries]}} (see `maintenance`)."""
    out = {}
    for taxon, by_met in found.items():
        for mid, entries in by_met.items():
            units = [e["unit"] for e in entries]
            unit = max(set(units), key=units.count)
            same = [e for e in entries if e["unit"] == unit]
            values = [e["value"] for e in same]
            sd = statistics.stdev(values) if len(values) > 1 else None
            out.setdefault(taxon, {})[mid] = {
                "coefficient": statistics.mean(values), "unit": f"mM per {unit or 'biomass unit'} per h",
                "biomass_unit": unit, "sd": sd, "se": None if sd is None else sd / math.sqrt(len(values)),
                "n": len(values), "start_h": statistics.mean(e["start"] for e in same),
                "end_h": statistics.mean(e["end"] for e in same), "studies": sorted({e["study"] for e in same}),
                "cultures": entries}
    return out


def initial_concentrations(value_cultures) -> dict:
    """{metabolite: {"name", "mean", "min", "max", "n"}}: the concentration at each value-set culture's first
    metabolite sample, in mM."""
    found = defaultdict(list)
    names = {}
    for c in value_cultures:
        for mid, met in c.metabolites.items():
            found[mid].append(met["series"][0][1])
            names[mid] = met["name"]
    return {mid: {"name": names[mid], "mean": statistics.mean(v), "min": min(v), "max": max(v), "n": len(v)}
            for mid, v in found.items()}


def biomass_changes(value_cultures, rows, phase: str) -> dict:
    """{taxon: {"change", "start", "unit", "hours", "n"}}: the growth of each taxon over `phase`, from the
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
        windows[culture] = max(set(found), key=lambda w: (found.count(w), w[1] - w[0]))
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
        technique, curve_unit = c.growth.get("technique") or "", c.growth.get("unit") or ""
        unit = technique if curve_unit.casefold() in ("", technique.casefold()) else f"{technique} {curve_unit}"
        by_taxon[c.taxon["id"]].append({"unit": unit, "experiment": c.experiment,
                                        "start": x0, "change": x1 - x0, "hours": end - start})
    out = {}
    for taxon, found in by_taxon.items():
        units = [f["unit"] for f in found]
        unit = max(sorted(set(units)), key=units.count)
        same = [f for f in found if f["unit"] == unit]

        def per_experiment(key, same=same):
            # each experiment counted once, as the values and their intervals are
            by_exp = {}
            for f in same:
                by_exp.setdefault(f["experiment"], []).append(f[key])
            return statistics.mean(statistics.mean(v) for v in by_exp.values())
        start, change, hours = per_experiment("start"), per_experiment("change"), per_experiment("hours")
        out[taxon] = {"change": change, "start": start, "unit": unit, "hours": hours, "n": len(same),
                      # the mean specific growth rate over the phase: a maximum rate below it cannot be right
                      "phase_rate": (math.log((start + change) / start) / hours
                                     if start > 0 and start + change > 0 and hours > 0 else None)}
    return out


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
