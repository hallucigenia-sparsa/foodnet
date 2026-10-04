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
"""
from __future__ import annotations

import statistics
from collections import defaultdict

from . import rates
from .phase import SPIKE_FACTOR, spike

FROM_METABOLITE_REPLICATES = "the replicates with metabolite data"
FROM_SAME_MEDIUM = "another monoculture in the same medium"


def _rate_values(cultures, method: str, spike_factor: float, skipped: list) -> list:
    feature = rates.feature(method)
    values = []
    for c in cultures:
        if c.growth is None:
            continue
        if spike(c.growth["values"], spike_factor):
            skipped.append((c.label, "growth curve spikes; no growth rate from it"))
            continue
        try:
            value = feature(c.growth["times"], c.growth["values"])
        except rates.RateUnavailable as e:
            skipped.append((c.label, f"no growth rate: {e}"))
            continue
        if value > 0:
            values.append(value)
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
            found[taxon] = {"rate": statistics.median(values), "unit": "1/h", "n": len(values), "source": source,
                            "studies": sorted({c.study for c in used}), "method": method}
        elif by_taxon.get(taxon) or others.get(taxon):
            missing[taxon] = "its growth curves in this medium give no rate (see the report)"
        else:
            missing[taxon] = "no growth curve of it in this medium"
    return found, missing, skipped


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
