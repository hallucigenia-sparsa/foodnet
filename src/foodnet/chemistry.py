"""The chemistry of each resource of a CRM, and each taxon's electron balance (Karoline, 2026-10-08).

**Formula and charge** come from ChEBI, by the id foodnet keys the compound under (foodnet.compounds). mGrowthDB
imports its metabolites from ChEBI but keeps neither (Karoline: "go ahead with ChEBI for now", until mGrowthDB
does). A ChEBI entry that is a class has no formula: succinate (CHEBI:26806) covers its 1- and 2- forms,
fructose (CHEBI:28757) every form of it. Its formula is then taken from a joined form (foodnet.compounds.
EQUIVALENT: succinic acid for succinate) or from FORMULA_FROM, and the resource says which. Nothing else is
filled in: a resource whose formula cannot be read stays blank, with the reason.

**Degree of reduction** (Roels, Energetics and Kinetics in Biotechnology, 1983), electrons per molecule
relative to CO2, H2O, NH3, H2SO4, H3PO4 and H+:

    gamma = 4 C + H - 2 O - 3 N + 6 S + 5 P - charge

The charge term makes it the same for every protonation state: acetate (C2H3O2, -1) and acetic acid (C2H4O2)
are both 8, as they must be, since foodnet joins them into one pool. Per C-mol it is gamma / C. A formula
with another element (Na, Cl, a metal), a variable part (n, R, X) or no charge is not evaluated.

**Electron balance**: per taxon, the electrons in its measured by-products over the electrons in what it
consumed of the measured resources, sum(produced * gamma) / sum(consumed * gamma), in the CRM's mM. It is a
check and changes no value. The electrons a culture takes up go to its biomass, to its by-products and to
compounds nobody measured (dihydrogen, a gas, is often one), so a share well below 1 is expected; a share
above 1 means more electrons out than in, from substrates nobody measured (peptides and amino acids of a
rich medium) or from a measurement problem.
"""
from __future__ import annotations

import datetime
import json
import re
import threading
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from . import __version__
from .compounds import EQUIVALENT

CHEBI_API = "https://www.ebi.ac.uk/chebi/backend/api/public/compound/{}/"
SOURCE = "ChEBI (www.ebi.ac.uk/chebi)"
REFERENCE = ("degree of reduction per molecule relative to CO2, H2O, NH3, H2SO4, H3PO4 and H+: 4C + H - 2O - 3N "
             "+ 6S + 5P - charge (Roels 1983)")
ELECTRONS = {"C": 4, "H": 1, "O": -2, "N": -3, "S": 6, "P": 5}
WORKERS = 6
TIMEOUT = 15

# a ChEBI class without a formula -> (the entry whose formula it takes, why); each checked against ChEBI
# (2026-10-08). A form joined in compounds.EQUIVALENT is found without being listed here.
FORMULA_FROM = {
    28757: (15824, "fructose is a class in ChEBI, with no formula; D-fructose's (C6H12O6) is every form's"),
    33984: (18287, "fucose is a class in ChEBI, with no formula; L-fucose's (C6H12O5) is every form's"),
}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(formula: str) -> dict | None:
    """{element: count} of a plain formula (C2H3O2), or None for one with a variable or repeated part."""
    if not formula or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", formula):
        return None
    counts = {}
    for element, n in _TOKEN.findall(formula):
        counts[element] = counts.get(element, 0) + (int(n) if n else 1)
    return counts


def degree_of_reduction(formula: str, charge) -> tuple:
    """(gamma, carbon atoms, None) or (None, carbon atoms or None, why not)."""
    counts = parse_formula(formula)
    if counts is None:
        return None, None, f"formula {formula!r} has a variable or repeated part"
    carbon = counts.get("C", 0)
    other = sorted(set(counts) - set(ELECTRONS))
    if other:
        return None, carbon, "formula holds " + ", ".join(other) + ", which the degree of reduction here does not cover"
    if charge is None:
        return None, carbon, "ChEBI gives no charge"
    return sum(ELECTRONS[e] * n for e, n in counts.items()) - int(charge), carbon, None


def fetch_compound(chebi_id) -> dict | None:
    """ChEBI's record of one id: {"name", "formula", "charge", "mass"}; None when ChEBI does not hold the id.
    Raises OSError when ChEBI cannot be reached."""
    request = urllib.request.Request(CHEBI_API.format(int(chebi_id)),
                                     headers={"Accept": "application/json", "User-Agent": f"foodnet/{__version__}"})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as r:
            record = json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise OSError(f"ChEBI returned HTTP {e.code}") from None
    except ValueError as e:
        raise OSError(f"ChEBI's answer was not JSON ({e})") from None
    data = record.get("chemical_data") or {}
    return {"name": record.get("name") or "", "formula": data.get("formula") or "", "charge": data.get("charge"),
            "mass": data.get("mass")}


class Memo:
    """fetch_compound, each id once for the life of the object (a client keeps one; nothing is stored on disk)."""

    def __init__(self, fetch=fetch_compound):
        self._fetch, self._seen, self._lock = fetch, {}, threading.Lock()

    def __call__(self, chebi_id):
        key = int(chebi_id)
        with self._lock:
            if key in self._seen:
                return self._seen[key]
        record = self._fetch(key)
        with self._lock:
            self._seen[key] = record
        return record


def _blank(chebi_id: str, note: str) -> dict:
    return {"chebi_id": chebi_id, "formula": "", "charge": None, "carbon": None, "degree_of_reduction": None,
            "per_cmol": None, "formula_from": "", "note": note}


def _one(chebi_id: str, fetch) -> dict:
    """What ChEBI gives one resource, with the formula of a joined form or of FORMULA_FROM for a class."""
    record = fetch(chebi_id)
    if record is None:
        return _blank(chebi_id, "ChEBI does not hold this id")
    source, why = "", ""
    if not record["formula"]:
        key = int(chebi_id)
        candidates = ([FORMULA_FROM[key]] if key in FORMULA_FROM else []) + \
            [(acid, f"ChEBI gives no formula for {record['name'] or 'this entry'} (a class); its form "
                    f"{{name}} (CHEBI:{acid}) has the same degree of reduction")
             for acid, (target, _) in sorted(EQUIVALENT.items()) if target == key and acid != key]
        for other, reason in candidates:
            found = fetch(other)
            if found and found["formula"]:
                record, source, why = found, f"CHEBI:{other}", reason.replace("{name}", found["name"])
                break
        else:
            return _blank(chebi_id, f"ChEBI gives no formula for {record['name'] or 'this entry'}")
    gamma, carbon, problem = degree_of_reduction(record["formula"], record["charge"])
    return {"chebi_id": chebi_id, "formula": record["formula"], "charge": record["charge"], "carbon": carbon,
            "degree_of_reduction": gamma,
            "per_cmol": gamma / carbon if gamma is not None and carbon else None,
            "formula_from": source, "note": problem or why}


def resolve(resources: dict, fetch) -> tuple:
    """({resource id: chemistry}, problem or None) for {resource id: ChEBI id or ""}.

    `fetch` reads one ChEBI id (the client's `chebi_compound`); None means no lookup is possible. When ChEBI
    cannot be reached every resource is blank and the problem says so; a resource without a ChEBI id is blank
    with that reason."""
    out = {mid: _blank("", "mGrowthDB gives no ChEBI id") for mid, chebi in resources.items() if not chebi}
    wanted = [(mid, str(chebi)) for mid, chebi in resources.items() if chebi]
    if not wanted:
        return out, None
    if fetch is None:
        out.update({mid: _blank(chebi, "no ChEBI lookup in this run") for mid, chebi in wanted})
        return out, "no ChEBI lookup in this run"
    failed = []

    def run(item):
        mid, chebi = item
        try:
            return mid, _one(chebi, fetch)
        except (OSError, ValueError) as e:
            failed.append(str(e))
            return mid, _blank(chebi, "ChEBI could not be reached")
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        out.update(dict(pool.map(run, wanted)))
    return out, (f"ChEBI could not be reached for {len(failed)} of {len(wanted)} resource(s) ({failed[0]})"
                 if failed else None)


def provenance(problem) -> dict:
    return {"source": SOURCE, "api": CHEBI_API.format("<id>"), "retrieved_at": datetime.date.today().isoformat(),
            "reference": REFERENCE, "problem": problem}


def electron_balance(consumed, produced, gammas, counted, bounds=None) -> dict:
    """One taxon's electron balance over one phase, from its rows of the consumed and produced matrices (mM,
    None for no number), the resources' degrees of reduction and which resources count (False: a compound of
    the second time window, measured over another window). `bounds` is {"consumed": (lower, upper),
    "produced": (lower, upper)}, the taxon's rows of the bound matrices.

    {"consumed_e_mM", "produced_e_mM", "share", "lower", "upper", "not_counted"}: the share None when nothing
    with electrons was consumed; lower and upper from the replicates' bounds (the fewest electrons out over
    the most in, and the reverse), None unless every number counted has bounds and no degree of reduction is
    negative. not_counted is [(resource index, why)]."""
    sums = {("consumed", "value"): 0.0, ("produced", "value"): 0.0}
    for d in ("consumed", "produced"):
        sums[(d, "lower")] = sums[(d, "upper")] = 0.0
    bounded = bounds is not None
    not_counted = []
    for j, gamma in enumerate(gammas):
        values = {"consumed": consumed[j], "produced": produced[j]}
        if not counted[j]:
            not_counted.append((j, "measured over the second time window"))
            continue
        if all(v is None for v in values.values()):
            not_counted.append((j, "no number"))
            continue
        if gamma is None:
            if any(values.values()):
                not_counted.append((j, "no degree of reduction"))
            continue
        for d, v in values.items():
            if v is None:
                not_counted.append((j, f"{d}: no number"))
                continue
            sums[(d, "value")] += v * gamma
            if not bounded:
                continue
            lo, hi = bounds[d][0][j], bounds[d][1][j]
            if lo is None or hi is None or gamma < 0:
                bounded = False
                continue
            sums[(d, "lower")] += lo * gamma
            sums[(d, "upper")] += hi * gamma
    e_in, e_out = sums[("consumed", "value")], sums[("produced", "value")]
    share = e_out / e_in if e_in > 0 else None
    low = high = None
    if share is not None and bounded:
        low = sums[("produced", "lower")] / sums[("consumed", "upper")] if sums[("consumed", "upper")] > 0 else None
        high = sums[("produced", "upper")] / sums[("consumed", "lower")] if sums[("consumed", "lower")] > 0 else None
    return {"consumed_e_mM": e_in, "produced_e_mM": e_out, "share": share, "lower": low, "upper": high,
            "not_counted": not_counted}
