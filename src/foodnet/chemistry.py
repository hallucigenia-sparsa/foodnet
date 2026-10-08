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
check and changes no value, and it sees one part of each side of the balance (Roels 1983): electrons also come
from substrates nobody measured (peptides and amino acids of a rich medium, dihydrogen or formate taken up by
hydrogen-using taxa) and go to biomass and to by-products nobody measured. So it has no expected side of 1,
a share near 1 is no proof of closure, and 1 minus the share is no biomass yield. Its range is the lowest and
highest culture's own share over the resources the share counts; it is withheld when a resource the taxon
consumed or produced has no degree of reduction, or when its uptake is within the detection limits (see
electron_balance).
"""
from __future__ import annotations

import datetime
import http.client
import json
import re
import threading
import time
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
TIMEOUT = 8
ATTEMPTS = 2           # a dropped or slow answer is tried once more, as mGrowthDB's are (foodnet.mgrowthdb)
_NETWORK = (OSError, http.client.HTTPException)

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
        if formula and "." in formula:
            return None, None, f"formula {formula} is a salt or a mixture of parts"
        return None, None, f"formula {formula!r} has a variable or repeated part"
    carbon = counts.get("C", 0)
    other = sorted(set(counts) - set(ELECTRONS))
    if other:
        return None, carbon, "formula holds " + ", ".join(other) + ", which the degree of reduction here does not cover"
    if charge is None:
        return None, carbon, "ChEBI gives no charge"
    return sum(ELECTRONS[e] * n for e, n in counts.items()) - int(charge), carbon, None


class ChebiUnreachable(OSError):
    """ChEBI could not be reached, or answered in a way foodnet does not read (a moved or changed API)."""


def _plain(name: str) -> str:
    """A ChEBI name without its HTML markup (<small>D</small>-fructose)."""
    return re.sub(r"<[^>]+>", "", name or "")


def _read(url: str) -> dict | None:
    request = urllib.request.Request(url, headers={"Accept": "application/json",
                                                   "User-Agent": f"foodnet/{__version__}"})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # ChEBI says so in JSON when it does not hold an id; a page that is not JSON is a moved API
            try:
                body = json.loads(e.read().decode("utf-8", "replace"))
            except (ValueError, *_NETWORK):
                body = None
            if isinstance(body, dict) and "detail" in body:
                return None
            raise ChebiUnreachable("ChEBI answered 404 without its JSON: its API may have moved") from None
        raise ChebiUnreachable(f"ChEBI returned HTTP {e.code}") from None
    except ValueError as e:
        raise ChebiUnreachable(f"ChEBI's answer was not JSON ({e})") from None
    except _NETWORK as e:
        raise ChebiUnreachable(f"{getattr(e, 'reason', e) or type(e).__name__}") from None


def fetch_compound(chebi_id) -> dict | None:
    """ChEBI's record of one id: {"id", "name", "formula", "charge", "mass"}; None when ChEBI does not hold the
    id. Tried ATTEMPTS times; raises ChebiUnreachable when ChEBI cannot be reached or its answer lacks the
    fields foodnet reads (its API documents no response body, so a change is caught here, not read as a
    compound without a formula)."""
    number = int(str(chebi_id).upper().removeprefix("CHEBI:"))
    last = None
    for attempt in range(ATTEMPTS):
        try:
            record = _read(CHEBI_API.format(number))
            break
        except ChebiUnreachable as e:
            last = e
            if attempt < ATTEMPTS - 1:
                time.sleep(0.5)
    else:
        raise last
    if record is None:
        return None
    if not isinstance(record, dict) or "chemical_data" not in record or "id" not in record:
        raise ChebiUnreachable("ChEBI's answer lacks the fields foodnet reads (chemical_data, id): has its API "
                               "changed?")
    data = record.get("chemical_data") or {}
    return {"id": int(record["id"]), "name": _plain(record.get("ascii_name") or record.get("name")),
            "formula": data.get("formula") or "", "charge": data.get("charge"), "mass": data.get("mass")}


class Memo:
    """fetch_compound, each id once for the life of the object (a client keeps one; nothing is stored on disk).
    A failure is not kept: the next lookup of the id asks again."""

    def __init__(self, fetch=fetch_compound):
        self._fetch, self._seen, self._lock = fetch, {}, threading.Lock()

    def known(self, chebi_id) -> bool:
        """Whether the id was read already, so asking for it sends nothing to ChEBI."""
        with self._lock:
            return int(str(chebi_id).upper().removeprefix("CHEBI:")) in self._seen

    def __call__(self, chebi_id):
        key = int(str(chebi_id).upper().removeprefix("CHEBI:"))
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
        # a secondary id answers with its primary record, which the tables are keyed by
        key = int(record.get("id") or str(chebi_id).upper().removeprefix("CHEBI:"))
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

    `fetch` reads one ChEBI id (the client's `chebi_compound`); None means no lookup is possible. A resource
    ChEBI could not be read for is blank and the problem says how many; a resource without a ChEBI id is blank
    with that reason."""
    out = {mid: _blank("", "mGrowthDB gives no ChEBI id") for mid, chebi in resources.items() if not chebi}
    wanted = [(mid, str(chebi)) for mid, chebi in resources.items() if chebi]
    if not wanted:
        return out, None
    if fetch is None:
        out.update({mid: _blank(chebi, "no ChEBI lookup in this run") for mid, chebi in wanted})
        return out, "no ChEBI lookup in this run"
    failed = []
    # once ChEBI could not be reached, the rest of this call's lookups fail at once, so an unreachable ChEBI
    # costs a search one wait, not one per resource; the next call asks again (reviews: a quiet period over the
    # process, and then over the page's client, which lives for an hour, made reruns fail as the warning advised
    # them)
    down = []
    lock = threading.Lock()

    owner = getattr(fetch, "__self__", fetch)        # the client's Memo, or a Memo itself
    known = getattr(getattr(owner, "_chebi", owner), "known", None)

    def once(chebi_id):
        if known is not None and known(chebi_id):
            return fetch(chebi_id)                   # read before: nothing to wait for (a review)
        with lock:
            if down:
                raise ChebiUnreachable(f"{down[0]} (earlier in this search)")
        try:
            return fetch(chebi_id)
        except ChebiUnreachable as e:
            with lock:
                down.append(str(e))
            raise

    def run(item):
        mid, chebi = item
        try:
            return mid, _one(chebi, once)
        except (*_NETWORK, ValueError) as e:
            failed.append(str(e))
            return mid, _blank(chebi, "ChEBI could not be reached")
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        out.update(dict(pool.map(run, wanted)))
    return out, (f"ChEBI could not be reached for {len(failed)} of {len(wanted)} resource(s) ({failed[0]})"
                 if failed else None)


def provenance(problem) -> dict:
    return {"source": SOURCE, "api": CHEBI_API.format("<id>"), "retrieved_at": datetime.date.today().isoformat(),
            "reference": REFERENCE, "licence": "CC BY 4.0", "problem": problem}


def electron_balance(consumed, produced, gammas, counted, cultures=(), in_medium=None, limits=None,
                     starts=None, names=None) -> dict:
    """One taxon's electron balance over one phase (Karoline, 2026-10-08: "Own gaps + name the rest",
    "Per-culture shares").

    `consumed` and `produced` are its rows of the CRM matrices (mM, None for no number), `gammas` the resources'
    degrees of reduction, `counted` which resources count (False: a compound of the second time window,
    measured over another window), `cultures` each value-medium culture's own changes ([mM per resource, None
    where it measured none]; negative for uptake), `in_medium` which resources the medium held above the
    detection limit when the phase began, `limits` each resource's detection limit (mM) and `starts` each
    resource's concentration when the phase began (mM), the most a resource without a number could hide.

    {"consumed_electrons_mM", "produced_electrons_mM", "share", "share_lower", "share_upper", "share_at_least",
    "cultures", "cultures_left_out", "withheld", "incomplete", "not_counted"}:
      - share: electrons out over electrons in; None when nothing with electrons was consumed, or withheld;
      - withheld: why there is no share although there are numbers: a resource it consumed or produced has no
        degree of reduction, so the share would leave out a substrate (too high) or a product (too low): a gap of
        foodnet's own (Karoline, 2026-10-08: "Withhold for products too"); or the electrons taken up are within
        the detection limits of the resources it counts and of those it has no number for (incomplete): the sum
        of each one's limit times its degree of reduction, the uptake they could hide, so the share would be
        undetermined (Karoline, 2026-10-08: "Withhold, say why"; share_at_least then holds what it is at least:
        "Report 'at least X'", electrons out over electrons in plus the most the resources could hide: for the
        pooled amounts, a counted one its detection limit and one without a number what the medium held when the
        phase began; for a culture that measured all of them, its own changes, each hiding its limit; the lowest
        of these); or the share lies outside its cultures' own range ("Withhold when outside");
      - incomplete: resources the medium held when the phase began, with no number for this taxon (not assayed,
        inconclusive, seen only elsewhere): the share stands, but leaves them out. `in_medium` is per resource,
        at the phase's start (Karoline, 2026-10-08: "At the phase's start");
      - share_lower, share_upper: the lowest and highest culture's own share, over the same resources as the
        share (Karoline, 2026-10-08: "Pooled resources"), over the cultures that measured all of them and took
        up more electrons than those limits ("cultures" says how many, "cultures_left_out" how many measured
        them all but took up no more); None without two such cultures, or with a negative degree of reduction
        counted (an electron acceptor, which makes the share no ratio of electrons out over in);
      - not_counted: [(resource index, why)]."""
    e_in = e_out = 0.0
    used, not_counted, incomplete, missing_gamma = [], [], [], []
    for j, gamma in enumerate(gammas):
        c, p = consumed[j], produced[j]
        if not counted[j]:
            not_counted.append((j, "measured over the second time window"))
            continue
        if c is None and p is None:
            not_counted.append((j, "no number"))
            if in_medium is not None and in_medium[j]:
                incomplete.append(j)
            continue
        if gamma is None:
            if c or p:
                not_counted.append((j, "no degree of reduction"))
                missing_gamma.append(j)
            continue
        for d, v in (("consumed", c), ("produced", p)):
            if v is None:
                not_counted.append((j, f"{d}: no number"))
                if d == "consumed" and in_medium is not None and in_medium[j]:
                    incomplete.append(j)
        e_in += (c or 0.0) * gamma
        e_out += (p or 0.0) * gamma
        used.append(j)
    # the uptake the resources could hide below their detection limits, in electrons: every one counted, and
    # every one the medium held without a number for the taxon (incomplete), charged its limit too (Karoline,
    # 2026-10-08, "Withhold, say why"); the floor below charges an incomplete one more, its phase-start amount,
    # since a bound must hold whatever it hid
    noise = sum((limits[j] if limits else 0.0) * gammas[j] for j in sorted(set(used) | set(incomplete))
                if gammas[j] is not None and gammas[j] > 0)
    withheld = None
    if missing_gamma:
        withheld = ("no degree of reduction for a resource it consumed or produced, so the share would leave it "
                    "out (see not_counted)")
    elif 0 < e_in <= noise:
        withheld = (f"the electrons taken up ({e_in:.3g} mM) are within the detection limits of the resources "
                    f"counted or missing ({noise:.3g} mM), so the share would be undetermined")
    small = withheld is not None and not missing_gamma
    share = e_out / e_in if e_in > 0 and withheld is None else None
    # each culture's own share, over the resources the share counts, from those that measured all of them
    measured = []
    if not missing_gamma and e_in > 0 and all(gammas[j] >= 0 for j in used):
        for changes in cultures:
            if any(changes[j] is None for j in used):
                continue
            measured.append((sum(max(0.0, -changes[j]) * gammas[j] for j in used),
                             sum(max(0.0, changes[j]) * gammas[j] for j in used)))
    shares = [c_out / c_in for c_in, c_out in measured if c_in > noise and c_in > 0]
    left_out = len(measured) - len(shares)
    ranged = len(shares) >= 2
    low, high = (min(shares), max(shares)) if ranged else (None, None)
    if share is not None and ranged and not low - 1e-9 <= share <= high + 1e-9:
        # the share pools every culture, the range only those that measured every resource counted and took up
        # enough to judge: outside it, the others drive it (Karoline, 2026-10-08, "Withhold when outside")
        withheld = (f"the share ({share:.3g}) lies outside its cultures' own ({low:.3g} to {high:.3g}), so cultures "
                    "outside the range (too little uptake to judge, or not every resource counted) drive it")
        share = None
    if small:
        # what the measured numbers still say (Karoline, 2026-10-08, "Report 'at least X'"): a counted resource
        # hides at most its detection limit, one without a number at most what the medium held when the phase
        # began; the lowest of the pooled amounts' and each culture's own (reviews: the limit alone understated
        # what an inconclusive resource hides, and the pooled means hid cultures below them)
        hidden = [j for j in incomplete if gammas[j] is not None and gammas[j] > 0]
        bounded = starts is not None and all(starts[j] is not None for j in hidden) and \
            all(gammas[j] is not None and gammas[j] >= 0 for j in incomplete)
        floor = None
        if bounded:
            def limit(j):
                return (limits[j] if limits else 0.0) * gammas[j]
            # a resource without a number for the taxon that the medium held below the limit hides at most what it
            # held; one nobody measured in the medium is taken as absent, and the reason says so (a review)
            unseen = [j for j in range(len(gammas)) if counted[j] and consumed[j] is None and produced[j] is None
                      and j not in incomplete and gammas[j] is not None and gammas[j] > 0]
            absent = [j for j in unseen if starts[j] is None]
            below = sum(min(starts[j], limits[j] if limits else 0.0) * gammas[j] for j in unseen
                        if starts[j] is not None)
            most = sum(limit(j) for j in used if gammas[j] > 0) + sum(starts[j] * gammas[j] for j in hidden) + below
            floors = [e_out / (e_in + most)]
            # a culture by its own numbers, where it measured every resource counted and every one the pooled
            # values lack, each hiding at most its limit (a review: the pooled phase-start amount bounds no one
            # culture, and a culture's own change is a measured number)
            every = sorted(set(used) | set(hidden))
            for changes in cultures:
                if any(changes[j] is None for j in every):
                    continue
                c_in = sum(max(0.0, -changes[j]) * gammas[j] for j in every)
                c_out = sum(max(0.0, changes[j]) * gammas[j] for j in every)
                floors.append(c_out / (c_in + sum(limit(j) for j in every if gammas[j] > 0) + below))
            floor = min(floors)
            # a second-window compound was measured over its own window, which can reach into this phase: the
            # floor leaves it out, and says so (a review)
            notes = (["leaving out the second time window's compounds"] if not all(counted) else []) + \
                (["taking " + ", ".join(names[j] if names else str(j) for j in absent) + ", which nobody measured "
                  "in this medium, as absent"] if absent else [])
            withheld += f"; it is at least {floor:.3g}" + (" (" + ", and ".join(notes) + ")" if notes else "")
        else:
            withheld += "; a resource without a number has no bound, so no floor either"
    return {"consumed_electrons_mM": e_in, "produced_electrons_mM": e_out, "share": share,
            "share_lower": low, "share_upper": high, "share_at_least": floor if small else None,
            "cultures": len(shares), "cultures_left_out": left_out, "withheld": withheld,
            "incomplete": incomplete, "not_counted": not_counted}
