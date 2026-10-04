"""Batch monocultures with metabolite data, read from mGrowthDB.

foodnet works with batch monocultures (Karoline, 2026-10-04): one strain, in batch, with metabolites
measured. Each independent bioreplicate becomes one `Culture`, holding a growth curve (for the phase
boundary and the growth rate) and one series per metabolite. What cannot be used is reported with its reason,
as in grownet; nothing is guessed.

Rules worth stating, because each is a choice:

  * A bioreplicate flagged `isAverage` is the mean of the real replicates (mGrowthDB) and is left out: it
    is not an independent replicate.
  * The growth curve: a per-strain measurement first, then the culture-level one, which in a monoculture
    measures that one strain; within each, cell counts before optical density (flow cytometry, qPCR,
    plate counts, 16S, then OD), since OD saturates first and so ends exponential growth early. When one
    replicate holds several curves of the same kind (study SMGDB00000007 records three culture-level flow
    cytometry traces), the per-strain one decides.
  * Times are converted to hours, concentrations to mM (`foodnet.compounds`). A series in a unit that is
    not a concentration is left out and reported.
  * Experiments that are not batch (chemostat, serial dilution) are left out and reported: production and
    consumption in a diluted culture are not net changes in the vessel.
  * Monocultures without any metabolite are kept apart, for one thing only: a growth rate in the same
    medium when the metabolite replicates give none (Karoline's order, 2026-10-04).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import compounds
from .model import genus_species
from .phase import hours

BATCH = "batch"
AVERAGE = "average of the replicates, not an independent replicate"
UNREAD = "could not read"
# growth curve techniques, most preferred first (mGrowthDB techniqueType)
GROWTH_TECHNIQUES = ("fc", "qpcr", "plates", "cfu", "16s", "od")
# what is never a growth curve
NOT_GROWTH = ("metabolite", "ph")
MIN_POINTS = 2
NCBI, BY_NAME = "ncbi", "name"


@dataclass
class Culture:
    """One bioreplicate of a batch monoculture: its growth curve and its metabolite series."""

    taxon: dict                    # {"id", "name", "taxon_id", "species", "identity"}
    study: str
    experiment: str
    experiment_name: str
    description: str
    medium: str                    # as mGrowthDB names it
    medium_key: str                # the name normalized, for grouping (`medium_key`)
    replicate: str
    growth: dict | None = None     # {"times", "values", "technique", "level", "unit"}, times in hours
    metabolites: dict = field(default_factory=dict)   # node id -> {"name", "chebi_id", "series", "recorded_as"}

    @property
    def label(self) -> str:
        return f"{self.experiment_name or self.experiment} ({self.study}), replicate {self.replicate}"


def cultivation(exp: dict) -> str:
    """The experiment's cultivation mode, lowercased, or "unspecified" (which is not batch)."""
    return (exp.get("cultivationMode") or "unspecified").strip().lower()


def members(exp: dict) -> list:
    return [s.get("name", "") for s in exp.get("communityStrains", [])]


def medium_of(exp: dict) -> str:
    """The medium an experiment ran in, as mGrowthDB names it: one name per compartment, joined."""
    names = [(c.get("mediumName") or "").strip() for c in exp.get("compartments", [])]
    return "; ".join(dict.fromkeys(n for n in names if n))


def medium_key(name: str) -> str:
    """A medium name reduced to what tells media apart, so spellings of one medium group together.

    Case, punctuation and a parenthesized abbreviation are dropped: "Wilkins-Chalgren Anaerobe Broth (WC)"
    (study SMGDB00000007) and "Wilkins-Chalgren Anaerobe Broth" (study SMGDB00000009) are one medium. A
    second compartment medium is kept, so WC with mucin stays apart from WC.
    """
    import re
    parts = []
    for part in (name or "").split(";"):
        text = re.sub(r"\([^)]*\)", " ", part.casefold())
        text = " ".join(re.sub(r"[^0-9a-z]+", " ", text).split())
        if text:
            parts.append(text)
    return "; ".join(sorted(set(parts))) or "unnamed medium"


def _designation(name: str) -> str:
    return " ".join((name or "").split()[2:])


def strain_identities(exps, skipped) -> dict:
    """strain name -> node identity, as grownet: by NCBI taxon id (`ncbi:<id>`), or by genus and species when
    a study gives one id to different strains or none at all (grownet #23)."""
    names_by_taxon = {}
    for exp in exps:
        for strain in exp.get("communityStrains", []):
            if strain.get("NCBId") is not None and strain.get("name"):
                names_by_taxon.setdefault(str(strain["NCBId"]), set()).add(strain["name"])
    ambiguous = set()
    for taxon, names in names_by_taxon.items():
        if len({_designation(n) for n in names if _designation(n)}) > 1:
            ambiguous.add(taxon)
            skipped.append((f"taxon id {taxon}", "given to different strains in this study "
                            f"({', '.join(sorted(names))}); they are identified by genus and species instead"))
    identities = {}
    for exp in exps:
        for strain in exp.get("communityStrains", []):
            name = strain.get("name", "")
            taxon = None if strain.get("NCBId") is None else str(strain["NCBId"])
            if name in identities:
                continue
            if taxon is not None and taxon not in ambiguous:
                identities[name] = {"id": f"ncbi:{taxon}", "taxon_id": taxon, "species": genus_species(name),
                                    "identity": NCBI, "name": name}
            else:
                identities[name] = {"id": genus_species(name), "taxon_id": taxon or "",
                                    "species": genus_species(name), "identity": BY_NAME, "name": name}
    return identities


def _growth_rank(context: dict) -> tuple | None:
    """How good a measurement context is as the growth curve, smaller is better; None when it is none."""
    subject = context.get("subject") or {}
    technique = (context.get("techniqueType") or "").lower()
    if subject.get("type") not in ("strain", "bioreplicate") or technique in NOT_GROWTH:
        return None
    level = 0 if subject.get("type") == "strain" else 1
    kind = GROWTH_TECHNIQUES.index(technique) if technique in GROWTH_TECHNIQUES else len(GROWTH_TECHNIQUES)
    return (level, kind)


def _series(client, context_id, factor: float, scale: float = 1.0) -> list:
    """[(time in hours, value * scale)] in time order, without points that repeat a time."""
    points, seen = [], set()
    for t, v, _ in client.get_measurement_series(context_id):
        th = t * factor
        if th in seen:
            continue
        seen.add(th)
        points.append((th, v * scale))
    return sorted(points)


def cultures_of_experiment(client, exp: dict, identities: dict, excluded_metabolites=(),
                           with_metabolites: bool = True) -> tuple:
    """(cultures, skipped) for one batch monoculture experiment.

    With `with_metabolites` False, the cultures carry a growth curve only (the rate fallback)."""
    cultures, skipped = [], []
    label = exp.get("name") or exp.get("id") or "experiment"
    name = members(exp)[0]
    taxon = identities.get(name) or {"id": genus_species(name), "taxon_id": "", "species": genus_species(name),
                                     "identity": BY_NAME, "name": name}
    medium = medium_of(exp)
    excluded = {" ".join(m.casefold().split()) for m in excluded_metabolites}
    for stub in exp.get("bioreplicates", []):
        try:
            bio = client.get_bioreplicate(stub["id"])
        except Exception as e:  # noqa: BLE001 - one unreadable replicate must not lose the others
            skipped.append((f"{label}: {stub.get('name', stub.get('id'))}", f"{UNREAD} it: {e}"))
            continue
        rep = bio.get("name") or str(bio.get("id"))
        if bio.get("isAverage"):
            skipped.append((f"{label}: {rep}", AVERAGE))
            continue
        factor = hours(bio.get("measurementTimeUnits") or "")
        if factor is None:
            skipped.append((f"{label}: {rep}", f"time unit {bio.get('measurementTimeUnits')!r} is not one foodnet "
                            "converts to hours"))
            continue
        culture = Culture(taxon=taxon, study=str(exp.get("studyId", "")), experiment=str(exp.get("id", "")),
                          experiment_name=exp.get("name") or "", description=exp.get("description") or "",
                          medium=medium, medium_key=medium_key(medium), replicate=rep)
        growth = sorted(((r, c) for c in bio.get("measurementContexts", []) if (r := _growth_rank(c)) is not None),
                        key=lambda rc: rc[0])
        for _, context in growth:
            try:
                points = _series(client, context["id"], factor)
            except Exception as e:  # noqa: BLE001
                skipped.append((f"{label}: {rep}", f"{UNREAD} its growth curve: {e}"))
                continue
            if len(points) >= MIN_POINTS:
                culture.growth = {"times": [t for t, _ in points], "values": [v for _, v in points],
                                  "technique": context.get("techniqueType") or "",
                                  "level": (context.get("subject") or {}).get("type", ""),
                                  "unit": context.get("techniqueUnits") or context.get("techniqueType") or ""}
                break
        if with_metabolites:
            for context in bio.get("measurementContexts", []):
                subject = context.get("subject") or {}
                if subject.get("type") != "metabolite":
                    continue
                compound = compounds.canonical(subject.get("chebiId"), subject.get("name", ""))
                if " ".join((subject.get("name") or "").casefold().split()) in excluded or \
                        compound["name"].casefold() in excluded:
                    continue
                unit = context.get("techniqueUnits") or context.get("techniqueOriginalUnits") or ""
                scale, why = compounds.to_mm(unit, compound["chebi_id"])
                where = f"{label}: {rep}, {subject.get('name', '')}"
                if scale is None:
                    skipped.append((where, why))
                    continue
                try:
                    points = _series(client, context["id"], factor, scale)
                except Exception as e:  # noqa: BLE001
                    skipped.append((where, f"{UNREAD} its series: {e}"))
                    continue
                if len(points) < MIN_POINTS:
                    skipped.append((where, f"{len(points)} time point(s); a change needs {MIN_POINTS}"))
                    continue
                if compound["id"] in culture.metabolites:
                    skipped.append((where, "a second series of this compound in one replicate (recorded as "
                                    f"{subject.get('name')}); the first is used"))
                    continue
                culture.metabolites[compound["id"]] = {
                    "name": compound["name"], "chebi_id": compound["chebi_id"], "series": points,
                    "recorded_as": f"{subject.get('name', '')} (CHEBI:{subject.get('chebiId')})"
                    if subject.get("chebiId") else subject.get("name", ""),
                    "joined": compound["joined"]}
            if not culture.metabolites:
                continue
        elif culture.growth is None:
            continue
        cultures.append(culture)
    return cultures, skipped


def read_cultures(client, study_ids, keep=None, excluded_metabolites=(), include_non_batch: bool = False,
                  progress=None) -> dict:
    """{"cultures", "growth_only", "skipped", "experiments"} over several studies.

    `keep(name, taxon_id)` says whether a strain was asked for; None keeps every strain. `cultures` are the
    replicates with metabolite data; `growth_only` the replicates of the kept strains' other batch
    monocultures (no metabolites), for the growth-rate fallback. `experiments` maps every experiment id read
    to its record, for the report and the medium rule.
    """
    from .mgrowthdb import MGrowthDBError
    cultures, growth_only, skipped, experiments = [], [], [], {}
    for i, sid in enumerate(study_ids):
        if progress:
            progress(i, len(study_ids), f"Reading {sid} ({i + 1} of {len(study_ids)})")
        try:
            study = client.get_study(sid)
            exps = [client.get_experiment(e["id"]) for e in study.get("experiments", [])]
        except MGrowthDBError as e:
            skipped.append((sid, f"{UNREAD} it: {e}"))
            continue
        identities = strain_identities(exps, skipped)
        for exp in exps:
            names = members(exp)
            if len(names) != 1:
                continue
            ident = identities.get(names[0], {})
            if keep is not None and not keep(names[0], ident.get("taxon_id", "")):
                continue
            experiments[str(exp.get("id"))] = exp
            if cultivation(exp) != BATCH and not include_non_batch:
                skipped.append((exp.get("name") or exp.get("id"), f"{cultivation(exp)}, not batch: production and "
                                "consumption under dilution are not net changes in the vessel; left out"))
                continue
            found, skips = cultures_of_experiment(client, exp, identities, excluded_metabolites)
            skipped += skips
            if found:
                cultures += found
            else:
                rates, _ = cultures_of_experiment(client, exp, identities, with_metabolites=False)
                growth_only += rates
    return {"cultures": cultures, "growth_only": growth_only, "skipped": skipped, "experiments": experiments}
