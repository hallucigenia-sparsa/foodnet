"""Resolve species names to NCBI taxon ids, from mGrowthDB's own records.

mGrowthDB queries take NCBI taxon ids (`MGrowthDBClient.search`), while people type species names. Every
strain in mGrowthDB carries its taxon id (`NCBId` on an experiment's community strains), so the mapping is
built from the database itself: no second external service, and only names that mGrowthDB actually holds,
which are the only ones with data to find.

`species_index` crawls the studies once (the client caches responses, on disk when a cache_dir is set) and
returns genus and species keys mapped to the taxon ids seen under them. `resolve_species` turns what a
person typed, names or numeric taxon ids, into ids, and reports what mGrowthDB does not hold.

The species list is kept on this machine for a day (Karoline, 2026-10-09: "you can keep the species list for
a day", after a search spent three of its minutes crawling a slow mGrowthDB): names, taxon ids and the studies
holding them, in the user's cache folder, never in the repository, and never a list the crawl could not read
whole. What was read of a study is kept beside it while the study is unchanged (foodnet.store;
docs/DATA_GOVERNANCE.md).
"""
from __future__ import annotations

import difflib
import json
import os
import re
import sys
import time

from .model import genus_name, genus_species

STUDY_ID = "SMGDB{:08d}"
# stop crawling after this many consecutive study ids are absent. Study 3 is already missing; a larger gap
# (withdrawn or unpublished studies) would have hidden every later study from name lookups at 5, and now that
# the crawl runs in parallel, looking 25 ids further costs well under a second
MISS_RUN = 25
# past the last study a previous crawl found, this many absent ids in a row end it; before it, absent ids end
# nothing, since a later study is known to exist. It was 10 (to query mGrowthDB less); a review found that a
# study after a longer gap would then be missed for good, on machines with a kept list only, for 15 requests a
# day, so it is MISS_RUN again
MISS_AFTER_KNOWN = MISS_RUN
MAX_STUDIES = 500     # a hard stop, so a crawl can never run away


def _strain_entries(exp: dict):
    """(name, taxon id) for every community strain of an experiment that carries a taxon id."""
    for strain in exp.get("communityStrains", []):
        name, taxon = strain.get("name"), strain.get("NCBId")
        if name and taxon is not None:
            yield name, int(taxon)


class SpeciesIndex(dict):
    """The species list: genus and species key -> {taxon id: a name seen for it}, as a plain dict, plus
    `current`: taxon id -> its current name, the one used by the most recently published study holding it
    (Karoline, grownet #24, 2026-09-18). Old names stay in the list, so they still resolve. `studies`: the ids of
    every study the crawl found, in id order, which the page's All button derives."""

    def __init__(self, *args, current=None, studies=(), where=None, failed=(), **kwargs):
        super().__init__(*args, **kwargs)
        # what the crawl could not read: a study or experiment missing here hides its taxa from every search,
        # so a search reports it as an error and the page does not keep an index that holds any
        self.failed = list(failed)
        self.current = dict(current or {})
        self.studies = list(studies)
        # taxon id -> the studies with an experiment holding it. foodnet finds studies here, not through
        # mGrowthDB's search, which matches per-strain measurements only and so misses a monoculture measured
        # at the culture level (study SMGDB00000009, found 2026-10-04)
        self.where = {t: list(s) for t, s in (where or {}).items()}
        # when the crawl ran (seconds since the epoch), so a kept list can be told its age, and whether this one
        # came from the cache folder
        self.built_at = time.time()
        self.kept = False

    def studies_of(self, taxon_ids) -> list:
        """The studies holding any of these taxon ids, in id order."""
        found = {sid for t in taxon_ids for sid in self.where.get(int(t), ())}
        return [sid for sid in self.studies if sid in found]


def species_index(client, max_studies: int = MAX_STUDIES, progress=None, known_last: int = 0) -> dict:
    """Map a genus and species key to {taxon id: a name seen for it}, crawled from mGrowthDB.

    Study ids are consecutive, so the crawl walks them and stops after MISS_RUN absent ids in a row. The
    studies and their experiments are read in parallel (`foodnet.fetch`), then walked in id order, so the
    first name seen for a taxon is the same as a one-by-one crawl would keep. A study or experiment that
    cannot be read is skipped, since a partial index is more useful than no index, and recorded in `failed`,
    since its taxa are then missing from every search.
    """
    from .fetch import _each, study_ids_in_order

    if not (hasattr(client, "get_study") and hasattr(client, "get_experiment")):
        return _species_index_one_by_one(client, max_studies)
    from . import store
    ids = study_ids_in_order(client, STUDY_ID, max_studies, MISS_RUN, known_last)
    failed = []
    experiment_ids = []
    read = {}
    for sid in ids:
        try:
            study = client.get_study(sid) or {}
            store.load(client, study, sid)             # the experiments of an unchanged study, kept before
            experiment_ids += [e["id"] for e in study.get("experiments", [])]
            read[sid] = study
        except Exception as e:  # noqa: BLE001 - recorded, and the study is skipped below
            failed.append((sid, str(e)))
    failures = []
    _each(client.get_experiment, experiment_ids, progress, "Reading the species list of mGrowthDB",
          failures=failures)
    for sid, study in read.items():                       # a study that could not be read is not asked again
        store.save(client, study, sid)
    failed += [(eid, str(e)) for eid, e in failures]
    index, published, where = {}, [], {}
    for order, sid in enumerate(ids):
        try:
            experiments = client.study_experiments(sid)       # from the cache the parallel reads filled
        except Exception as e:  # noqa: BLE001 - one unreadable study is skipped, the rest still count
            if not any(sid == x for x, _ in failed):
                failed.append((sid, str(e)))
            continue
        for exp in experiments:
            for name, taxon in _strain_entries(exp):
                index.setdefault(genus_species(name), {}).setdefault(taxon, name)
                where.setdefault(taxon, set()).add(sid)
        published.append(((client.get_study(sid) or {}).get("publishedAt") or "", order, experiments))
    current = {}
    for _, _, experiments in sorted(published, key=lambda p: (p[0], p[1])):    # the latest publication last
        for exp in experiments:
            for name, taxon in _strain_entries(exp):
                current[taxon] = name
    return SpeciesIndex(index, current=current, studies=ids, where=where, failed=failed)


# how long a species list is kept (seconds), and the format of the kept file
KEEP_FOR = 24 * 3600
KEPT_FORMAT = "foodnet.species_list/v1"


def cache_folder() -> str:
    """Where foodnet keeps its species list and what it read of each study: FOODNET_CACHE_DIR, or the system's
    cache folder for programs."""
    if os.environ.get("FOODNET_CACHE_DIR"):
        return os.environ["FOODNET_CACHE_DIR"]
    if sys.platform == "darwin":
        return os.path.expanduser("~/Library/Caches/foodnet")
    if sys.platform.startswith("win"):
        return os.path.join(os.environ.get("LOCALAPPDATA") or os.path.expanduser("~"), "foodnet", "Cache")
    return os.path.join(os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache"), "foodnet")


def _kept_path() -> str:
    return os.path.join(cache_folder(), "species_list.json")


def _load_kept(source: str, keep_for: float):
    """The kept species list for this mGrowthDB, or None when there is none, it is older than `keep_for`, or
    it cannot be read."""
    try:
        with open(_kept_path(), encoding="utf-8") as f:
            kept = json.load(f)
        if kept.get("format") != KEPT_FORMAT or kept.get("source") != source:
            return None
        if not 0 <= time.time() - float(kept["built_at"]) < keep_for:
            return None
        index = SpeciesIndex({key: {int(t): n for t, n in ids.items()} for key, ids in kept["index"].items()},
                             current={int(t): n for t, n in kept["current"].items()}, studies=kept["studies"],
                             where={int(t): s for t, s in kept["where"].items()})
        index.built_at = float(kept["built_at"])
        index.kept = True
        return index
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None


def _keep(index, source: str) -> None:
    """Write the species list to the cache folder; a list the crawl could not read whole is never kept."""
    if index.failed:
        return
    kept = {"format": KEPT_FORMAT, "source": source, "built_at": index.built_at,
            "index": {key: {str(t): n for t, n in ids.items()} for key, ids in index.items()},
            "current": {str(t): n for t, n in index.current.items()}, "studies": index.studies,
            "where": {str(t): s for t, s in index.where.items()}}
    from .store import write_whole
    write_whole(_kept_path(), kept)           # a list that cannot be kept is crawled again next time


def kept_species_index(client, progress=None, refresh: bool = False, keep_for: float = KEEP_FOR):
    """The species list, from the cache folder when one younger than `keep_for` is there for this mGrowthDB,
    else crawled (`species_index`) and kept. `refresh` crawls whatever is kept. A client without a base url
    (a test double) is always crawled and never kept."""
    source = getattr(client, "base_url", None)
    if refresh and hasattr(client, "refresh"):
        client.refresh = True             # before the crawl, which would otherwise load the kept studies (a review)
    if source is None:
        return species_index(client, progress=progress)
    if not refresh:
        kept = _load_kept(source, keep_for)
        if kept is not None:
            return kept
    # the last study the previous list found, however old: past it, fewer absent ids end the crawl
    before = _load_kept(source, float("inf"))
    known_last = max((int(s[5:]) for s in (before.studies if before else []) if s[5:].isdigit()), default=0)
    index = species_index(client, progress=progress, known_last=known_last)
    _keep(index, source)
    return index


def _species_index_one_by_one(client, max_studies: int) -> dict:
    """The crawl for a client that only lists a study's experiments (as test doubles do)."""
    index, misses, found, where = {}, 0, [], {}
    for n in range(1, max_studies + 1):
        if misses >= MISS_RUN:
            break
        try:
            experiments = client.study_experiments(STUDY_ID.format(n))
        except Exception:      # noqa: BLE001 - an absent id is expected; any other failure skips one study
            misses += 1
            continue
        misses = 0
        found.append(STUDY_ID.format(n))
        for exp in experiments:
            for name, taxon in _strain_entries(exp):
                index.setdefault(genus_species(name), {}).setdefault(taxon, name)
                where.setdefault(taxon, set()).add(STUDY_ID.format(n))
    return SpeciesIndex(index, studies=found, where=where)


# a list marker or numbering someone pasted along with a name: "- ", "* ", "1. ", "2) "
_BULLET = re.compile(r"^\s*(?:[-*\u2022]+|\d+[.)])\s+")
# an NCBI taxon id written the ways NCBI and papers write it: 476272, txid476272, NCBI:txid476272, taxid:476272
TAXON_ID = re.compile(r"(?i)^(?:ncbi\s*[:_-]?\s*)?(?:txid|taxid|taxon(?:\s*id)?)?\s*[:=]?\s*(\d+)$")
# what a species or strain name can hold: it starts with a letter (or "[" as in "[Clostridium] scindens")
_NAME = re.compile(r"^[\[A-Za-z][A-Za-z0-9 .,'\-\[\]()/:+]*$")


def split_entries(lines) -> list:
    """The entries in what someone typed or pasted: one per line, and also split at commas and semicolons,
    which no species or strain name contains; list markers are dropped."""
    out = []
    for line in lines:
        for part in re.split(r"[,;]", line or ""):
            part = _BULLET.sub("", part).strip()
            if part:
                out.append(part)
    return out


def _display(key: str) -> str:
    return key[:1].upper() + key[1:]


def _suggestions(text: str, index: dict) -> list:
    """Names mGrowthDB holds that the person may have meant: an abbreviated genus ("B. hydrogenotrophica"),
    or a close spelling."""
    key = genus_species(text)
    words = key.split()
    if len(words) >= 2 and (words[0].endswith(".") or len(words[0]) == 1):
        initial, epithet = words[0].rstrip(".")[:1], words[1]
        return [_display(k) for k in sorted(index)
                if len(k.split()) > 1 and k.split()[0].startswith(initial) and k.split()[1] == epithet][:8]
    return [_display(k) for k in difflib.get_close_matches(key, list(index), n=3, cutoff=0.8)]


def _genus_suggestions(key: str, index: dict) -> list:
    """Genera mGrowthDB holds that are spelled like the one typed."""
    genera = sorted({k.split()[0] for k in index})
    return [_display(g) for g in difflib.get_close_matches(key, genera, n=3, cutoff=0.8)]


def resolve_species(entries, index: dict) -> dict:
    """Resolve typed species names and taxon ids against an index from `species_index`.

    entries: strings, each a species name ("Faecalibacterium prausnitzii", strain designations and case
    ignored), a genus ("Blautia": every species of it that mGrowthDB holds, listed in "genera"), or a
    numeric NCBI taxon id ("853"), which is taken as given.

    One species name often carries several taxon ids in mGrowthDB: a species-level id and strain-level ids
    below it. They are the same species, so a name resolves to all of them and the caller shows which
    strains were used rather than asking the person to choose.

    An entry may carry a list marker ("- ", "1. ") and an id may be written as NCBI writes it ("txid476272").
    Returns {"taxon_ids": [ids in the order first seen], "resolved": [(entry, {taxon id: name})],
    "unresolved": [entries that gave nothing], "reasons": {entry: why}, "suggestions": {entry: [names]}}:
    unreadable text, a taxon id no strain in mGrowthDB carries, or a genus or name mGrowthDB does not hold,
    with the names it does hold that the person may have meant; "genera": {entry: [species]} for each
    genus entered.
    """
    out = {"taxon_ids": [], "resolved": [], "unresolved": [], "reasons": {}, "suggestions": {}, "genera": {}}

    def fail(text, reason, suggestions=()):
        out["unresolved"].append(text)
        out["reasons"][text] = reason
        if suggestions:
            out["suggestions"][text] = list(suggestions)

    for entry in entries:
        text = _BULLET.sub("", (entry or "")).strip()
        if not text:
            continue
        taxon_id = TAXON_ID.match(text)
        if taxon_id:
            taxon = int(taxon_id.group(1))
            name = next((names[taxon] for names in index.values() if taxon in names), None)
            if name is None:
                # an id no strain in mGrowthDB carries can give no result, so it is not taken as found
                fail(text, f"no strain in mGrowthDB has NCBI taxon id {taxon}")
                continue
            matches = {taxon: name}
        elif not _NAME.match(text):
            fail(text, "not readable as a species or strain name, nor as an NCBI taxon id")
            continue
        else:
            key = genus_species(text)
            matches = index.get(key)
            if not matches and len(key.split()) == 1:
                # a genus alone stands for every strain of it in mGrowthDB (Karoline, 2026-09-28), by the rule
                # the genus merge uses: "unclassified Bacteroides" is Bacteroides, "[Clostridium]" is not
                # Clostridium
                species = sorted(k for k in index if genus_name(k).lower() == genus_name(key).lower())
                if species:
                    matches = {t: n for k in species for t, n in index[k].items()}
                    out["genera"][text] = [_display(k) for k in species]
            if not matches:
                genus = len(key.split()) == 1
                fail(text, "no genus or species of this name in mGrowthDB" if genus
                     else "no species or strain of this name in mGrowthDB",
                     _genus_suggestions(key, index) if genus else _suggestions(text, index))
                continue
            matches = dict(matches)
        out["resolved"].append((text, matches))
        for taxon in matches:
            if taxon not in out["taxon_ids"]:
                out["taxon_ids"].append(taxon)
    return out
