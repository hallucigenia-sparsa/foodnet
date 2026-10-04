"""Fetch from mGrowthDB in parallel, so a search waits on the network once instead of hundreds of times.

grownet's prefetch, narrowed to what foodnet reads: the studies, their experiments, the bioreplicates of the
batch monocultures of the strains asked for, and every measurement series of those replicates (growth
curves and metabolites). The derivation then reads the same records from the client's in-memory cache, so
its result cannot differ from reading one request at a time. Nothing is kept between sessions: mGrowthDB
shows only the latest version of a study, and a stored copy could not tell it had gone stale (Karoline, for
grownet, 2026-09-27).
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

from .mgrowthdb import MGrowthDBError
from .reading import BATCH, cultivation, members

WORKERS = 6


def _each(fn, items, progress=None, message="", failures=None):
    """fn over items, WORKERS at a time; a failed item gives None (the derivation meets and reports the same
    failure when it asks for that item itself). `failures`, when given, collects (item, error)."""
    items = list(items)
    results = [None] * len(items)

    def run(i):
        try:
            results[i] = fn(items[i])
        except (MGrowthDBError, OSError, ValueError, KeyError) as e:
            results[i] = None
            if failures is not None:
                failures.append((items[i], e))
        return i

    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for done, _ in enumerate(pool.map(run, range(len(items))), start=1):
            if progress and (done == len(items) or done % 10 == 0):
                progress(done, len(items), f"{message} ({done} of {len(items)})")
    return results


def can_prefetch(client) -> bool:
    return all(hasattr(client, name) for name in
               ("get_study", "get_experiment", "get_bioreplicate", "get_measurement_series"))


def prefetch_studies(client, study_ids, keep=None, include_non_batch: bool = False, progress=None) -> None:
    """Load what reading these studies' monocultures needs, a few requests at a time."""
    if not can_prefetch(client):
        return
    studies = [s for s in _each(client.get_study, study_ids) if s]
    experiment_ids = [e["id"] for s in studies for e in s.get("experiments", [])]
    experiments = [e for e in _each(client.get_experiment, experiment_ids, progress, "Reading experiments") if e]
    wanted = []
    for e in experiments:
        names = members(e)
        if len(names) != 1 or (cultivation(e) != BATCH and not include_non_batch):
            continue
        ids = {str(s.get("NCBId")) for s in e.get("communityStrains", []) if s.get("NCBId") is not None}
        if keep is not None and not keep(names[0], next(iter(ids), "")):
            continue
        wanted.append(e)
    stubs = [b["id"] for e in wanted for b in e.get("bioreplicates", [])]
    bioreplicates = [b for b in _each(client.get_bioreplicate, stubs, progress, "Reading replicates") if b]
    contexts = [c["id"] for b in bioreplicates if not b.get("isAverage")
                for c in b.get("measurementContexts", []) if (c.get("techniqueType") or "") != "ph"]
    _each(client.get_measurement_series, contexts, progress, "Reading growth curves and metabolites")


def study_ids_in_order(client, study_id_format: str, max_studies: int, miss_run: int) -> list:
    """The ids of the studies mGrowthDB holds, in id order, stopping after `miss_run` absent ids in a row
    (grownet's crawl). Only "no such study" (HTTP 404) ends it; anything else is raised."""
    found, misses, n = [], 0, 1
    while n <= max_studies and misses < miss_run:
        batch = list(range(n, min(n + WORKERS, max_studies + 1)))
        failures = []
        studies = _each(client.get_study, [study_id_format.format(i) for i in batch], failures=failures)
        down = [(sid, e) for sid, e in failures if getattr(e, "status", None) != 404]
        if down:
            sid, e = down[0]
            raise MGrowthDBError(f"mGrowthDB could not be read (asking for {sid}): {e}",
                                 status=getattr(e, "status", None))
        for i, study in zip(batch, studies, strict=True):
            if misses >= miss_run:
                break
            if study is None:
                misses += 1
            else:
                misses = 0
                found.append(study_id_format.format(i))
        n = batch[-1] + 1
    return found
