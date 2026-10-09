"""What foodnet read of a study, kept on this machine while the study is unchanged (Karoline, 2026-10-09).

mGrowthDB shows only the latest version of a study, so foodnet used to read everything fresh. A study's
`uploadedAt` changes whenever it is submitted again (Karoline: "I know that uploadedAt is reliable; it's tied
to the submission process"), so what was read of it can be used again while its record, read fresh, gives the
same `uploadedAt`: its experiments, its replicates and their series. Since mGrowthDB also changes what it
serves without a submission (it converts units with shared metabolite masses, and may recreate numeric ids),
a kept study is read again after MAX_AGE as well (Karoline, 2026-10-09: "30 days"). The study record itself is
never kept, nor a study without an `uploadedAt`. `refresh` on the client (`foodnet derive --refresh`, the page's
option) uses nothing kept and keeps what it reads.

The responses are kept as mGrowthDB gave them, per study, in the cache folder (`foodnet.taxonomy.cache_folder`),
never in the repository, with when they were first read; the data versions of a result say which studies came
from here and since when (`foodnet.mgrowthdb.data_versions`).
"""
from __future__ import annotations

import json
import os
import re
import tempfile
import threading
import time

from .taxonomy import cache_folder

FORMAT = "foodnet.study_store/v1"
MAX_AGE = 30 * 24 * 3600
# what a study id may look like before it becomes a file name (a review: an id from a mirror could hold "../")
_STUDY_ID = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
_lock = threading.Lock()


def _usable(client) -> bool:
    return getattr(client, "base_url", None) is not None and hasattr(client, "_cache_get") \
        and hasattr(client, "_cache_put")


def _path(study_id: str) -> str:
    return os.path.join(cache_folder(), "studies", f"{study_id}.json")


def _read(study_id: str) -> dict | None:
    """The kept file of a study, or None when there is none or it is not one foodnet wrote whole."""
    try:
        with open(_path(study_id), encoding="utf-8") as f:
            kept = json.load(f)
        if not isinstance(kept, dict) or kept.get("format") != FORMAT or \
                not isinstance(kept.get("responses"), dict) or not isinstance(kept.get("read_at"), (int, float)):
            return None
        return kept
    except (OSError, ValueError):
        return None


def _matches(kept, client, study: dict) -> bool:
    return bool(kept) and kept.get("source") == client.base_url and kept.get("uploadedAt") == study["uploadedAt"] \
        and 0 <= time.time() - kept["read_at"] < MAX_AGE


def write_whole(path: str, data) -> bool:
    """Write JSON to `path` through a temporary file of its own in the same folder, so another process never
    reads half a file nor replaces it with one (a review: one shared ".part" name was truncated under a reader)."""
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        handle, part = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".part")
        try:
            with os.fdopen(handle, "w", encoding="utf-8") as f:
                json.dump(data, f)
            os.replace(part, path)
        finally:
            if os.path.exists(part):
                os.remove(part)
        return True
    except OSError:
        return False                           # not kept: read again next time


def load(client, study: dict, study_id: str | None = None) -> int:
    """Put what is kept of a study (`study` is its fresh record) into the client's memory, when it was kept from
    the same mGrowthDB under the same `uploadedAt` less than MAX_AGE ago and the client is not refreshing; the
    number of responses it gave. Anything unexpected gives 0: a search never fails for the store."""
    try:
        sid = study_id or (study or {}).get("id", "")
        if not _usable(client) or getattr(client, "refresh", False) or not study or \
                not study.get("uploadedAt") or not _STUDY_ID.match(str(sid)):
            return 0
        kept = _read(sid)
        if not _matches(kept, client, study):
            return 0
        loaded = getattr(client, "loaded_urls", None)
        # a copy written after this client started (a refresh or re-read elsewhere) is newer than anything it
        # holds of the study; otherwise only what came from an older copy is replaced (reviews: a page's memory
        # outlived a refresh, in another process too)
        newer = kept["read_at"] > getattr(client, "started_at", float("inf")) and \
            getattr(client, "written", {}).get(sid) != kept["read_at"]
        series = getattr(client, "_series", None)
        if newer and hasattr(client, "_mem") and isinstance(client._mem, dict):
            # what this client holds of the study that the newer copy lacks is forgotten and read again, so one
            # study is never part old and part new (a review: a refresh of some replicates)
            for url in _urls(client, study):
                if url not in kept["responses"] and url in client._mem:
                    del client._mem[url]
                    if series is not None and "/bioreplicate/" in url and url.endswith(".csv"):
                        series.pop(url.rsplit("/", 1)[-1][:-4], None)
        for url, data in kept["responses"].items():
            if client._cache_get(url) is None or newer or (loaded is not None and url in loaded):
                if client._cache_get(url) != data and series is not None and "/bioreplicate/" in url \
                        and url.endswith(".csv"):
                    series.pop(url.rsplit("/", 1)[-1][:-4], None)    # its parsed series go with it
                client._cache_put(url, data)
                if loaded is not None:
                    loaded.add(url)
        if hasattr(client, "kept_studies") and sid not in getattr(client, "saved_studies", ()):
            client.kept_studies[sid] = kept["read_at"]
        return len(kept["responses"])
    except Exception:  # noqa: BLE001 - the store is a convenience; the search reads mGrowthDB instead
        return 0


def _urls(client, study: dict) -> list:
    """The responses of this study in the client's memory: its experiments, their replicates, and the
    replicates' series (one CSV each, or one per series where that was read)."""
    base = client.base_url
    urls = []
    for e in study.get("experiments", []):
        url = f"{base}/experiment/{e['id']}.json"
        record = client._cache_get(url)
        if not isinstance(record, dict):
            continue
        urls.append(url)
        for b in record.get("bioreplicates", []):
            for kind in ("json", "csv"):
                urls.append(f"{base}/bioreplicate/{b['id']}.{kind}")
            bio = client._cache_get(f"{base}/bioreplicate/{b['id']}.json")
            if isinstance(bio, dict):
                urls += [f"{base}/measurement-context/{c['id']}.csv" for c in bio.get("measurementContexts", [])]
    return urls


def save(client, study: dict, study_id: str | None = None) -> int:
    """Keep what the client read of a study under its `uploadedAt`, added to what was kept under the same one
    (which keeps its first read time); the number of responses kept. Nothing is written when nothing new was
    read. Anything unexpected keeps nothing."""
    try:
        sid = study_id or (study or {}).get("id", "")
        if not _usable(client) or not study or not study.get("uploadedAt") or not _STUDY_ID.match(str(sid)):
            return 0
        found = {}
        for url in _urls(client, study):
            data = client._cache_get(url)
            if data is not None:
                found[url] = data
        if not found:
            return 0
        with _lock:
            kept = _read(sid)
            loaded = getattr(client, "loaded_urls", set())
            # a copy another process wrote after this client started (a refresh, a re-read after 30 days, a newer
            # uploadedAt) is left alone: what this client read may be older (reviews: a page's hour-old reads were
            # written back over a command-line refresh)
            written = getattr(client, "written", {})
            if kept and kept.get("source") == client.base_url and not getattr(client, "refresh", False) and \
                    kept["read_at"] > getattr(client, "started_at", float("inf")) and \
                    written.get(sid) != kept["read_at"]:
                return 0
            if _matches(kept, client, study) and not getattr(client, "refresh", False):
                # only what this client read from mGrowthDB is added: what it took from a copy never overwrites
                # the kept one (a review: a page's memory wrote the pre-refresh values back)
                live = {url: data for url, data in found.items() if url not in loaded}
                if all(kept["responses"].get(url) == data for url, data in live.items()):
                    return len(kept["responses"])       # nothing new: no rewrite
                responses, read_at = {**kept["responses"], **live}, kept["read_at"]
            else:
                # a new copy: only what this client read from mGrowthDB, never what it loaded from an older copy
                # (a review: an expired copy still in the page's memory was saved again as newly read); a refresh
                # replaces the copy whole (a review: it read everything again but kept the old file)
                responses = {url: data for url, data in found.items() if url not in loaded}
                if not responses:
                    return 0
                read_at = time.time()
            record = {"format": FORMAT, "source": client.base_url, "study": sid, "uploadedAt": study["uploadedAt"],
                      "read_at": read_at, "responses": responses}
            if not write_whole(_path(sid), record):
                return 0
            if hasattr(client, "saved_studies"):
                client.saved_studies.add(sid)
            if hasattr(client, "written"):
                client.written[sid] = read_at            # this client may add to the copy it wrote, while it stands
            return len(responses)
    except Exception:  # noqa: BLE001 - the store is a convenience; nothing kept is the safe outcome
        return 0
