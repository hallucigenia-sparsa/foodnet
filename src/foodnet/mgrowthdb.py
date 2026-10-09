"""mGrowthDB API client.

Wired against the live REST API, base `https://mgrowthdb.gbiomed.kuleuven.be/api/v1/`, documented at
https://mgrowthdb.readthedocs.io/en/latest/api.html .

The API serves raw data (study -> experiments -> bioreplicates -> measurement contexts). A metabolite is
a measurement context whose subject has type "metabolite" (with its name and ChEBI id), recorded per
bioreplicate beside the growth curves; foodnet derives production and consumption from those series (see
foodnet.derive). Public data needs no auth. mGrowthDB is open (see docs/DATA_GOVERNANCE.md); foodnet pulls
from it but never commits raw or pulled data.

The client is grownet's, with foodnet's additions: a fallback past an address that does not answer
(open_socket), slowness counting, and a replicate's series read from its one CSV.

The client works for ANY study id (get_study, get_experiment, study_experiments). It caches responses
in memory for the life of the client (and optionally on disk via `cache_dir`) so repeated pulls of the
same study do not re-hit the API, and it retries transient network failures and 5xx responses with a
short backoff. Nothing pulled is committed to the repository. What was read of a study is kept on this machine
while the study is unchanged (foodnet.store), and the species list for a day (foodnet.taxonomy).
"""
from __future__ import annotations

import csv
import datetime
import hashlib
import http.client
import io
import json
import os
import socket
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

from . import __version__
from .brand import NAME

MGROWTHDB_API = "https://mgrowthdb.gbiomed.kuleuven.be/api/v1"
API_DOCS = "https://mgrowthdb.readthedocs.io/en/latest/api.html"


class MGrowthDBError(RuntimeError):
    """A failed mGrowthDB request: a bad id, a network error, or the API being unreachable.

    `status` is the HTTP status when mGrowthDB answered (404 for an id it does not hold), and None when it
    could not be reached, so a caller can tell "no such study" from "mGrowthDB is down" (step 7 of the
    audit, 2026-09-28: an unreachable mGrowthDB read as an empty database)."""

    def __init__(self, message: str, status: int | None = None):
        super().__init__(message)
        self.status = status


class _Status(Exception):
    """An HTTP error status from mGrowthDB, raised by `_send` for the retry logic to judge."""

    def __init__(self, code: int, reason: str = ""):
        super().__init__(f"{code} {reason}".strip())
        self.code = code


# what a dropped or refused connection raises below urllib: retried like a network error
_NETWORK = (OSError, http.client.HTTPException)


# a request that takes longer than SLOW seconds counts as slow, and the progress text says so (Karoline,
# 2026-10-09, after a search took over 8 minutes)
SLOW = 10.0
# how long one address may take to accept a connection before the next is tried. On 2026-10-09 mGrowthDB's IPv6
# address did not accept connections over the KU Leuven VPN, which routes KU Leuven's IPv6 prefix into its
# tunnel, while its IPv4 address answered at once; every new connection waited 17 s for the operating system to
# give up on IPv6. This is a sequential fallback, simpler than the overlapping attempts of RFC 8305 (Happy
# Eyeballs): each address gets CONNECT_TIMEOUT, and each that timed out is tried once more with the full
# timeout, so a slow but working path still connects. The family that answered is tried first afterwards.
CONNECT_TIMEOUT = 3.0


def open_socket(host: str, port: int, timeout, preferred: list):
    """A connected socket to `host`: its addresses in turn with CONNECT_TIMEOUT each, the family in `preferred`
    first, then each one that timed out (not refused) once more with the full `timeout` (None: no limit), so a
    slow but working address connects wherever it is in the list (a review); `preferred` is updated to the
    family that answered. Every attempt's error is in the one raised when none answered."""
    infos = socket.getaddrinfo(host, port, type=socket.SOCK_STREAM)
    if preferred:
        infos.sort(key=lambda info: info[0] != preferred[0])
    quick = CONNECT_TIMEOUT if timeout is None else min(CONNECT_TIMEOUT, timeout)
    errors = []

    def attempt(info, limit):
        family, kind, proto, _, address = info
        sock = None
        try:
            # created inside the try, as socket.create_connection does: a machine without IPv6 refuses the
            # socket itself, and the next address is then tried (a review)
            sock = socket.socket(family, kind, proto)
            sock.settimeout(limit)
            sock.connect(address)
        except OSError as e:
            errors.append(f"{address[0]}: {e}")
            if sock is not None:
                sock.close()
            return None, isinstance(e, TimeoutError)
        sock.settimeout(timeout)
        preferred[:] = [family]
        return sock, False

    slow = []
    for info in infos:
        sock, timed_out = attempt(info, quick)
        if sock is not None:
            return sock
        if timed_out:
            slow.append(info)
    if timeout is None or timeout > quick:
        for info in slow:
            sock, _ = attempt(info, timeout)
            if sock is not None:
                return sock
    raise OSError(f"could not connect to {host} ({'; '.join(errors) or 'no address'})")


class MGrowthDBClient:
    """A thin read-only client over the mGrowthDB REST API (public endpoints, no auth needed).

    Args:
      base_url: API base; override to point at a mirror or a test server.
      timeout:  per-request timeout in seconds.
      retries:  attempts on a transient failure (network error or HTTP 5xx) before giving up.
      backoff:  base seconds between retries (grows linearly with the attempt).
      cache:    keep an in-memory response cache for the life of the client.
      cache_dir: optional directory for an on-disk JSON cache across runs (never inside the repo).
    """

    def __init__(self, base_url: str = MGROWTHDB_API, timeout: int = 30, retries: int = 3,
                 backoff: float = 0.5, cache: bool = True, cache_dir: str | None = None):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.retries = max(1, int(retries))
        self.backoff = max(0.0, float(backoff))
        self._mem = {} if cache else None
        self.cache_dir = cache_dir
        # one kept-open connection per thread: a new HTTPS connection per request cost about 120 ms against
        # 55 ms on an open one, and the parallel prefetch (foodnet.fetch) gives each worker its own
        self._local = threading.local()
        self._chebi = None
        # what slowness was seen, and the address family that answered
        self._stats = {"slow": 0, "retried": 0}
        self._family = []
        self._stats_lock = threading.Lock()
        # measurement context id -> its bioreplicate, and a bioreplicate's series read from its one CSV
        self._owner = {}
        self._series = {}
        # read nothing kept on disk (foodnet.store), and which studies came from there, since when
        self.refresh = False
        self.kept_studies = {}
        self.saved_studies = set()        # studies this client read from mGrowthDB and kept: not "from the copy"
        self.loaded_urls = set()          # responses this client took from the copy, not from mGrowthDB
        self.written = {}                 # study id -> read_at of the copy this client wrote
        self.started_at = time.time()     # what it read is reused while it lives (the page keeps one an hour)
        self.searches = 0                 # how many searches it served
        parsed = urllib.parse.urlsplit(self.base_url)
        self._scheme, self._host, self._prefix = parsed.scheme, parsed.netloc, parsed.path
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _disk_path(self, key: str) -> str:
        h = hashlib.sha256(key.encode("utf-8")).hexdigest()[:32]
        return os.path.join(self.cache_dir, f"{h}.json")

    def _cache_get(self, key: str):
        if self._mem is not None and key in self._mem:
            return self._mem[key]
        if self.cache_dir:
            p = self._disk_path(key)
            if os.path.exists(p):
                try:
                    with open(p, encoding="utf-8") as f:
                        data = json.load(f)
                    if self._mem is not None:
                        self._mem[key] = data
                    return data
                except (OSError, ValueError):
                    return None
        return None

    def _cache_put(self, key: str, data) -> None:
        if self._mem is not None:
            self._mem[key] = data
        if self.cache_dir:
            try:
                with open(self._disk_path(key), "w", encoding="utf-8") as f:
                    json.dump(data, f)
            except OSError:
                pass

    def _connection(self, fresh: bool = False):
        conn = getattr(self._local, "conn", None)
        if fresh or conn is None:
            if conn is not None:
                conn.close()
            kind = http.client.HTTPSConnection if self._scheme == "https" else http.client.HTTPConnection
            conn = self._local.conn = kind(self._host, timeout=self.timeout)
            conn._create_connection = lambda address, timeout, source=None: open_socket(
                address[0], address[1], timeout, self._family)
        return conn

    def _send(self, url: str, accept: str) -> bytes:
        """One request over this thread's open connection, reopened once if the server dropped it."""
        path = url[len(f"{self._scheme}://{self._host}"):] if url.startswith(f"{self._scheme}://") else url
        headers = {"Accept": accept, "User-Agent": f"foodnet/{__version__}", "Connection": "keep-alive"}
        for attempt in (0, 1):
            conn = self._connection(fresh=attempt == 1)
            was_open = conn.sock is not None
            try:
                conn.request("GET", path, headers=headers)
                response = conn.getresponse()
                body = response.read()
            except _NETWORK:
                # only a kept-open connection the server had closed is reopened; a connection that could not be
                # made at all would only try every address again (a review: 5 minutes to report an outage)
                if attempt == 1 or not was_open:
                    raise
                continue
            if response.status >= 400:
                raise _Status(response.status, response.reason)
            return body
        raise OSError("unreachable")               # not reached

    def _timed_send(self, url: str, accept: str) -> bytes:
        """`_send`, counting slow answers for the progress text."""
        start = time.monotonic()
        try:
            return self._send(url, accept)
        finally:
            if time.monotonic() - start > SLOW:
                with self._stats_lock:
                    self._stats["slow"] += 1

    def slowness(self) -> str:
        """What the progress text adds while mGrowthDB answers slowly, or "" (Karoline, 2026-10-09)."""
        with self._stats_lock:
            slow, retried = self._stats["slow"], self._stats["retried"]
        if not slow and not retried:
            return ""
        parts = ([f"{slow} request(s) took over {SLOW:g} s"] if slow else []) + \
            ([f"{retried} retried"] if retried else [])
        return f"mGrowthDB is answering slowly: {', '.join(parts)}"

    def _request(self, url: str, accept: str) -> bytes:
        """`_send` with the retry rule: a 4xx is a real error (a bad id) and is not retried; a 5xx or a
        network failure may be transient and is retried with a growing pause. JSON and CSV alike."""
        last = None
        for attempt in range(self.retries):
            if attempt:
                with self._stats_lock:
                    self._stats["retried"] += 1
            try:
                return self._timed_send(url, accept)
            except _Status as e:
                if e.code < 500:
                    raise MGrowthDBError(
                        f"mGrowthDB returned HTTP {e.code} for {url} (check the id; API docs: {API_DOCS})",
                        status=e.code) from None
                last = MGrowthDBError(f"mGrowthDB returned HTTP {e.code} for {url} (server error)", status=e.code)
            except _NETWORK as e:
                last = MGrowthDBError(
                    f"could not reach mGrowthDB at {url}: {getattr(e, 'reason', e) or type(e).__name__} "
                    "(check your network; the API may be temporarily down)"
                )
            if attempt < self.retries - 1:
                time.sleep(self.backoff * (attempt + 1))
        raise last

    def _get_text(self, path: str) -> str:
        """Fetch a non-JSON representation (the CSV of a measurement context), cached like the rest."""
        url = f"{self.base_url}/{path.lstrip('/')}"
        cached = self._cache_get(url)
        if cached is not None:
            return cached
        text = self._request(url, "text/csv").decode("utf-8")
        self._cache_put(url, text)
        return text

    def _get(self, path: str, params: dict | None = None):
        url = f"{self.base_url}/{path.lstrip('/')}"
        if params:
            clean = {k: v for k, v in params.items() if v is not None}
            if clean:
                url += "?" + urllib.parse.urlencode(clean, doseq=True)
        cached = self._cache_get(url)
        if cached is not None:
            return cached
        data = json.loads(self._request(url, "application/json").decode("utf-8"))
        self._cache_put(url, data)
        return data

    def get_study(self, study_id: str) -> dict:
        """Study metadata: id, name, projectId, description, publishedAt, experiments[{id, name}]."""
        return self._get(f"study/{study_id}.json")

    def get_experiment(self, experiment_id: str) -> dict:
        """Experiment: cultivationMode, communityStrains[], bioreplicates[] (each with measurementContexts)."""
        return self._get(f"experiment/{experiment_id}.json")

    def get_bioreplicate(self, bioreplicate_id) -> dict:
        record = self._get(f"bioreplicate/{bioreplicate_id}.json")
        for context in record.get("measurementContexts", []):
            self._owner[str(context.get("id"))] = bioreplicate_id
        return record

    def get_replicate_series(self, bioreplicate_id) -> dict:
        """{measurement context id: [(time, value, std or None), ...]} for every series of one bioreplicate,
        from its one CSV (`bioreplicate/<id>.csv`), in time order: one request where reading each series
        alone takes one per series (Karoline, 2026-10-09, after a slow mGrowthDB: 53 requests instead of 542
        for the CRM example; a review found the values those of the series' own CSVs in every unit served)."""
        key = str(bioreplicate_id)
        if key in self._series:
            return self._series[key]
        text = self._get_text(f"bioreplicate/{bioreplicate_id}.csv")
        series = {}
        for row in csv.DictReader(io.StringIO(text)):
            point = _point(row)
            if point is not None and row.get("measurementContextId"):
                series.setdefault(str(row["measurementContextId"]).strip(), []).append(point)
        series = {cid: sorted(points) for cid, points in series.items()}
        self._series[key] = series
        return series

    def get_measurement_context(self, context_id) -> dict:
        """A single measurement context: techniqueType, subject, auc, growthRate, units, ..."""
        return self._get(f"measurement-context/{context_id}.json")

    def get_measurement_series(self, context_id) -> list:
        """The measured time series of one measurement context: [(time, value, std or None), ...].

        mGrowthDB serves the points as CSV (`measurement-context/<id>.csv`, columns time, value, std); the
        JSON representation carries only the summarized growthRate and auc. Rows without a readable time
        or value are dropped, and the points are returned in time order.
        """
        owner = self._owner.get(str(context_id))
        if owner is not None and self._series.get(str(owner), {}) is not None:
            try:
                series = self.get_replicate_series(owner)
            except MGrowthDBError:
                series = self._series[str(owner)] = None     # this replicate's series one by one, asked once
            if series and str(context_id) in series:
                return list(series[str(context_id)])
        text = self._get_text(f"measurement-context/{context_id}.csv")
        points = [p for p in (_point(row) for row in csv.DictReader(io.StringIO(text))) if p is not None]
        return sorted(points)

    def search(self, strain_ncbi_ids=None, metabolite_chebi_ids=None) -> dict:
        return self._get("search.json", {
            "strainNcbiIds": strain_ncbi_ids, "metaboliteChebiIds": metabolite_chebi_ids,
        })

    def chebi_compound(self, chebi_id) -> dict | None:
        """ChEBI's formula, charge and mass of one compound (foodnet.chemistry), each id read once per client."""
        if self._chebi is None:
            from .chemistry import Memo
            self._chebi = Memo()
        return self._chebi(chebi_id)

    def study_experiments(self, study_id: str) -> list:
        """Full experiment records for a study (study metadata lists experiment ids only)."""
        study = self.get_study(study_id)
        return [self.get_experiment(e["id"]) for e in study.get("experiments", [])]


def _point(row: dict):
    """(time, value, std or None) of a CSV row, or None for a row without a readable time or value."""
    try:
        time_point, value = float(row["time"]), float(row["value"])
    except (TypeError, ValueError, KeyError):
        return None
    try:
        std = float(row.get("std") or "")
    except ValueError:
        std = None
    return time_point, value, std


def noting_slowness(progress, client):
    """`progress` with the client's slowness, if any, after each message, so a slow mGrowthDB does not look like
    a stalled search."""
    if progress is None or not hasattr(client, "slowness"):
        return progress

    def say(done, total, message):
        note = client.slowness()
        progress(done, total, f"{message}. {note}" if note else message)
    return say


# ---- what a network was made from ------------------------------------------------------------

def provenance(today: datetime.date | None = None, now: datetime.datetime | None = None) -> dict:
    """What made a network and when: the tool, its version, and the derivation date and time (local, with
    its offset). mGrowthDB changes over time, so the same version can derive a different network later;
    which data it saw is in `data_versions`: each study's uploadedAt, and since when what was used was read."""
    now = now or datetime.datetime.now().astimezone()
    return {"tool": NAME, "tool_version": __version__,
            "derived_on": (today or now.date()).isoformat(),
            "derived_at": now.isoformat(timespec="seconds")}


# mGrowthDB publishes no version of the database as a whole (its API has no version endpoint), so the
# version of the data behind a network is when it was read plus each study's own upload and publication
# dates, which change when a study is corrected.
NO_DATABASE_VERSION = "mGrowthDB publishes no database version; each study's upload and publication dates are given"


def data_versions(client, study_ids, retrieved_at: str, species_list=None) -> dict:
    """The data a search read: the API, when, and for each study its uploadedAt and publishedAt, and, for a study
    whose experiments, replicates and series came from foodnet's local copy (foodnet.store), when they were
    first read (`kept_since`); and when the species list that found the studies was read (a review: a result
    said "live" and gave the search's time for data read days before)."""
    kept = getattr(client, "kept_studies", {}) or {}
    studies = {}
    for sid in study_ids:
        try:
            study = client.get_study(sid)          # cached by the client, so no second request
        except MGrowthDBError:
            continue
        studies[sid] = {"uploaded_at": study.get("uploadedAt", ""), "published_at": study.get("publishedAt", "")}
        if sid in kept:
            studies[sid]["kept_since"] = _iso(kept[sid])
    out = {"api": MGROWTHDB_API, "retrieved_at": retrieved_at, "database_version": NO_DATABASE_VERSION,
           "studies": studies}
    if getattr(client, "searches", 0) > 1:
        # what the client read for earlier searches is used again: on the page up to an hour of them (reviews)
        out["read_since"] = _iso(client.started_at)
    if species_list is not None:
        out["species_list"] = {"read_at": _iso(species_list.built_at),
                               "kept": bool(getattr(species_list, "kept", False))}
    return out


def _iso(seconds: float) -> str:
    return datetime.datetime.fromtimestamp(seconds).astimezone().isoformat(timespec="seconds")
