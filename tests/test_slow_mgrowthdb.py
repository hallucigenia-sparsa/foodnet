"""Reaching mGrowthDB less and faster (Karoline, 2026-10-09, after a search took over 8 minutes: every new
connection waited 17 s for an IPv6 address that did not answer): the address fallback, the slowness note, a
replicate's series from one CSV, the species list kept a day and what was read of a study kept while it is
unchanged."""
import json
import time

import pytest

from foodnet import mgrowthdb, taxonomy


def test_an_address_that_does_not_answer_is_left_after_a_moment(monkeypatch):
    # 2026-10-09: mGrowthDB's IPv6 address did not answer, and each connection waited 17 s for it
    import socket
    tried = []

    class Sock:
        def __init__(self, family, kind, proto):
            self.family = family

        def settimeout(self, t):
            self.timeout = t

        def connect(self, address):
            tried.append((address[0], self.timeout))
            if address[0] == "v6":
                raise TimeoutError("timed out")

        def close(self):
            pass
    monkeypatch.setattr(socket, "getaddrinfo", lambda host, port, type: [
        (socket.AF_INET6, socket.SOCK_STREAM, 0, "", ("v6", port)),
        (socket.AF_INET, socket.SOCK_STREAM, 0, "", ("v4", port))])
    monkeypatch.setattr(socket, "socket", Sock)
    preferred = []
    sock = mgrowthdb.open_socket("host", 443, 30, preferred)
    quick = mgrowthdb.CONNECT_TIMEOUT
    assert sock.family == socket.AF_INET and tried == [("v6", quick), ("v4", quick)]
    assert preferred == [socket.AF_INET]
    tried.clear()
    mgrowthdb.open_socket("host", 443, 30, preferred)
    assert tried == [("v4", quick)]                                # the family that answered goes first
    tried.clear()
    assert mgrowthdb.open_socket("host", 443, None, []).family == socket.AF_INET   # no limit: no crash


def test_a_machine_without_ipv6_and_a_slow_path_still_connect(monkeypatch):
    import socket
    tried = []

    class Sock:
        def __init__(self, family, kind, proto):
            if family == socket.AF_INET6:
                raise OSError(47, "Address family not supported by protocol")
            self.family = family

        def settimeout(self, t):
            self.timeout = t

        def connect(self, address):
            tried.append(self.timeout)
            if self.timeout <= mgrowthdb.CONNECT_TIMEOUT:
                raise TimeoutError("slow path")                    # answers only with more time

        def close(self):
            pass
    monkeypatch.setattr(socket, "getaddrinfo", lambda host, port, type: [
        (socket.AF_INET6, socket.SOCK_STREAM, 0, "", ("v6", port)),
        (socket.AF_INET, socket.SOCK_STREAM, 0, "", ("v4", port))])
    monkeypatch.setattr(socket, "socket", Sock)
    sock = mgrowthdb.open_socket("host", 443, 30, [socket.AF_INET])
    assert sock.family == socket.AF_INET and tried == [mgrowthdb.CONNECT_TIMEOUT, 30]
    # the slow working address second, after a dead one (a review): it still gets the full time
    tried.clear()
    sock = mgrowthdb.open_socket("host", 443, 30, [])
    assert sock.family == socket.AF_INET and tried == [mgrowthdb.CONNECT_TIMEOUT, 30]
    monkeypatch.setattr(socket, "getaddrinfo", lambda host, port, type: [
        (socket.AF_INET6, socket.SOCK_STREAM, 0, "", ("v6", port))])
    with pytest.raises(OSError, match="v6: .*not supported"):
        mgrowthdb.open_socket("host", 443, 30, [])


class _Slow(mgrowthdb.MGrowthDBClient):
    """A client whose server answers after `delay` seconds, failing with 503 the first `fail` times."""

    def __init__(self, fail=0):
        super().__init__(backoff=0)
        self.fail = fail

    def _send(self, url, accept):
        if self.fail:
            self.fail -= 1
            raise mgrowthdb._Status(503, "Service Unavailable")
        return b'{"id": "x"}'


def test_slowness_is_counted_and_said(monkeypatch):
    client = _Slow(fail=1)
    assert client.slowness() == ""
    client._get("study/x.json")
    assert client.slowness() == "mGrowthDB is answering slowly: 1 retried"
    clock = iter([0.0, mgrowthdb.SLOW + 2])
    monkeypatch.setattr(mgrowthdb.time, "monotonic", lambda: next(clock))
    client._get("study/y.json")
    assert "1 request(s) took over 10 s" in client.slowness()


def test_the_progress_text_carries_the_slowness():
    client = _Slow(fail=1)
    seen = []
    say = mgrowthdb.noting_slowness(lambda d, t, m: seen.append(m), client)
    say(0, None, "Reading studies")
    client._get("study/x.json")
    say(1, None, "Reading studies")
    assert seen == ["Reading studies", "Reading studies. mGrowthDB is answering slowly: 1 retried"]


class _Crawl:
    """The reads the species-list crawl makes, counted; `base_url` makes the list keepable."""

    base_url = "https://example.test/api/v1"

    def __init__(self, failed=False):
        self.calls = 0
        self.failed = failed

    def get_study(self, sid):
        self.calls += 1
        if sid == "SMGDB00000001":
            return {"id": sid, "publishedAt": "2026", "experiments": [{"id": "E1"}]}
        raise mgrowthdb.MGrowthDBError("no such study", status=404)

    def get_experiment(self, eid):
        if self.failed:
            raise mgrowthdb.MGrowthDBError("down", status=503)
        return {"id": eid, "communityStrains": [{"name": "Alpha alpha A1", "NCBId": 11}]}

    def study_experiments(self, sid):
        return [self.get_experiment(e["id"]) for e in self.get_study(sid)["experiments"]]


@pytest.fixture
def cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FOODNET_CACHE_DIR", str(tmp_path))
    return tmp_path


def test_the_species_list_is_kept_for_a_day_and_read_back(cache):
    first = _Crawl()
    index = taxonomy.kept_species_index(first)
    assert index["alpha alpha"] == {11: "Alpha alpha A1"} and first.calls > 0
    again = _Crawl()
    kept = taxonomy.kept_species_index(again)
    assert again.calls == 0                                          # read from the cache folder
    assert kept["alpha alpha"] == {11: "Alpha alpha A1"} and kept.studies_of([11]) == ["SMGDB00000001"]
    assert kept.current == {11: "Alpha alpha A1"}
    refreshed = _Crawl()
    taxonomy.kept_species_index(refreshed, refresh=True)
    assert refreshed.calls > 0


def test_an_old_or_foreign_or_incomplete_list_is_not_used(cache):
    taxonomy.kept_species_index(_Crawl())
    path = cache / "species_list.json"
    kept = json.loads(path.read_text())
    kept["built_at"] -= taxonomy.KEEP_FOR + 1
    path.write_text(json.dumps(kept))
    old = _Crawl()
    taxonomy.kept_species_index(old)
    assert old.calls > 0                                             # a day old: crawled again
    other = _Crawl()
    other.base_url = "https://mirror.test/api/v1"
    taxonomy.kept_species_index(other)
    assert other.calls > 0                                           # another mGrowthDB: its own list
    path.unlink()
    broken = taxonomy.kept_species_index(_Crawl(failed=True))
    assert broken.failed and not path.exists()                       # an incomplete crawl is never kept


def test_a_client_without_a_base_url_is_never_kept(cache):
    client = _Crawl()
    del_url = type("NoUrl", (_Crawl,), {"base_url": None})()
    taxonomy.kept_species_index(del_url)
    assert not (cache / "species_list.json").exists() and client.calls == 0


class _Recorded(mgrowthdb.MGrowthDBClient):
    """A client answered from a table of url -> body, counting what it asked for."""

    def __init__(self, answers):
        super().__init__(base_url="https://example.test/api/v1", backoff=0)
        self.answers, self.asked = answers, []

    def _send(self, url, accept):
        self.asked.append(url.split("/api/v1/")[-1])
        if url not in self.answers:
            raise mgrowthdb._Status(404, "Not Found")
        body = self.answers[url]
        return body.encode() if isinstance(body, str) else json.dumps(body).encode()


BASE = "https://example.test/api/v1"
REPLICATE_CSV = ("measurementContextId,subjectType,subjectName,subjectExternalId,time,value,std\n"
                 "7,bioreplicate,r1,,0,0.1,\n7,bioreplicate,r1,,4,0.5,\n8,metabolite,acetate,30089,4,2.5,0.1\n"
                 "8,metabolite,acetate,30089,0,1.0,\n")


def _answers(uploaded="2026-01-01"):
    return {f"{BASE}/study/S1.json": {"id": "S1", "uploadedAt": uploaded, "experiments": [{"id": "E1"}]},
            f"{BASE}/experiment/E1.json": {"id": "E1", "bioreplicates": [{"id": 5}]},
            f"{BASE}/bioreplicate/5.json": {"id": 5, "measurementContexts": [{"id": 7}, {"id": 8}]},
            f"{BASE}/bioreplicate/5.csv": REPLICATE_CSV}


def test_a_series_is_read_from_its_replicates_one_csv():
    client = _Recorded(_answers())
    client.get_bioreplicate(5)
    assert client.get_measurement_series(8) == [(0.0, 1.0, None), (4.0, 2.5, 0.1)]
    assert client.get_measurement_series(7) == [(0.0, 0.1, None), (4.0, 0.5, None)]
    assert client.asked == ["bioreplicate/5.json", "bioreplicate/5.csv"]          # one CSV for both series


def test_what_was_read_of_an_unchanged_study_is_used_again(cache):
    from foodnet import store
    first = _Recorded(_answers())
    study = first.get_study("S1")
    first.get_experiment("E1")
    first.get_bioreplicate(5)
    first.get_replicate_series(5)
    assert store.save(first, study) == 3
    again = _Recorded(_answers())
    assert store.load(again, again.get_study("S1")) == 3
    again.get_experiment("E1")
    again.get_bioreplicate(5)
    assert again.get_measurement_series(8)[-1] == (4.0, 2.5, 0.1)
    assert again.asked == ["study/S1.json"]                                       # the study record only
    changed = _Recorded(_answers(uploaded="2026-02-02"))
    assert store.load(changed, changed.get_study("S1")) == 0                      # resubmitted: read again
    bare = _Recorded({**_answers(), f"{BASE}/study/S1.json": {"id": "S1", "experiments": [{"id": "E1"}]}})
    assert store.save(bare, bare.get_study("S1")) == 0                            # no uploadedAt: never kept


def test_past_the_last_known_study_fewer_absent_ids_end_the_crawl():
    from foodnet.fetch import study_ids_in_order
    from foodnet.taxonomy import MISS_AFTER_KNOWN

    class Studies:
        asked = []

        def get_study(self, sid):
            Studies.asked.append(sid)
            if sid in ("S1", "S4"):
                return {"id": sid}
            raise mgrowthdb.MGrowthDBError("no", status=404)
    assert study_ids_in_order(Studies(), "S{}", 500, 25, known_last=4) == ["S1", "S4"]
    assert len(Studies.asked) <= 4 + MISS_AFTER_KNOWN + 6                         # a batch of six past the stop


def test_the_about_page_lists_every_release_of_the_changelog():
    # Karoline, 2026-10-09: "include the release history in the About page"
    import pathlib
    import re as _re

    from foodnet import help as help_page
    changelog = (pathlib.Path(__file__).resolve().parents[1] / "CHANGELOG.md").read_text()
    released = _re.findall(r"^## \[(\d+\.\d+\.\d+)\] \((\d{4}-\d{2}-\d{2})\)", changelog, _re.M)
    assert [(v, d) for v, d, _ in help_page.RELEASES] == released
    page = help_page.render_about()
    assert "Release history" in page and all(f"<b>{v}</b> ({d})" in page for v, d in released)


def test_taxa_on_the_command_line_can_be_separated_by_commas_or_semicolons():
    # Karoline, 2026-10-09: "allow taxon names on command line to be separated with a non-blank delimiter since
    # strain names include blanks"
    from foodnet.__main__ import taxa_entries
    assert taxa_entries(["Escherichia", "coli", "LF82,", "Bacteroides", "fragilis"]) == \
        ["Escherichia coli LF82", "Bacteroides fragilis"]
    assert taxa_entries(["Escherichia coli LF82;Roseburia intestinalis L1-82"]) == \
        ["Escherichia coli LF82", "Roseburia intestinalis L1-82"]
    assert taxa_entries(["Blautia", "Roseburia"]) == ["Blautia", "Roseburia"]           # as before
    assert taxa_entries(["Escherichia coli LF82"]) == ["Escherichia coli LF82"]         # a quoted name
    # a quoted name stays one of its own beside commas (a review: it was joined to its neighbours)
    assert taxa_entries(["Escherichia coli LF82", "Bacteroides fragilis,", "Roseburia intestinalis"]) == \
        ["Escherichia coli LF82", "Bacteroides fragilis", "Roseburia intestinalis"]



def test_the_store_works_through_a_search_read(cache):
    # the path a search takes (foodnet.fetch.prefetch_studies), not only the store's own functions (a review)
    from foodnet.fetch import prefetch_studies
    answers = {**_answers(), f"{BASE}/experiment/E1.json": {
        "id": "E1", "cultivationMode": "batch", "bioreplicates": [{"id": 5}],
        "communityStrains": [{"id": 1, "NCBId": 11, "name": "Alpha alpha A1"}]}}
    first = _Recorded(answers)
    prefetch_studies(first, ["S1"])
    assert "bioreplicate/5.csv" in first.asked
    again = _Recorded(answers)
    prefetch_studies(again, ["S1"])
    assert again.asked == ["study/S1.json"] and "S1" in again.kept_studies
    assert again.get_measurement_series(8)[-1] == (4.0, 2.5, 0.1)
    fresh = _Recorded(answers)
    fresh.refresh = True
    prefetch_studies(fresh, ["S1"])
    assert "bioreplicate/5.csv" in fresh.asked                      # a refresh reads it all again


def test_a_kept_study_is_read_again_after_thirty_days(cache, monkeypatch):
    from foodnet import store
    first = _Recorded(_answers())
    study = first.get_study("S1")
    first.get_experiment("E1")
    assert store.save(first, study, "S1") == 1
    later = time.time() + store.MAX_AGE + 1
    monkeypatch.setattr(store.time, "time", lambda: later)
    assert store.load(_Recorded(_answers()), study, "S1") == 0


def test_a_broken_or_strange_store_never_breaks_a_search(cache):
    from foodnet import store
    (cache / "studies").mkdir(parents=True)
    (cache / "studies" / "S1.json").write_text(json.dumps(
        {"format": store.FORMAT, "source": BASE, "uploadedAt": "2026-01-01", "read_at": time.time(),
         "responses": ["a"]}))
    client = _Recorded(_answers())
    study = client.get_study("S1")
    assert store.load(client, study, "S1") == 0                     # read as absent, no crash
    client.get_experiment("E1")
    assert store.save(client, study, "S1") == 1                     # and replaced by a sound file
    assert store.save(client, {**study, "id": "../x"}) == 0          # an id that is no file name is never used
    assert not (cache / "x.json").exists()


def test_a_replicate_csv_that_fails_falls_back_to_each_series_once():
    answers = {**_answers(), f"{BASE}/measurement-context/8.csv": "time,value,std\n0,1.0,\n4,2.5,0.1\n"}
    del answers[f"{BASE}/bioreplicate/5.csv"]
    client = _Recorded(answers)
    client.get_bioreplicate(5)
    assert client.get_measurement_series(8) == [(0.0, 1.0, None), (4.0, 2.5, 0.1)]
    assert client.asked.count("bioreplicate/5.csv") == 1


def test_the_data_versions_say_what_came_from_the_copy_and_how_old_the_species_list_is(cache):
    from foodnet import store
    first = _Recorded(_answers())
    study = first.get_study("S1")
    first.get_experiment("E1")
    store.save(first, study, "S1")
    again = _Recorded(_answers())
    store.load(again, again.get_study("S1"), "S1")
    index = taxonomy.SpeciesIndex({})
    index.kept = True
    versions = mgrowthdb.data_versions(again, ["S1"], "now", index)
    assert versions["studies"]["S1"]["kept_since"] and versions["species_list"]["kept"] is True



def test_a_refresh_in_the_pages_order_reads_everything_again(cache):
    # the page makes the species list before the search sets refresh on the client (a review: the crawl then
    # loaded every kept study into memory, and the refresh read them from there)
    from foodnet.fetch import prefetch_studies

    class Crawlable(_Recorded):
        def study_experiments(self, sid):
            return [self.get_experiment(e["id"]) for e in self.get_study(sid)["experiments"]]
    answers = {**_answers(), f"{BASE}/experiment/E1.json": {
        "id": "E1", "cultivationMode": "batch", "bioreplicates": [{"id": 5}],
        "communityStrains": [{"id": 1, "NCBId": 11, "name": "Alpha alpha A1"}]}}
    answers[f"{BASE}/study/SMGDB00000001.json"] = answers.pop(f"{BASE}/study/S1.json")
    answers[f"{BASE}/study/SMGDB00000001.json"]["id"] = "SMGDB00000001"
    first = Crawlable(answers)
    taxonomy.kept_species_index(first)
    prefetch_studies(first, ["SMGDB00000001"])
    page = Crawlable(answers)
    taxonomy.kept_species_index(page, refresh=True)                  # the page's order: the list first
    prefetch_studies(page, ["SMGDB00000001"])
    assert "bioreplicate/5.csv" in page.asked and not page.kept_studies



def test_a_refresh_replaces_the_copy_and_an_expired_copy_is_not_saved_again(cache, monkeypatch):
    # reviews: a refresh read everything again but kept the old file; an expired copy still in memory was saved
    # again with today's date
    from foodnet import store
    old = _Recorded(_answers())
    study = old.get_study("S1")
    old.get_experiment("E1")
    store.save(old, study, "S1")
    corrected = {**_answers(), f"{BASE}/experiment/E1.json": {"id": "E1", "bioreplicates": [{"id": 6}]}}
    fresh = _Recorded(corrected)
    fresh.refresh = True
    fresh.get_experiment("E1")
    store.save(fresh, fresh.get_study("S1"), "S1")
    after = _Recorded(_answers())
    store.load(after, after.get_study("S1"), "S1")
    assert after.get_experiment("E1")["bioreplicates"] == [{"id": 6}]          # the refreshed values
    page = _Recorded(_answers())
    store.load(page, page.get_study("S1"), "S1")                                 # in memory from the copy
    later = time.time() + store.MAX_AGE + 1
    monkeypatch.setattr(store.time, "time", lambda: later)
    assert store.save(page, page.get_study("S1"), "S1") == 0                      # not stamped as newly read



def test_a_page_client_never_writes_back_what_it_took_from_an_older_copy(cache):
    # a review: after a refresh, the page's client still held the old values and its next search saved them back
    from foodnet import store
    old = _Recorded(_answers())
    study = old.get_study("S1")
    old.get_experiment("E1")
    store.save(old, study, "S1")
    page = _Recorded(_answers())
    store.load(page, page.get_study("S1"), "S1")                    # the page's client holds the old copy
    corrected = {**_answers(), f"{BASE}/experiment/E1.json": {"id": "E1", "bioreplicates": [{"id": 6}]}}
    fresh = _Recorded(corrected)
    fresh.refresh = True
    fresh.get_experiment("E1")
    store.save(fresh, fresh.get_study("S1"), "S1")
    store.load(page, page.get_study("S1"), "S1")                    # its next search takes the new copy
    assert page.get_experiment("E1")["bioreplicates"] == [{"id": 6}]
    store.save(page, page.get_study("S1"), "S1")
    after = _Recorded(_answers())
    store.load(after, after.get_study("S1"), "S1")
    assert after.get_experiment("E1")["bioreplicates"] == [{"id": 6}]   # the refresh stands


def test_a_connection_that_could_not_be_made_is_not_made_twice_per_try():
    # a review: the reopen meant for a kept-open connection the server closed doubled the wait of an outage
    opened = []

    class Conn:
        sock = None

        def request(self, *a, **k):
            raise TimeoutError("could not connect")

    class Client(mgrowthdb.MGrowthDBClient):
        def _connection(self, fresh=False):
            opened.append(fresh)
            return Conn()
    with pytest.raises(TimeoutError):
        Client()._send("https://example.test/api/v1/study/S1.json", "application/json")
    assert opened == [False]



def test_a_refresh_by_another_process_is_not_undone_by_an_older_client(cache):
    # a review: the page's hour-old client merged what it had read over a command-line refresh
    from foodnet import store
    page = _Recorded(_answers())
    study = page.get_study("S1")
    page.get_experiment("E1")
    store.save(page, study, "S1")                                   # the page's client keeps what it read
    time.sleep(0.01)
    corrected = {**_answers(), f"{BASE}/experiment/E1.json": {"id": "E1", "bioreplicates": [{"id": 6}]}}
    other = _Recorded(corrected)                                    # another process, started later, refreshes
    other.refresh = True
    other.get_experiment("E1")
    store.save(other, other.get_study("S1"), "S1")
    page.get_bioreplicate(5)                                        # the page reads something more, and saves
    assert store.save(page, study, "S1") == 0
    after = _Recorded(_answers())
    store.load(after, after.get_study("S1"), "S1")
    assert after.get_experiment("E1")["bioreplicates"] == [{"id": 6}]



def test_after_a_refresh_elsewhere_the_older_client_takes_the_new_copy(cache):
    # a review: the page kept serving its own older reads and parsed series after a terminal refresh
    from foodnet import store
    page = _Recorded(_answers())
    study = page.get_study("S1")
    page.get_experiment("E1")
    page.get_bioreplicate(5)
    assert page.get_measurement_series(8)[-1] == (4.0, 2.5, 0.1)
    store.save(page, study, "S1")
    time.sleep(0.01)
    newer = REPLICATE_CSV.replace("4,2.5,0.1", "4,9.9,0.1")
    other = _Recorded({**_answers(), f"{BASE}/bioreplicate/5.csv": newer})
    other.refresh = True
    other.get_experiment("E1")
    other.get_bioreplicate(5)
    other.get_replicate_series(5)
    store.save(other, other.get_study("S1"), "S1")
    store.load(page, study, "S1")
    assert page.get_measurement_series(8)[-1] == (4.0, 9.9, 0.1)



def test_a_partial_refresh_elsewhere_leaves_no_old_replicate_in_memory(cache):
    # a review: a refresh that read one of two replicates left the other old in the page's memory
    from foodnet import store
    two = {**_answers(), f"{BASE}/experiment/E1.json": {"id": "E1", "bioreplicates": [{"id": 5}, {"id": 6}]},
           f"{BASE}/bioreplicate/6.json": {"id": 6, "measurementContexts": [{"id": 9}]},
           f"{BASE}/bioreplicate/6.csv": REPLICATE_CSV.replace("7,", "9,").replace("8,", "9,")}
    page = _Recorded(two)
    study = page.get_study("S1")
    page.get_experiment("E1")
    for rid in (5, 6):
        page.get_bioreplicate(rid)
        page.get_replicate_series(rid)
    store.save(page, study, "S1")
    time.sleep(0.01)
    other = _Recorded(two)
    other.refresh = True
    other.get_experiment("E1")
    other.get_bioreplicate(5)
    other.get_replicate_series(5)
    store.save(other, other.get_study("S1"), "S1")                   # a copy without replicate 6
    store.load(page, study, "S1")
    asked = len(page.asked)
    page.get_replicate_series(6)
    assert page.asked[asked:] == ["bioreplicate/6.csv"]              # read again, not served old
