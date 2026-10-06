"""A result built on records that could not be read says so: in its files, and in the exit status."""
import json

from conftest import FakeClient

from foodnet import matrix
from foodnet.mgrowthdb import MGrowthDBError
from foodnet.search import run_query
from foodnet.taxonomy import species_index

TAXA = ["Alpha alpha", "Beta beta", "Gamma gamma"]


class _Broken(FakeClient):
    """One experiment cannot be read."""

    def get_experiment(self, eid):
        if eid == "EMGDB000000004":
            raise MGrowthDBError("503 Service Unavailable", status=503)
        return super().get_experiment(eid)


def test_a_species_list_that_could_not_be_read_whole_says_so():
    index = species_index(_Broken())
    assert any(x == "EMGDB000000004" for x, _ in index.failed)
    result = run_query(_Broken(), TAXA, {}, index=index)
    assert result["errors"] and "incomplete" in result["errors"][0]
    meta = json.loads(result["network"].to_json())["meta"]
    assert meta["incomplete"] is True and meta["errors"]


def test_a_complete_result_says_it_is_complete(client):
    meta = run_query(client, TAXA, {})["network"].meta
    assert meta["incomplete"] is False and meta["errors"] == []


def test_the_command_line_exits_with_an_error_on_an_incomplete_result(monkeypatch, tmp_path):
    from foodnet import __main__ as cli
    from foodnet import mgrowthdb
    monkeypatch.setattr(mgrowthdb, "MGrowthDBClient", _Broken)
    out = tmp_path / "n.json"
    assert cli.main(["derive", "--taxa", *TAXA, "--out", str(out)]) == cli.INCOMPLETE
    assert json.loads(out.read_text())["meta"]["incomplete"] is True
    assert cli.main(["derive", "--taxa", *TAXA, "--out", str(out), "--allow-incomplete"]) == 0


def test_a_name_a_spreadsheet_would_run_is_written_inert():
    text = matrix._csv(["taxon", "=HYPERLINK(1)"], [["@SUM(A1)", "-1.5", "+x"]])
    assert "'=HYPERLINK(1)" in text and "'@SUM(A1)" in text and ",-1.5," in text and "'+x" in text


def test_an_incomplete_result_says_so_in_every_file(monkeypatch):
    import io
    import zipfile
    result = run_query(_Broken(), TAXA, {"report_rates": True}, index=species_index(_Broken()))
    assert matrix.signed_csv(result["network"], result).startswith("taxon [INCOMPLETE; ")
    z = zipfile.ZipFile(io.BytesIO(matrix.pair_package(result)))
    assert "INCOMPLETE" in z.read("README.txt").decode()
    payload = matrix.crm_payload(result)
    assert payload["caveats"]["incomplete"] is True and payload["caveats"]["errors"]
    assert "INCOMPLETE" in payload["readme"]


def test_an_incomplete_result_is_not_sent_on(monkeypatch, tmp_path):
    from foodnet import __main__ as cli
    from foodnet import mgrowthdb
    monkeypatch.setattr(mgrowthdb, "MGrowthDBClient", _Broken)
    assert cli.main(["derive", "--taxa", *TAXA, "--out", str(tmp_path / "n.json"), "--crm-mode", "--to-r"]) \
        == cli.INCOMPLETE
