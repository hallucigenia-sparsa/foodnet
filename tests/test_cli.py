"""The command line, on the synthetic world of conftest.py (the client is replaced, nothing is read live)."""
import io
import json
import zipfile

import pytest
from conftest import FakeClient

from foodnet import mgrowthdb
from foodnet.__main__ import main

TAXA = ["--taxa", "Alpha alpha", "Beta beta", "Gamma gamma"]


@pytest.fixture(autouse=True)
def fake_mgrowthdb(monkeypatch):
    monkeypatch.setattr(mgrowthdb, "MGrowthDBClient", FakeClient)


def test_derive_writes_a_valid_network(tmp_path, capsys):
    out = tmp_path / "net.json"
    assert main(["derive", *TAXA, "--out", str(out)]) == 0
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert doc["schema"] == "foodnet.metabolite_network/v1" and doc["edges"]
    assert main(["validate", str(out)]) == 0
    assert "valid" in capsys.readouterr().out


def test_every_format_and_the_side_files(tmp_path):
    assert main(["derive", *TAXA, "--phase", "both", "--format", "matrices", "--out", str(tmp_path / "m.zip"),
                 "--crm-mode", "--crm", str(tmp_path / "crm.zip"), "--rates", str(tmp_path / "r.csv"),
                 "--report", str(tmp_path / "report.txt")]) == 0
    assert "consumed.csv" in zipfile.ZipFile(tmp_path / "m.zip").namelist()
    assert "crm.json" in zipfile.ZipFile(io.BytesIO((tmp_path / "crm.zip").read_bytes())).namelist()
    assert (tmp_path / "r.csv").read_text().startswith("taxon,growth_rate")
    assert main(["derive", *TAXA, "--format", "matrix", "--out", str(tmp_path / "m.csv")]) == 0
    assert (tmp_path / "m.csv").read_text().startswith('taxon [values from Wilkins')
    assert main(["derive", *TAXA, "--format", "graphml", "--out", str(tmp_path / "n.graphml")]) == 0


def test_a_window_on_the_command_line(tmp_path):
    out = tmp_path / "w.json"
    assert main(["derive", *TAXA, "--window", "0", "6", "--out", str(out)]) == 0
    assert {e["phase"] for e in json.loads(out.read_text())["edges"]} == {"window"}


def test_crm_files_need_the_growth_rates(tmp_path, capsys):
    assert main(["derive", *TAXA, "--crm", str(tmp_path / "c.zip")]) == 2
    assert "--crm-mode" in capsys.readouterr().err


def test_a_path_that_cannot_be_written_is_one_line(tmp_path, capsys):
    missing = str(tmp_path / "no_such_folder" / "out.json")
    assert main(["derive", *TAXA, "--out", missing]) == 2
    err = capsys.readouterr().err
    assert f"foodnet: {missing}: No such file or directory" in err and "Traceback" not in err


def test_validate_reports_problems(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": "foodnet.metabolite_network/v1", "nodes": [], "studies": [],
                               "edges": [{"source": "a", "target": "b", "direction": "produced",
                                          "phase": "exponential", "evidence": "measured"}]}))
    assert main(["validate", str(bad)]) == 1
    assert "study_ids" in capsys.readouterr().out
    (tmp_path / "x.json").write_text("not json")
    assert main(["validate", str(tmp_path / "x.json")]) == 2
