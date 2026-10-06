"""GraphML, the Cytoscape payload and style, the report, and the schema contract."""
import json
import os
import xml.etree.ElementTree as ET

from conftest import run

from foodnet import cytoscape
from foodnet.export import to_graphml
from foodnet.report import report_text
from foodnet.schema import SCHEMA_DOC, SCHEMA_FILE, schema_json, validate_document

NS = "{http://graphml.graphdrawing.org/xmlns}"


def test_graphml_keeps_kinds_directions_and_leaves_missing_numbers_out(client):
    net = run(client)["network"]
    root = ET.fromstring(to_graphml(net))
    keys = {k.get("id"): k.get("attr.name") for k in root.iter(f"{NS}key")}
    edges = root.find(f"{NS}graph").findall(f"{NS}edge")
    assert len(edges) == len(net.edges)
    presence = [e for e in edges if any(d.text == "presence_only" for d in e)]
    assert presence and not any(keys[d.get("key")] == "amount" for e in presence for d in e)


def test_cytoscape_gets_no_null_numbers_and_one_interaction_per_kind(client):
    payload = cytoscape.network_json(run(client, phase="both")["network"])
    for edge in payload["elements"]["edges"]:
        assert None not in edge["data"].values()
        assert edge["data"]["interaction"].split()[0] in ("produced", "consumed")
    kinds = {n["data"]["kind"] for n in payload["elements"]["nodes"]}
    assert kinds == {"taxon", "metabolite"}


def test_the_style_draws_the_two_kinds_and_the_two_directions():
    style = cytoscape.style()
    maps = {m["visualProperty"]: m for m in style["mappings"]}
    assert {e["key"]: e["value"] for e in maps["NODE_SHAPE"]["map"]} == {"taxon": "ELLIPSE",
                                                                         "metabolite": "ROUND_RECTANGLE"}
    colors = {e["key"]: e["value"] for e in maps["EDGE_STROKE_UNSELECTED_PAINT"]["map"]}
    assert colors == {"produced": "#2160A8", "consumed": "#B45309"}
    assert "<vizmap" in cytoscape.style_xml()


def test_the_report_holds_settings_media_and_reasons(client):
    text = report_text(run(client))
    assert "Values from the medium holding data for the most taxa" in text
    assert "detection_limit: 0.2" in text and "Sources (cited at the arc level" in text


def test_the_shipped_schema_is_the_code_and_a_network_validates(client):
    with open(SCHEMA_FILE, encoding="utf-8") as f:
        assert json.load(f) == SCHEMA_DOC, "regenerate it: make schema"
    doc = json.loads(run(client, phase="both")["network"].to_json())
    assert validate_document(doc) == []
    assert schema_json().startswith("{")


def test_the_validator_finds_what_is_wrong(client):
    doc = json.loads(run(client)["network"].to_json())
    doc["edges"][0]["direction"] = "sideways"
    doc["edges"][1]["study_ids"] = []
    problems = validate_document(doc)
    assert any("sideways" in p for p in problems) and any("at least 1" in p for p in problems)


def test_the_example_fixture_is_valid():
    path = os.path.join(os.path.dirname(__file__), "fixtures", "example_network.json")
    with open(path, encoding="utf-8") as f:
        assert validate_document(json.load(f)) == []


def test_a_taxon_with_a_growth_rate_but_no_arc_is_reported_by_name(client):
    # found in the pre-release audit: with values from mMCB only and the stationary phase, Beta beta has a
    # growth rate and no arc, and the report looked its name up among the network's nodes
    r = run(client, phase="stationary", booleans=True, conditions="mMCB", report_rates=True, rate_window=3)
    assert "ncbi:2" not in r["network"].nodes
    text = report_text(r)
    assert "Beta beta B1:" in text
