"""Reading batch monocultures from mGrowthDB records (foodnet.reading)."""
from conftest import FakeClient

from foodnet import reading


def test_spellings_of_one_medium_share_a_key():
    # study 7 writes the abbreviation, study 9 does not; mucin as a second compartment medium stays apart
    assert reading.medium_key("Wilkins-Chalgren Anaerobe Broth (WC)") == reading.medium_key(
        "Wilkins-Chalgren Anaerobe Broth") == "wilkins chalgren anaerobe broth"
    assert reading.medium_key("Wilkins-Chalgren Anaerobe Broth (WC); Mucin") != \
        reading.medium_key("Wilkins-Chalgren Anaerobe Broth")


def test_cultures_hold_a_growth_curve_and_metabolites_in_hours_and_mm():
    found = reading.read_cultures(FakeClient(), ["SMGDB00000001"])
    cultures = found["cultures"]
    assert [c.replicate for c in cultures] == ["a1", "a2", "b1"]
    a1 = cultures[0]
    assert a1.taxon["id"] == "ncbi:1" and a1.medium_key == "wilkins chalgren anaerobe broth"
    assert a1.growth["values"] == [1, 10, 100, 1000, 1000]
    assert sorted(a1.metabolites) == ["chebi:17234", "chebi:30089"]     # acetic acid joined into acetate
    assert a1.metabolites["chebi:30089"]["recorded_as"] == "acetic acid (CHEBI:15366)"


def test_a_strain_not_asked_for_is_not_read():
    found = reading.read_cultures(FakeClient(), ["SMGDB00000001"], keep=lambda name, taxon: taxon == "2")
    assert {c.taxon["id"] for c in found["cultures"]} == {"ncbi:2"}


def test_an_average_replicate_and_a_non_batch_experiment_are_left_out_with_reasons():
    client = FakeClient()
    client.bioreplicates["a2"]["isAverage"] = True
    from conftest import EXPERIMENTS
    old = EXPERIMENTS["EMGDB000000002"]["cultivationMode"]
    EXPERIMENTS["EMGDB000000002"]["cultivationMode"] = "chemostat"
    try:
        found = reading.read_cultures(client, ["SMGDB00000001"])
    finally:
        EXPERIMENTS["EMGDB000000002"]["cultivationMode"] = old
    assert [c.replicate for c in found["cultures"]] == ["a1"]
    reasons = " ".join(r for _, r in found["skipped"])
    assert "average of the replicates" in reasons and "chemostat, not batch" in reasons


def test_a_series_in_a_unit_that_is_not_a_concentration_is_reported():
    client = FakeClient()
    for context in client.bioreplicates["a1"]["measurementContexts"]:
        if context["subject"].get("name") == "glucose":
            context["techniqueUnits"] = "AUC"
    found = reading.read_cultures(client, ["SMGDB00000001"])
    assert "chebi:17234" not in found["cultures"][0].metabolites
    assert any("'AUC'" in r for _, r in found["skipped"])


def test_a_metabolite_can_be_left_out_by_name():
    client = FakeClient()
    found = reading.read_cultures(client, ["SMGDB00000001"], excluded_metabolites=["Glucose"])
    assert all("chebi:17234" not in c.metabolites for c in found["cultures"])


def test_no_metabolite_is_left_out_by_default():
    # Karoline, 2026-10-04: "we don't want to skip any metabolites by default"
    from foodnet.search import DEFAULTS
    assert DEFAULTS["exclude_metabolites"] == ""
