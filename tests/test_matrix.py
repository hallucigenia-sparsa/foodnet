"""The two matrix formats and the CRM parameters (foodnet.matrix), on the synthetic world of conftest.py.

Karoline, 2026-10-04: "one matrix with taxa as rows and metabolites as columns and another with 2 matrices:
1 for consumption and the other for production."
"""
import csv
import io
import json
import zipfile

from conftest import run

from foodnet import matrix

A = "ncbi:1"
GLC, AC = "chebi:17234", "chebi:30089"


def _rows(text):
    """The CSV's rows, with the first header cell (which names the value medium and the version) as "taxon"."""
    rows = list(csv.reader(io.StringIO(text)))
    assert rows[0][0].startswith("taxon [values from ")
    rows[0][0] = "taxon"
    return rows


def test_the_signed_matrix_has_production_positive_consumption_negative_and_na_for_no_value(client):
    r = run(client)
    assert _rows(matrix.signed_csv(r["network"], r)) == [
        ["taxon", "acetate", "butyrate", "formate", "glucose"],
        ["Alpha alpha A1", "4", "NA", "NA", "-8"],
        # B's formate was seen only in mMCB: NA here, never 0
        ["Beta beta B1", "NA", "3", "NA", "-6"],
        # C was measured only in mMCB: nothing gives it a value
        ["Gamma gamma C1", "NA", "NA", "NA", "NA"]]


def test_a_measured_change_below_the_limit_is_zero_not_na(client):
    r = run(client, phase="stationary", judge_spread=False)
    rows = _rows(matrix.signed_csv(r["network"], r))
    # judged by the mean alone, A's stationary acetate (+1 and -1) is measured and no change, so 0; glucose -0.5
    assert rows[1] == ["Alpha alpha A1", "0", "NA", "NA", "-0.5"]


def test_a_spread_across_the_limit_is_inconclusive_never_zero(client):
    r = run(client, phase="stationary")
    rows = _rows(matrix.signed_csv(r["network"], r))
    # the default judges the spread too: acetate +1 and -1 (sd 1.41) and glucose 0 and -1 (sd 0.71) both reach
    # across the 0.2 mM limit, so they are NA, and their evidence says why
    assert rows[1] == ["Alpha alpha A1", "NA", "NA", "NA", "NA"]
    assert matrix.entry(r, A, AC, "stationary", "produced") == (None, "inconclusive")
    assert r["cells"][(A, GLC, "stationary")]["state"] == "inconclusive"
    assert not [e for e in r["network"].edges if e.taxon == A and e.phase == "stationary"]


def test_the_pair_holds_magnitudes_and_says_why_a_cell_is_na(client):
    r = run(client)
    z = zipfile.ZipFile(io.BytesIO(matrix.pair_package(r)))
    assert sorted(z.namelist()) == ["README.txt", "consumed.csv", "evidence_consumed.csv",
                                    "evidence_produced.csv", "matrices.svg", "produced.csv", "signed.csv"]
    consumed = _rows(z.read("consumed.csv").decode())
    produced = _rows(z.read("produced.csv").decode())
    assert consumed[1] == ["Alpha alpha A1", "0", "NA", "NA", "8"]       # acetate measured, went up: 0 here
    assert produced[1] == ["Alpha alpha A1", "4", "NA", "NA", "0"]
    evidence = _rows(z.read("evidence_produced.csv").decode())
    assert evidence[2] == ["Beta beta B1", "not_assayed", "measured", "presence_only", "below_limit"]


def test_booleans_count_presence_in_another_medium_as_one(client):
    r = run(client, booleans=True)
    assert _rows(matrix.signed_csv(r["network"], r))[1:] == [
        ["Alpha alpha A1", "1", "NA", "NA", "-1"],
        ["Beta beta B1", "NA", "1", "1", "-1"],
        ["Gamma gamma C1", "NA", "NA", "NA", "-1"]]


def test_both_phases_give_a_column_per_phase(client):
    r = run(client, phase="both")
    header = _rows(matrix.signed_csv(r["network"], r))[0]
    assert header[:3] == ["taxon", "acetate (exponential)", "acetate (stationary)"]


def test_the_crm_payload_carries_the_caveats_as_data(client):
    r = run(client, report_rates=True, rate_window=3)
    p = matrix.crm_payload(r)
    assert p["format"] == "foodnet.crm/v1" and p["phase"] == "exponential"
    assert p["taxa"] == ["Alpha alpha A1", "Beta beta B1", "Gamma gamma C1"]
    assert p["resources"] == ["acetate", "butyrate", "formate", "glucose"]
    assert p["consumed"][0] == [0, None, None, 8.0]
    assert p["initial_concentrations"] == [0.0, 0.0, None, 10.0]
    assert p["caveats"]["without_a_rate"] == ["Gamma gamma C1"]
    assert {(x["taxon"], x["resource"], x["direction"]) for x in p["caveats"]["presence_only"]} == {
        ("Beta beta B1", "formate", "produced"), ("Gamma gamma C1", "glucose", "consumed")}
    json.dumps(p)                                                     # it travels as JSON


def test_with_both_phases_the_crm_takes_the_exponential_one(client):
    r = run(client, phase="both", report_rates=True)
    p = matrix.crm_payload(r)
    assert p["phase"] == "exponential" and p["caveats"]["searched_both_phases"]
    assert "exponential phase is used" in p["readme"]


def test_the_crm_package_holds_every_file(client):
    r = run(client, report_rates=True)
    z = zipfile.ZipFile(io.BytesIO(matrix.crm_package(r)))
    assert sorted(z.namelist()) == ["README.txt", "biomass.csv", "consumed.csv", "crm.json", "evidence_consumed.csv",
                                    "evidence_produced.csv", "growth_rates.csv", "initial_concentrations.csv",
                                    "produced.csv"]


def test_the_crm_payload_carries_each_taxons_growth_over_the_phase(client):
    p = matrix.crm_payload(run(client, report_rates=True))
    a = p["taxa"].index("Alpha alpha A1")
    # A's flow cytometry counts rise from 1 at 0 h to 1000 at the 12 h boundary (conftest.py)
    assert p["biomass_change"][a] == 999 and p["biomass_start"][a] == 1 and p["phase_hours"][a] == 12
    # C has no culture in the value medium, so no growth there
    assert p["biomass_change"][p["taxa"].index("Gamma gamma C1")] is None
    assert p["caveats"]["mixed_media"] is False and p["caveats"]["stationary_phase"] is False
    assert matrix.crm_payload(run(client, report_rates=True, ignore_media=True))["caveats"]["mixed_media"]


def test_the_readme_lists_presence_the_matrices_lack(client):
    r = run(client)
    text = matrix.readme(r)
    assert "NA is never zero" in text
    assert "Beta beta B1 produced formate (exponential phase) in mMCB" in text


def test_second_window_columns_say_what_they_were_measured_over(client):
    r = run(client, second_window_metabolites="glucose")
    header = _rows(matrix.signed_csv(r["network"], r))[0]
    assert header == ["taxon", "acetate (exponential)", "butyrate (exponential)", "formate (exponential)",
                      "glucose (0 h to the last sample)"]
    assert _rows(matrix.signed_csv(r["network"], r))[1] == ["Alpha alpha A1", "4", "NA", "NA", "-8.5"]


def test_entries_seen_only_in_another_medium_are_na_true_or_their_value(client):
    # Karoline, 2026-10-06: NA "is a cautious default. I'd like to have an option to set them to TRUE and another
    # option to show the value but with a different background color in the image"
    # B's formate rose 5 mM in mMCB only; C's glucose fell 5 mM in mMCB only (conftest.py)
    na = run(client)
    assert _rows(matrix.signed_csv(na["network"], na))[2][3] == "NA"
    true = run(client, presence_entries="true")
    assert _rows(matrix.signed_csv(true["network"], true))[2][3] == "TRUE"
    assert _rows(matrix.signed_csv(true["network"], true))[3] == ["Gamma gamma C1", "NA", "NA", "NA", "TRUE"]
    value = run(client, presence_entries="value")
    assert _rows(matrix.signed_csv(value["network"], value))[2][3] == "5"
    assert _rows(matrix.signed_csv(value["network"], value))[3][4] == "-5"
    pair = matrix.pair_rows(value["network"], value)
    assert pair["consumed"][2][3] == 5.0 and pair["evidence_consumed"][2][3] == "presence_only"


def test_true_entries_reach_a_crm_as_missing(client):
    r = run(client, presence_entries="true", report_rates=True)
    p = matrix.crm_payload(r)
    assert p["produced"][1][2] is None and p["evidence_produced"][1][2] == "presence_only"
    assert "set to TRUE in the matrices" in p["readme"]
