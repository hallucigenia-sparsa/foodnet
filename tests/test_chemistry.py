"""The resources' chemistry from ChEBI, and each taxon's electron balance (Karoline, 2026-10-08)."""
import csv
import io
import zipfile

import pytest
from conftest import CHEBI_RECORDS, FakeClient, run

from foodnet import chemistry, compounds, matrix


def test_the_degree_of_reduction_is_the_same_for_an_acid_and_its_base():
    assert chemistry.degree_of_reduction("C2H3O2", -1) == (8, 2, None)           # acetate
    assert chemistry.degree_of_reduction("C2H4O2", 0) == (8, 2, None)            # acetic acid
    assert chemistry.degree_of_reduction("C4H4O4", -2)[0] == chemistry.degree_of_reduction("C4H6O4", 0)[0] == 14
    assert chemistry.degree_of_reduction("C6H12O6", 0)[0] == 24
    assert chemistry.degree_of_reduction("H2", 0) == (2, 0, None)
    assert chemistry.degree_of_reduction("CO2", 0)[0] == 0
    assert chemistry.degree_of_reduction("C3H7NO2", 0)[0] == 12                  # alanine: N as NH3
    assert chemistry.degree_of_reduction("H2S", 0)[0] == 8                       # S as H2SO4


def test_a_formula_that_cannot_be_evaluated_is_said_not_guessed():
    assert "Na" in chemistry.degree_of_reduction("C6H11NaO7", 0)[2]
    assert "variable" in chemistry.degree_of_reduction("(C6H10O5)n", 0)[2]
    assert "holds R," in chemistry.degree_of_reduction("C2H4O2R", 0)[2]          # a residue, not an element
    assert chemistry.degree_of_reduction("C2H3O2", None) == (None, 2, "ChEBI gives no charge")


def _fake(records):
    calls = []

    def fetch(chebi_id):
        calls.append(int(chebi_id))
        return records.get(int(chebi_id))
    return fetch, calls


def test_a_class_without_a_formula_takes_a_joined_forms_and_says_so():
    fetch, _ = _fake({26806: {"name": "succinate", "formula": "", "charge": None, "mass": None},
                      15741: {"name": "succinic acid", "formula": "C4H6O4", "charge": 0, "mass": "118.088"},
                      28757: {"name": "fructose", "formula": "", "charge": None, "mass": None},
                      15824: {"name": "D-fructose", "formula": "C6H12O6", "charge": 0, "mass": "180.156"}})
    succinate = chemistry._one("26806", fetch)
    assert (succinate["degree_of_reduction"], succinate["carbon"], succinate["formula_from"]) == (14, 4, "CHEBI:15741")
    assert "succinic acid" in succinate["note"]
    fructose = chemistry._one("28757", fetch)
    assert fructose["degree_of_reduction"] == 24 and fructose["formula_from"] == "CHEBI:15824"


def test_a_class_with_no_form_to_take_stays_blank():
    fetch, _ = _fake({28017: {"name": "starch", "formula": "", "charge": None, "mass": None}})
    starch = chemistry._one("28017", fetch)
    assert starch["degree_of_reduction"] is None and "no formula for starch" in starch["note"]
    assert chemistry._one("1", fetch)["note"] == "ChEBI does not hold this id"


def test_chebi_out_of_reach_leaves_every_resource_blank_with_the_reason():
    def down(_):
        raise OSError("timed out")
    found, problem = chemistry.resolve({"chebi:17234": "17234", "metabolite:x": ""}, down)
    assert found["chebi:17234"]["note"] == "ChEBI could not be reached" and "timed out" in problem
    assert found["metabolite:x"]["note"] == "mGrowthDB gives no ChEBI id"
    found, problem = chemistry.resolve({"chebi:17234": "17234"}, None)
    assert found["chebi:17234"]["degree_of_reduction"] is None and problem


def test_isovaleric_acid_joins_isovalerate_under_its_chebi_id():
    # foodnet 0.2.0 joined it to CHEBI:50128, which is biflavonoid; mGrowthDB records isovalerate as 48942
    node = compounds.canonical(28484, "isovaleric acid")
    assert node["id"] == "chebi:48942" and node["name"] == "isovalerate"
    assert compounds.canonical(48942, "isovalerate")["id"] == "chebi:48942"


def test_the_electron_balance_counts_what_has_a_number_and_a_degree_of_reduction():
    # 10 mM glucose (24) in; 15 mM acetate (8) and 5 mM butyrate (20) out; a resource without a number and one
    # without a degree of reduction are named, not counted
    b = chemistry.electron_balance([10, 0, 0, None, 3], [0, 15, 5, None, 0], [24, 8, 20, 8, None],
                                   [True] * 5)
    assert b["consumed_e_mM"] == 240 and b["produced_e_mM"] == 220
    assert b["share"] == pytest.approx(220 / 240)
    assert b["not_counted"] == [(3, "no number"), (4, "no degree of reduction")]
    assert b["lower"] is None and b["upper"] is None


def test_a_second_window_compound_is_not_counted():
    b = chemistry.electron_balance([10, 4], [12, 0], [24, 24], [True, False])
    assert b["consumed_e_mM"] == 240 and b["not_counted"] == [(1, "measured over the second time window")]


def test_the_range_takes_the_fewest_electrons_out_over_the_most_in_and_the_reverse():
    bounds = {"consumed": ([8, 0], [12, 0.2]), "produced": ([0, 10], [0.2, 20])}
    b = chemistry.electron_balance([10, 0], [0, 15], [24, 8], [True, True], bounds)
    assert b["lower"] == pytest.approx(10 * 8 / (12 * 24 + 0.2 * 8))
    assert b["upper"] == pytest.approx((0.2 * 24 + 20 * 8) / (8 * 24))
    # a number without bounds (one replicate) gives no range
    bounds["produced"] = ([0, None], [0.2, None])
    assert chemistry.electron_balance([10, 0], [0, 15], [24, 8], [True, True], bounds)["lower"] is None


def test_crm_mode_reads_the_chemistry_and_sends_it_with_the_balance(client):
    r = run(client, report_rates=True)
    p = matrix.crm_payload(r)
    chem = dict(zip(p["resources"], p["resource_chemistry"], strict=True))
    assert chem["glucose"]["degree_of_reduction"] == 24 and chem["acetate"]["degree_of_reduction"] == 8
    assert chem["butyrate"]["degree_of_reduction"] == 20 and chem["formate"]["degree_of_reduction"] == 2
    assert p["chemistry_source"]["source"].startswith("ChEBI") and p["chemistry_source"]["problem"] is None
    gammas = [c["degree_of_reduction"] for c in p["resource_chemistry"]]
    for i, b in enumerate(p["electron_balance"]):
        consumed = sum((v or 0) * g for v, g in zip(p["consumed"][i], gammas, strict=True))
        produced = sum((v or 0) * g for v, g in zip(p["produced"][i], gammas, strict=True))
        assert b["consumed_e_mM"] == pytest.approx(consumed) and b["produced_e_mM"] == pytest.approx(produced)
        if consumed > 0:
            assert b["share"] == pytest.approx(produced / consumed)
    assert any(b["share"] is not None for b in p["electron_balance"])
    z = zipfile.ZipFile(io.BytesIO(matrix.crm_package(r)))
    rows = list(csv.DictReader(io.StringIO(z.read("chemistry.csv").decode())))
    assert {row["resource"]: row["degree_of_reduction"] for row in rows}["glucose"] == "24"
    balance = list(csv.DictReader(io.StringIO(z.read("electron_balance.csv").decode())))
    assert balance and {row["phase"] for row in balance} == {"exponential"}
    assert "Electron balance:" in p["readme"] and "degree of reduction" in p["readme"]


def test_without_crm_mode_chebi_is_not_asked():
    class Counting(FakeClient):
        asked = 0

        def chebi_compound(self, chebi_id):
            Counting.asked += 1
            return CHEBI_RECORDS.get(int(chebi_id))
    r = run(Counting())
    assert Counting.asked == 0 and r["chemistry"] == {}


def test_chebi_out_of_reach_is_a_warning_and_the_parameters_still_come():
    class Down(FakeClient):
        def chebi_compound(self, chebi_id):
            raise OSError("ChEBI timed out")
    r = run(Down(), report_rates=True)
    assert any("ChEBI could not be reached" in w for w in r["warnings"])
    p = matrix.crm_payload(r)
    assert all(c["degree_of_reduction"] is None for c in p["resource_chemistry"])
    assert all(b is None for b in p["electron_balance"]) and p["consumed"]
