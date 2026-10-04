"""One compound under two ChEBI ids, and units to mM (foodnet.compounds)."""
import pytest

from foodnet import compounds


def test_an_acid_joins_its_conjugate_base():
    # study 9 records "acetic acid" (CHEBI:15366), study 7 "acetate" (CHEBI:30089): one node
    acid = compounds.canonical(15366, "acetic acid")
    base = compounds.canonical(30089, "acetate")
    assert acid["id"] == base["id"] == "chebi:30089"
    assert acid["name"] == "acetate" and "conjugate acid" in acid["joined"] and base["joined"] == ""


def test_lactic_acid_joins_lactate():
    assert compounds.canonical(422, "(S)-lactic acid")["id"] == "chebi:24996"


def test_a_compound_without_a_chebi_id_is_keyed_by_name():
    assert compounds.canonical(None, "Mucin  Sugars")["id"] == "metabolite:mucin sugars"


def test_units_to_millimolar():
    assert compounds.to_mm("mM") == (1.0, None)
    assert compounds.to_mm("µM") == (1e-3, None)
    # 1 g/L of glucose (180.156 g/mol) is 5.5508 mM
    factor, why = compounds.to_mm("g/L", "17234")
    assert why is None and factor == pytest.approx(1000 / 180.156)


def test_a_unit_that_is_not_a_concentration_is_refused():
    factor, why = compounds.to_mm("AUC")
    assert factor is None and "not a concentration" in why
    factor, why = compounds.to_mm("g/L", "99999999")
    assert factor is None and "molar mass" in why
