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
    # 10 mM glucose (24) in; 15 mM acetate (8) and 5 mM butyrate (20) out; a resource without a number is named,
    # and so is one the medium holds
    b = chemistry.electron_balance([10, 0, 0, None], [0, 15, 5, None], [24, 8, 20, 8], [True] * 4,
                                   in_medium=[True, True, False, True])
    assert b["consumed_electrons_mM"] == 240 and b["produced_electrons_mM"] == 220
    assert b["share"] == pytest.approx(220 / 240) and b["withheld"] is None
    assert b["not_counted"] == [(3, "no number")] and b["incomplete"] == [3]


def test_a_resource_without_a_degree_of_reduction_withholds_the_share():
    # Karoline, 2026-10-08, "Own gaps + name the rest": left out, a substrate would make the share too high;
    # and "Withhold for products too": a product left out makes it too low
    b = chemistry.electron_balance([10, 3], [12, 0], [8, None], [True, True])
    assert b["share"] is None and "no degree of reduction" in b["withheld"]
    assert (1, "no degree of reduction") in b["not_counted"]
    b = chemistry.electron_balance([10, 0], [12, 3], [8, None], [True, True])
    assert b["share"] is None and b["withheld"]
    # a resource without a degree of reduction that it neither took up nor made changes nothing
    b = chemistry.electron_balance([10, 0], [12, 0], [8, None], [True, True])
    assert b["share"] == pytest.approx(96 / 80) and b["withheld"] is None


def test_uptake_within_the_detection_limits_withholds_the_share():
    # Karoline, 2026-10-08, "Withhold, say why": 0.25 mM glucose (6 electrons) taken up, while glucose and acetate
    # could each hide 0.2 mM (0.2 x 24 + 0.2 x 8 = 6.4)... and a little more uptake is a share
    b = chemistry.electron_balance([0.25, 0], [0, 2], [24, 8], [True, True], limits=[0.2, 0.2])
    assert b["share"] is None and "within the detection limits" in b["withheld"]
    b = chemistry.electron_balance([1, 0], [0, 2], [24, 8], [True, True], limits=[0.2, 0.2])
    assert b["share"] == pytest.approx(16 / 24)
    # a culture whose own uptake is that small is left out of the range, and counted
    b = chemistry.electron_balance([1, 0], [0, 2], [24, 8], [True, True],
                                   [[-1, 2], [-1.2, 2.1], [-0.1, 3]], limits=[0.2, 0.2])
    assert (b["cultures"], b["cultures_left_out"]) == (2, 1)


def test_a_withheld_small_uptake_says_what_the_share_is_at_least():
    # Karoline, 2026-10-08, "Report 'at least X'": out over in plus what the limits hide (6 + 6.4)
    b = chemistry.electron_balance([0.25, 0], [0, 2], [24, 8], [True, True], limits=[0.2, 0.2], starts=[1, 1])
    assert b["share"] is None and b["share_upper"] is None
    assert b["share_at_least"] == pytest.approx(16 / (6 + 6.4)) and "at least 1.29" in b["withheld"]
    assert b["share_lower"] is None
    # a resource without a number hides what the medium held when the phase began, not just its limit (a review:
    # "at least 2.01" where the cultures allowed 0.46)
    b = chemistry.electron_balance([0.25, 0, None], [0, 2, None], [24, 8, 20], [True] * 3,
                                   in_medium=[True, True, True], limits=[0.2, 0.2, 0.2], starts=[1, 1, 3])
    assert b["share_at_least"] == pytest.approx(16 / (6 + 6.4 + 60))
    # and the lowest culture's floor, not only the pooled amounts' (a review)
    b = chemistry.electron_balance([0.25, 0], [0, 2], [24, 8], [True, True], [[-0.25, 0.5], [-0.25, 3.5]],
                                   limits=[0.2, 0.2], starts=[1, 1])
    assert b["share_at_least"] == pytest.approx(4 / (6 + 6.4)) and b["cultures_left_out"] == 2
    # a culture that measured the resource the pooled values lack is bounded by its own change, not by the
    # pooled phase-start amount (a review: a culture can start with more than the mean)
    b = chemistry.electron_balance([0.25, 0, None], [0, 2, None], [24, 8, 20], [True] * 3,
                                   [[-0.25, 2, -0.5], [-0.25, 2, None]],
                                   in_medium=[True, True, True], limits=[0.2, 0.2, 0.2], starts=[1, 1, 3])
    pooled = 16 / (6 + 6.4 + 60)
    own = 16 / (6 + 10 + 0.2 * (24 + 8 + 20))
    assert b["share_at_least"] == pytest.approx(min(pooled, own))


def test_the_floor_charges_what_an_unmeasured_resource_held_and_names_what_nobody_measured():
    # a review: a resource without a number, held below the limit, could still hide what it held
    b = chemistry.electron_balance([0.25, 0, None, None], [0, 2, None, None], [24, 8, 26, 24], [True] * 4,
                                   in_medium=[True, True, False, False], limits=[0.2] * 4,
                                   starts=[1, 1, 0.05, None], names=["glucose", "acetate", "isovalerate", "fructose"])
    assert b["share_at_least"] == pytest.approx(16 / (6 + 6.4 + 0.05 * 26))
    assert "taking fructose, which nobody measured in this medium, as absent" in b["withheld"]
    # and a culture's floor charges it too (a review)
    b = chemistry.electron_balance([0.25, 0, None], [0, 2, None], [24, 8, 26], [True] * 3, [[-0.25, 0.5, None]],
                                   in_medium=[True, True, False], limits=[0.2] * 3, starts=[1, 1, 0.05])
    assert b["share_at_least"] == pytest.approx(min(16 / (6 + 6.4 + 1.3), 4 / (6 + 6.4 + 1.3)))


def test_the_floor_says_it_leaves_out_the_second_window():
    # a review: a second-window compound, measured over its own window, can be taken up within the phase
    b = chemistry.electron_balance([0.25, 0, 3], [0, 2, 0], [24, 8, 48], [True, True, False], limits=[0.2] * 3,
                                   starts=[1, 1, 1])
    assert "leaving out the second time window's compounds" in b["withheld"]


def test_a_share_outside_its_cultures_range_is_withheld():
    # Karoline, 2026-10-08, "Withhold when outside": a culture that took up almost nothing, and made much, drives
    # the pooled share but is left out of the range
    cultures = [[-10, 10], [-10, 11], [-0.1, 60]]
    b = chemistry.electron_balance([6.7, 0], [0, 27], [24, 8], [True, True], cultures, limits=[0.2, 0.2])
    assert b["cultures_left_out"] == 1 and b["share"] is None and "outside its cultures" in b["withheld"]
    assert (b["share_lower"], b["share_upper"]) == (pytest.approx(80 / 240), pytest.approx(88 / 240))


def test_a_second_window_compound_is_not_counted():
    b = chemistry.electron_balance([10, 4], [12, 0], [24, 24], [True, False])
    assert b["consumed_electrons_mM"] == 240 and b["not_counted"] == [(1, "measured over the second time window")]


def test_the_range_is_the_lowest_and_highest_cultures_own_share():
    # Karoline, 2026-10-08, "Per-culture shares": each culture is a closed balance; a culture that missed a
    # resource counted is left out of the range
    cultures = [[-10, 14], [-12, 15], [-8, 13], [None, 20]]
    b = chemistry.electron_balance([10, 0], [0, 14], [24, 8], [True, True], cultures)
    assert b["cultures"] == 3
    assert b["share_lower"] == pytest.approx(15 * 8 / (12 * 24))
    assert b["share_upper"] == pytest.approx(13 * 8 / (8 * 24))
    # one culture is no range
    one = chemistry.electron_balance([10, 0], [0, 14], [24, 8], [True, True], cultures[:1])
    assert one["share_lower"] is None and one["cultures"] == 1
    # an electron acceptor counted gives no range
    acc = chemistry.electron_balance([10, 1], [0, 0], [24, -8], [True, True], cultures[:3])
    assert acc["share_lower"] is None


def _chebi_answer(monkeypatch, answers):
    """fetch_compound against a stand-in for urllib: `answers` is a list of what each request gets (a dict
    for JSON, an exception to raise)."""
    import json as _json
    calls = []

    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def urlopen(request, timeout):
        calls.append(request.full_url)
        answer = answers.pop(0)
        if isinstance(answer, BaseException):
            raise answer
        return Response(_json.dumps(answer).encode())
    monkeypatch.setattr(chemistry.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(chemistry.time, "sleep", lambda s: None)
    return calls


ACETATE = {"id": 30089, "name": "acetate", "ascii_name": "acetate",
           "chemical_data": {"formula": "C2H3O2", "charge": -1, "mass": "59.044"}}


def test_a_dropped_answer_is_tried_again(monkeypatch):
    import http.client
    calls = _chebi_answer(monkeypatch, [http.client.IncompleteRead(b"{"), ACETATE])
    assert chemistry.fetch_compound("30089")["formula"] == "C2H3O2" and len(calls) == 2


def test_an_unreachable_chebi_is_left_alone_for_the_rest_of_the_search_only(monkeypatch):
    import http.client
    calls = _chebi_answer(monkeypatch, [OSError("timed out"), http.client.IncompleteRead(b""), ACETATE])
    client = chemistry.Memo()                         # one client, as the page keeps for an hour
    found, problem = chemistry.resolve({"a": "30089", "g": "17234"}, client)
    assert found["g"]["note"] == "ChEBI could not be reached" and "2 of 2" in problem
    assert len(calls) == chemistry.ATTEMPTS           # the second resource did not wait again
    # the next search with the same client asks again (a review: a rerun, as the warning advises, failed)
    found, problem = chemistry.resolve({"a": "30089"}, client)
    assert found["a"]["formula"] == "C2H3O2" and problem is None
    # and one that fails does not blank what the client read before (a review)
    _chebi_answer(monkeypatch, [OSError("timed out"), OSError("timed out")])
    found, problem = chemistry.resolve({"g": "17234", "a": "30089"}, client)
    assert found["a"]["formula"] == "C2H3O2" and "1 of 2" in problem


def test_a_changed_or_moved_api_is_unreachable_not_a_compound_without_formula(monkeypatch):
    _chebi_answer(monkeypatch, [{"id": 30089, "name": "acetate", "structure": {}}])
    with pytest.raises(chemistry.ChebiUnreachable, match="lacks the fields"):
        chemistry.fetch_compound(30089)
    import urllib.error

    def html_404():
        return urllib.error.HTTPError("u", 404, "Not Found", {}, io.BytesIO(b"<html>Not Found</html>"))
    _chebi_answer(monkeypatch, [html_404(), html_404()])
    with pytest.raises(chemistry.ChebiUnreachable, match="moved"):
        chemistry.fetch_compound(30089)
    _chebi_answer(monkeypatch, [urllib.error.HTTPError("u", 404, "Not Found", {},
                                                       io.BytesIO(b'{"detail": "No Compound matches"}'))])
    assert chemistry.fetch_compound(999999999) is None

    class Cut(io.BytesIO):
        def read(self, *a):
            import http.client
            raise http.client.IncompleteRead(b"{")
    _chebi_answer(monkeypatch, [urllib.error.HTTPError("u", 404, "Not Found", {}, Cut()),
                                urllib.error.HTTPError("u", 404, "Not Found", {}, Cut())])
    with pytest.raises(chemistry.ChebiUnreachable):
        chemistry.fetch_compound(30089)


def test_a_prefixed_or_secondary_id_is_read_as_its_compound(monkeypatch):
    calls = _chebi_answer(monkeypatch, [ACETATE])
    assert chemistry.fetch_compound("CHEBI:30089")["id"] == 30089 and calls[0].endswith("/30089/")
    assert compounds.canonical("CHEBI:15366", "acetic acid")["id"] == "chebi:30089"
    # a secondary id of fructose answers with fructose's record, whose class has a named form
    fetch, _ = _fake({5172: {"id": 28757, "name": "fructose", "formula": "", "charge": None, "mass": None},
                      15824: {"id": 15824, "name": "D-fructose", "formula": "C6H12O6", "charge": 0, "mass": None}})
    assert chemistry._one("5172", fetch)["degree_of_reduction"] == 24


def test_a_salt_says_it_is_one():
    assert "salt" in chemistry.degree_of_reduction("C3H3O3.Na", 0)[2]


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
        assert b["consumed_electrons_mM"] == pytest.approx(consumed)
        assert b["produced_electrons_mM"] == pytest.approx(produced)
        if consumed > 0:
            assert b["share"] == pytest.approx(produced / consumed)
            # every culture's own share lies around the pooled one
            if b["share_lower"] is not None:
                assert b["share_lower"] <= b["share_upper"] and b["cultures"] >= 2
    assert any(b["share"] is not None for b in p["electron_balance"])
    assert any(b["share_lower"] is not None for b in p["electron_balance"])
    assert len(p["culture_changes"]) == len(p["taxa"])
    z = zipfile.ZipFile(io.BytesIO(matrix.crm_package(r)))
    rows = list(csv.DictReader(io.StringIO(z.read("chemistry.csv").decode())))
    assert {row["resource"]: row["degree_of_reduction"] for row in rows}["glucose"] == "24"
    balance = list(csv.DictReader(io.StringIO(z.read("electron_balance.csv").decode())))
    assert balance and {row["phase"] for row in balance} == {"exponential"}
    assert list(balance[0])[2:] == list(matrix.BALANCE_COLUMNS)
    assert "Electron balance:" in p["readme"] and "no expected side of 1" in p["readme"]


def test_both_phases_carry_a_balance_each(client):
    p = matrix.crm_payload(run(client, report_rates=True, phase="both"))
    assert p["other_phases"]["stationary"]["electron_balance"] is not None
    assert len(p["other_phases"]["stationary"]["culture_changes"]) == len(p["taxa"])


def test_booleans_have_no_balance(client):
    p = matrix.crm_payload(run(client, report_rates=True, booleans=True))
    assert all(b is None for b in p["electron_balance"])


def test_the_report_gives_the_share_the_cultures_and_what_is_missing(client):
    from foodnet.report import report_text
    text = report_text(run(client, report_rates=True))
    assert "Electron balance" in text and "1 is no target" in text
    assert "below 1 expected" not in text


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
    assert any("ChEBI could not be reached" in w and "withheld" in w for w in r["warnings"])
    p = matrix.crm_payload(r)
    assert all(c["degree_of_reduction"] is None for c in p["resource_chemistry"])
    assert all(b is None for b in p["electron_balance"]) and p["consumed"]


def test_a_cultures_change_inside_the_limit_is_zero_as_in_the_matrices():
    # a review: acetate -0.155 mM, below the 0.2 mM limit, was one culture's only substrate (share 39.7)
    from foodnet import crm

    class C:                                         # the label is all culture_changes reads
        label = "c1"
    rows = [{"culture": 0, "taxon": "t", "metabolite": "m", "phase": "exponential", "change": -0.15,
             "cautions": []},
            {"culture": 0, "taxon": "t", "metabolite": "n", "phase": "exponential", "change": -3.0, "cautions": []}]
    assert crm.culture_changes(rows, [C()], 0.2) == {("t", "exponential"): {"c1": {"m": 0.0, "n": -3.0}}}
    assert crm.culture_changes(rows, [C()], 0.2, {"m": 0.1})[("t", "exponential")]["c1"]["m"] == -0.15


def test_the_medium_is_judged_at_the_phases_start(client):
    # Karoline, 2026-10-08, "At the phase's start": a substrate used up before the phase is not missing from it
    r = run(client, report_rates=True, phase="both")
    taxa = matrix.taxa_rows(r["network"], r)
    cols = matrix.crm_columns(r["network"], r, "stationary")
    starts = matrix.phase_start(r, taxa, cols)
    glucose = [m.name for m, _, _ in cols].index("glucose")
    alpha = [t.name for t in taxa].index("Alpha alpha A1")
    # Alpha's cultures held about 2 mM glucose at the end of their exponential phase, not the medium's 10
    assert 1 < starts[alpha][glucose] < 3


def test_one_resource_out_of_reach_withholds_the_share_of_whoever_consumed_it():
    # a review: a transient failure for glucose alone turned shares of 0.7-1.1 into 2.4-4.3
    class Glucose(FakeClient):
        def chebi_compound(self, chebi_id):
            if int(chebi_id) == 17234:
                raise OSError("ChEBI returned HTTP 503")
            return CHEBI_RECORDS.get(int(chebi_id))
    p = matrix.crm_payload(run(Glucose(), report_rates=True))
    glucose = p["resources"].index("glucose")
    assert any((row[glucose] or 0) > 0 for row in p["consumed"])
    for i, b in enumerate(p["electron_balance"]):
        if (p["consumed"][i][glucose] or 0) > 0:
            assert b["share"] is None and b["withheld"]


FIXTURE = "r/tests/testthat/crm_synthetic.json"


def _steady_series() -> dict:
    """The synthetic world with Alpha still taking up glucose and making acetate after its exponential phase, the
    same way in each culture, so the stationary phase has a share and a range."""
    import copy

    from conftest import SERIES
    series = copy.deepcopy(SERIES)
    for rid, (glucose, acetate) in {"a1": (0.9, 5.0), "a2": (1.0, 4.9), "a3": (1.1, 4.8)}.items():
        series[rid]["glucose"][-1] = (24, glucose)
        series[rid]["acetic acid"][-1] = (24, acetate)
    return series


class _NoButyrate(FakeClient):
    """ChEBI without butyrate's formula: Beta, which makes butyrate, has its share withheld."""

    def chebi_compound(self, chebi_id):
        record = CHEBI_RECORDS.get(int(chebi_id))
        return {**record, "formula": ""} if record and record["name"] == "butyrate" else record


def synthetic_payload() -> dict:
    """The CRM parameters of the synthetic world with both phases, formate in a second time window and Beta's
    share withheld, without what changes from run to run: the R package's tests recompute the electron balance
    from it and compare with Python's (a review: nothing tied the two together)."""
    import json
    p = matrix.crm_payload(run(_NoButyrate(_steady_series()), report_rates=True, phase="both",
                               second_window_metabolites="formate",
                               second_window_start=0.0, second_window_end=24.0))
    for key in ("derived_at", "tool_version", "readme", "settings"):
        p.pop(key, None)
    p["chemistry_source"]["retrieved_at"] = ""
    p["caveats"]["warnings"] = []
    return json.loads(json.dumps(p))


def test_the_r_package_fixture_is_what_foodnet_writes():
    import json
    import pathlib
    path = pathlib.Path(__file__).resolve().parents[1] / FIXTURE
    assert json.loads(path.read_text()) == synthetic_payload(), \
        f"regenerate {FIXTURE}: python -c 'import sys; sys.path.insert(0, \"tests\"); import json, test_chemistry " \
        f"as t; open(\"{FIXTURE}\", \"w\").write(json.dumps(t.synthetic_payload(), indent=1))'"


def test_a_second_search_on_the_same_client_reads_chebi_again():
    # a review: the page keeps one client for up to an hour, and one outage blanked every later search
    class Flaky(FakeClient):
        down = True

        def __init__(self):
            super().__init__()
            self._chebi = chemistry.Memo(self._read)

        def _read(self, chebi_id):
            if Flaky.down:
                raise chemistry.ChebiUnreachable("timed out")
            return CHEBI_RECORDS.get(int(chebi_id))

        def chebi_compound(self, chebi_id):
            return self._chebi(chebi_id)
    client = Flaky()
    first = run(client, report_rates=True)
    assert any("ChEBI could not be reached" in w for w in first["warnings"])
    Flaky.down = False
    second = run(client, report_rates=True)
    assert not any("ChEBI" in w for w in second["warnings"])
    assert all(c["degree_of_reduction"] is not None for c in second["chemistry"].values())
