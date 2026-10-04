"""From monoculture series to values, presence and arcs, on the synthetic world of conftest.py.

Every expected number is worked out in conftest.py's docstring.
"""
import pytest
from conftest import FakeClient, run

from foodnet import derive
from foodnet.reading import Culture

A, B, C = "ncbi:1", "ncbi:2", "ncbi:3"
GLC, AC, BUT, FOR = "chebi:17234", "chebi:30089", "chebi:17968", "chebi:15740"


def test_exponential_phase_values_are_the_mean_net_change(client):
    r = run(client)
    assert r["cells"][(A, GLC, "exponential")]["mean"] == pytest.approx(-8.0)
    assert r["cells"][(A, GLC, "exponential")]["direction"] == "consumed"
    assert r["cells"][(A, AC, "exponential")]["mean"] == pytest.approx(4.0)       # acetic acid, joined
    assert r["cells"][(B, BUT, "exponential")]["mean"] == pytest.approx(3.0)
    assert r["cells"][(B, BUT, "exponential")]["cautions"] == ["single_replicate"]


def test_the_stationary_phase_starts_at_the_boundary(client):
    r = run(client, phase="stationary")
    # A's glucose from 12 h to 24 h: 0 and -1, mean -0.5, beyond the 0.2 mM limit
    assert r["cells"][(A, GLC, "stationary")]["mean"] == pytest.approx(-0.5)
    assert r["cells"][(A, GLC, "stationary")]["direction"] == "consumed"
    # A's acetate: +1 and -1, mean 0: measured, and no change
    assert r["cells"][(A, AC, "stationary")]["direction"] is None
    assert (A, AC, "exponential") not in r["cells"]


def test_both_phases_give_both_sets_of_cells(client):
    r = run(client, phase="both")
    assert {ph for (_, _, ph) in r["cells"]} == {"exponential", "stationary"}


def test_the_detection_limit_is_a_setting(client):
    # A's stationary glucose, mean -0.5, is no change with a 1 mM limit
    r = run(client, phase="stationary", detection_limit=1.0)
    assert r["cells"][(A, GLC, "stationary")]["direction"] is None


def test_a_time_window_replaces_the_phases(client):
    # 0 to 6 h: A's glucose 10 to 8 and 10 to 7, mean -2.5
    r = run(client, window_start=0.0, window_end=6.0)
    assert r["cells"][(A, GLC, "window")]["mean"] == pytest.approx(-2.5)
    assert {ph for (_, _, ph) in r["cells"]} == {"window"}


def test_values_come_from_the_majority_medium_and_a_tie_is_reported(client):
    r = run(client)
    rule = r["value_rule"]
    assert rule["rule"] == "majority" and rule["keys"] == ["wilkins chalgren anaerobe broth"]
    assert rule["tie"] == ["wilkins chalgren anaerobe broth", "mmcb"]
    assert any("tied" in w for w in r["warnings"])


def test_other_media_give_presence_only(client):
    r = run(client)
    # B's formate rose 5 mM in mMCB: presence, no value; C's glucose fell in mMCB: presence, no value
    assert r["presence"][(B, FOR, "exponential")]["produced"][0]["medium"] == "mMCB"
    assert (B, FOR, "exponential") not in r["cells"]
    arcs = {(e.taxon, e.metabolite, e.direction): e for e in r["network"].edges}
    assert arcs[(B, FOR, "produced")].evidence == "presence_only" and arcs[(B, FOR, "produced")].amount is None
    assert arcs[(C, GLC, "consumed")].evidence == "presence_only"


def test_the_second_box_chooses_the_value_medium(client):
    r = run(client, conditions="mMCB")
    assert r["value_rule"]["rule"] == "selected"
    assert r["cells"][(C, GLC, "exponential")]["mean"] == pytest.approx(-5.0)
    # now WC gives presence only: A's glucose uptake is a presence arc
    assert r["presence"][(A, GLC, "exponential")]["consumed"][0]["medium"].startswith("Wilkins")


def test_ignoring_media_pools_every_medium(client):
    r = run(client, ignore_media=True)
    assert r["value_rule"]["rule"] == "all" and r["presence"] == {}
    assert r["cells"][(B, FOR, "exponential")]["mean"] == pytest.approx(5.0)
    assert all(e.evidence == "measured" for e in r["network"].edges)


def test_booleans_carry_no_amounts(client):
    r = run(client, booleans=True)
    assert all(e.amount is None and e.change is None for e in r["network"].edges)


def test_arcs_point_from_the_taxon_for_production_and_to_it_for_consumption(client):
    net = run(client)["network"]
    for e in net.edges:
        assert net.nodes[e.taxon].kind == "taxon" and net.nodes[e.metabolite].kind == "metabolite"
        assert (e.source == e.taxon) == (e.direction == "produced")
    assert net.validate() == []


def test_one_arc_per_study_by_default_and_one_across_studies_merged():
    # B's glucose in two WC studies: a second WC study is added for B, consuming 4 mM
    import conftest
    studies = dict(conftest.STUDIES)
    conftest.STUDIES["SMGDB00000003"] = {"id": "SMGDB00000003", "name": "three", "url": "",
                                         "experiments": [{"id": "EMGDB000000005", "name": "B_WC2"}]}
    conftest.EXPERIMENTS["EMGDB000000005"] = conftest._experiment(
        "EMGDB000000005", "SMGDB00000003", 2, "Wilkins-Chalgren Anaerobe Broth", ["b3"])
    series = dict(conftest.SERIES)
    series["b3"] = {"growth": ("od", conftest.SLOW), "glucose": [(0, 10), (8, 9), (16, 6), (24, 6)]}
    try:
        client = FakeClient(series)
        per_study = [e for e in run(client)["network"].edges if e.taxon == B and e.metabolite == GLC]
        assert sorted((e.study_ids, e.amount) for e in per_study) == [(("SMGDB00000001",), 6.0),
                                                                      (("SMGDB00000003",), 4.0)]
        merged = [e for e in run(client, merge_arcs=True)["network"].edges if e.taxon == B and e.metabolite == GLC]
        # one arc, the pooled value (-6 and -4, mean -5), resting on both studies
        assert [(e.study_ids, e.amount, e.merged_arcs) for e in merged] == [
            (("SMGDB00000001", "SMGDB00000003"), 5.0, 2)]
        assert run(client, merge_arcs=True, min_studies=2)["network"].edges[0].merged_arcs == 2
    finally:
        conftest.STUDIES.clear()
        conftest.STUDIES.update(studies)
        del conftest.EXPERIMENTS["EMGDB000000005"]


def test_merging_to_genus_gives_genus_nodes(client):
    net = run(client, merge_genera=True)["network"]
    assert {n.id for n in net.taxa()} == {"genus:Alpha", "genus:Beta", "genus:Gamma"}
    assert net.validate() == []


def _culture(experiment, study, replicate, series, medium="WC"):
    return Culture(taxon={"id": "ncbi:9", "name": "X x"}, study=study, experiment=experiment,
                   experiment_name=experiment, description="", medium=medium, medium_key=medium.lower(),
                   replicate=replicate,
                   metabolites={m: {"name": m, "chebi_id": "", "series": s} for m, s in series.items()})


def test_a_deposit_repeated_under_another_study_is_counted_once():
    # E2 holds E1's series rounded to two decimals (the R. intestinalis case of studies 2 and 7)
    g1, a1 = [(0, 7.37), (8, 6.315), (16, 0.094)], [(0, 1.708), (8, 2.839), (16, 5.965)]
    g2, a2 = [(0, 7.37), (8, 6.32), (16, 0.09)], [(0, 1.71), (8, 2.84), (16, 5.97)]
    cultures = [_culture("E1", "S7", "r1", {"glc": g1, "ac": a1}),
                _culture("E2", "S2", "q1", {"glc": g2, "ac": a2, "extra": [(0, 1), (8, 2), (16, 3)]})]
    dropped, found = derive.duplicates(cultures)
    assert dropped == {("E2", "glc"), ("E2", "ac")} or dropped == {("E1", "glc"), ("E1", "ac")}
    assert len(found) == 1 and "counted once" in found[0]


def test_different_series_are_not_duplicates():
    cultures = [_culture("E1", "S7", "r1", {"glc": [(0, 7), (8, 6), (16, 0)], "ac": [(0, 1), (8, 2), (16, 5)]}),
                _culture("E2", "S2", "q1", {"glc": [(0, 7), (8, 4), (16, 0)], "ac": [(0, 1), (8, 2), (16, 5)]})]
    assert derive.duplicates(cultures) == (set(), [])


def test_pooled_experiments_that_disagree_are_a_conflict():
    # E1 consumed 2 mM, E2 produced 2 mM: pooled mean 0, and the conflict is named
    cultures = [_culture("E1", "S1", "r1", {"m": [(0, 5), (10, 3)]}),
                _culture("E2", "S2", "r1", {"m": [(0, 5), (10, 7)]})]
    rows = [{"culture": i, "taxon": "ncbi:9", "metabolite": "m", "phase": "window", "change": d, "start": 0,
             "end": 10, "initial": 5, "cautions": []} for i, d in ((0, -2.0), (1, 2.0))]
    cell = derive.pool(rows, cultures)[("ncbi:9", "m", "window")]
    assert cell["mean"] == 0 and cell["direction"] is None and "conflict" in cell["cautions"]
    assert "E1 consumed (-2.00 mM)" in cell["notes"][0] and "E2 produced (+2.00 mM)" in cell["notes"][0]


def test_growth_rates_come_from_the_metabolite_replicates_first(client):
    from foodnet import rates
    r = run(client, report_rates=True, rate_window=3)
    expected = rates.easylinear([0, 4, 8, 12, 24], [1, 10, 100, 1000, 1000], 3)
    assert r["rates"][A]["rate"] == pytest.approx(expected)
    assert r["rates"][A]["source"] == "the replicates with metabolite data"
    # C has no culture in the value medium (WC), and a rate is never taken from another medium
    assert C in r["without_a_rate"]


def test_initial_concentrations_are_the_value_medium_first_samples(client):
    r = run(client)
    # glucose starts at 10 in all three WC replicates
    assert r["initial"][GLC]["mean"] == pytest.approx(10.0) and r["initial"][GLC]["n"] == 3


def test_a_per_study_arc_that_is_the_whole_value_carries_its_q_value(client):
    net = run(client)["network"]
    arc = next(e for e in net.edges if e.taxon == A and e.metabolite == GLC)
    # A's glucose rests on one study, so its arc is the pooled cell and has the cell's q-value
    assert arc.q_value is not None and arc.p_value is not None


def test_the_length_of_the_exponential_phase_reaches_the_arcs(client):
    # A grows from its first sample at 0 h until 12 h, B until 16 h (conftest.py)
    net = run(client)["network"]
    by = {(e.taxon, e.metabolite): e.exponential_h for e in net.edges if e.evidence == "measured"}
    assert by[(A, GLC)] == 12.0 and by[(B, GLC)] == 16.0
    # with a time window the length is still reported: it describes the culture
    net = run(client, window_start=0.0, window_end=6.0)["network"]
    assert {e.exponential_h for e in net.edges if e.taxon == A} == {12.0}
