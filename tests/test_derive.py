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
    assert r["cells"][(B, BUT, "exponential")]["cautions"] == ["pair_decided"]          # two replicates
    assert r["cells"][(B, BUT, "exponential")]["n"] == 2


def test_the_stationary_phase_starts_at_the_boundary(client):
    r = run(client, phase="stationary", judge_confidence=False)
    # A's glucose from 12 h to 24 h: 0 and -1, mean -0.5, beyond the 0.2 mM limit (judged by the mean alone;
    # by default its spread reaches across the limit, test_matrix.py)
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


def test_a_filled_second_box_uses_only_the_data_matching_it(client):
    # Karoline, 2026-10-05: "by default, when something is entered in the 2nd field, only data matching what
    # was entered are shown"
    r = run(client, conditions="mMCB")
    assert r["value_rule"]["rule"] == "selected"
    assert r["cells"][(C, GLC, "exponential")]["mean"] == pytest.approx(-5.0)
    assert r["presence"] == {} and (A, GLC, "exponential") not in r["cells"]
    assert {s for e in r["network"].edges for s in e.study_ids} == {"SMGDB00000002"}


def test_outside_evidence_is_an_advanced_option(client):
    # "... but in advanced settings, we can switch on showing supporting evidence from other studies"
    r = run(client, conditions="mMCB", outside_evidence=True)
    assert r["cells"][(C, GLC, "exponential")]["mean"] == pytest.approx(-5.0)
    # WC gives presence only: A's glucose uptake is a presence arc
    assert r["presence"][(A, GLC, "exponential")]["consumed"][0]["medium"].startswith("Wilkins")


def test_an_empty_second_box_considers_all_data(client):
    # "By default, when the 2nd field is left empty, always all data are considered"
    r = run(client)
    assert set(r["studies"]) == {"SMGDB00000001", "SMGDB00000002"} and r["presence"]


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
        "EMGDB000000005", "SMGDB00000003", 2, "Wilkins-Chalgren Anaerobe Broth", ["b3", "b3b"])
    series = dict(conftest.SERIES)
    series["b3"] = {"growth": ("od", conftest.SLOW), "glucose": [(0, 10), (8, 9), (16, 6.1), (24, 6.1)]}
    series["b3b"] = {"growth": ("od", conftest.SLOW), "glucose": [(0, 10), (8, 9), (16, 5.9), (24, 5.9)]}
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


def test_no_phase_from_one_experiment_is_no_caution_on_a_value_of_another():
    # a review: with Both, E2 reached no end of growth (no_phase on its stationary rows), E1 has three stationary
    # values; the cell's value carried no_phase, which is no caution, and its arc failed validation
    cultures = [_culture("E1", "S1", f"r{i}", {"m": [(0, 5), (10, 3)]}) for i in range(3)] + \
               [_culture("E2", "S2", f"r{i}", {"m": [(0, 5), (10, 3)]}) for i in range(3)]
    rows = [{"culture": i, "taxon": "ncbi:9", "metabolite": "m", "phase": "stationary", "change": -2.0 - i / 10,
             "start": 0, "end": 10, "initial": 5, "cautions": []} for i in range(3)] + \
           [{"culture": i, "taxon": "ncbi:9", "metabolite": "m", "phase": "stationary", "change": None,
             "start": None, "end": None, "initial": 5, "cautions": ["no_phase"]} for i in range(3, 6)]
    cell = derive.pool(rows, cultures)[("ncbi:9", "m", "stationary")]
    assert cell["direction"] == "consumed" and "no_phase" not in cell["cautions"]
    # and it rests on E1 alone, naming E2
    assert cell["experiments"] == ["E1"] and cell["studies"] == ["S1"]
    assert any("no end of growth" in n and "E2" in n for n in cell["notes"])
    # a conflict in the same cell, and in a later one, still pools (a review: the note's local name replaced the
    # experiment names the conflict notes use)
    cultures += [_culture("E3", "S3", f"r{i}", {"m": [(0, 5), (10, 7)]}) for i in range(3)]
    rows += [{**r, "culture": r["culture"] + 6, "change": 2.0} for r in rows[:3]]
    rows += [{**r, "metabolite": "z"} for r in rows[:3]] + [{**r, "metabolite": "z"} for r in rows[6:9]]
    pooled = derive.pool(rows, cultures)
    for met in ("m", "z"):
        assert "conflict" in pooled[("ncbi:9", met, "stationary")]["cautions"], met


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
    # glucose starts at 10 in all five WC replicates
    assert r["initial"][GLC]["mean"] == pytest.approx(10.0) and r["initial"][GLC]["n"] == 5


def test_every_tested_arc_carries_its_q_value_and_identical_replicates_none(client):
    net = run(client, phase="stationary", judge_confidence=False)["network"]
    arc = next(e for e in net.edges if e.taxon == A and e.metabolite == GLC)
    # A's stationary glucose (0 and -1) is tested, and corrected within the family of the per-study arcs
    assert arc.q_value is not None and arc.p_value is not None
    from foodnet import stats
    # identical values are no certainty: no test
    assert stats.paired([0.21, 0.21], [0.0, 0.0])["p"] is None


def test_the_length_of_the_exponential_phase_reaches_the_arcs(client):
    # A grows from its first sample at 0 h until 12 h, B until 16 h (conftest.py)
    net = run(client)["network"]
    by = {(e.taxon, e.metabolite): e.exponential_h for e in net.edges if e.evidence == "measured"}
    assert by[(A, GLC)] == 12.0 and by[(B, GLC)] == 16.0
    # with a time window the length is still reported: it describes the culture
    net = run(client, window_start=0.0, window_end=6.0)["network"]
    assert {e.exponential_h for e in net.edges if e.taxon == A} == {12.0}


def test_a_second_box_that_matches_nothing_says_so_and_names_the_media(client):
    r = run(client, conditions="Db-MM")
    assert r["value_cultures"] == 0 and r["network"].edges == []
    warning = next(w for w in r["warnings"] if "Nothing in the second box" in w)
    assert "Wilkins-Chalgren Anaerobe Broth (WC)" in warning and "mMCB" in warning


def test_study_ids_in_the_second_box_limit_what_is_read(client):
    # Karoline, 2026-10-05: "when I gave a list of studies, the results also included studies that were not
    # in my list. This is not desired behavior."
    r = run(client, conditions="SMGDB00000002")
    assert r["studies"] == ["SMGDB00000002"]
    assert {s for e in r["network"].edges for s in e.study_ids} == {"SMGDB00000002"}
    assert r["presence"] == {}                          # nothing from study 1, not even as presence
    assert r["cells"][(B, FOR, "exponential")]["mean"] == pytest.approx(5.0)


def test_an_experiment_id_limits_the_search_to_that_experiment(client):
    r = run(client, conditions="EMGDB000000004")
    assert r["studies"] == ["SMGDB00000002"]
    assert {e.taxon for e in r["network"].edges} == {C}   # E3 (Beta beta in the same study) is left out
    assert any("limits the search" in reason for _, reason in r["skipped"])


def test_study_ids_with_outside_evidence_read_the_other_studies_too(client):
    r = run(client, conditions="SMGDB00000002", outside_evidence=True)
    assert set(r["studies"]) == {"SMGDB00000001", "SMGDB00000002"}
    assert (A, GLC, "exponential") in r["presence"]


def test_ids_set_the_scope_and_the_majority_rule_runs_inside_it(client):
    # Karoline, 2026-10-05, reproducing a hand-checked reference from four study ids: the ids say which data,
    # and within them the majority medium gives the values while the other media give presence
    r = run(client, conditions="SMGDB00000001\nSMGDB00000002")
    assert r["value_rule"]["rule"] == "majority_in_scope"
    assert r["value_rule"]["keys"] == ["wilkins chalgren anaerobe broth"]
    assert r["cells"][(A, AC, "exponential")]["mean"] == pytest.approx(4.0)
    assert "produced" in r["presence"][(B, FOR, "exponential")]          # mMCB, inside the scope: presence


def test_a_medium_name_chooses_the_value_medium_inside_the_ids(client):
    # the four studies and "Wilkins-Chalgren": values from WC only, the studies' other media as presence
    r = run(client, conditions="SMGDB00000001\nSMGDB00000002\nmMCB")
    assert r["value_rule"]["rule"] == "selected" and r["value_rule"]["media"] == ["mMCB"]
    assert r["cells"][(C, GLC, "exponential")]["mean"] == pytest.approx(-5.0)
    assert "consumed" in r["presence"][(A, GLC, "exponential")]


def test_a_strain_is_shown_by_its_current_name():
    # the species list says taxon 2's current name is another than the one its studies record
    from foodnet.search import run_query
    from foodnet.taxonomy import SpeciesIndex
    index = SpeciesIndex({"alpha alpha": {1: "Alpha alpha A1"}, "beta beta": {2: "Beta beta B1"}},
                         current={2: "Betanova beta B1"}, studies=["SMGDB00000001", "SMGDB00000002"],
                         where={1: {"SMGDB00000001"}, 2: {"SMGDB00000001", "SMGDB00000002"}})
    r = run_query(FakeClient(), ["Alpha alpha", "Beta beta"], {}, index)
    assert {n.name for n in r["taxa_nodes"]} == {"Alpha alpha A1", "Betanova beta B1"}


def test_a_second_window_measures_the_metabolites_it_names_over_its_own_interval(client):
    # Karoline, 2026-10-05: "2 time windows, with the option to add metabolites by name"; an empty end is
    # the last sample. A's glucose over 0 h to the end (24 h): 10 to 2 and 10 to 1, mean -8.5; acetate
    # keeps the exponential phase (+4)
    r = run(client, second_window_metabolites="Glucose")
    assert r["cells"][(A, GLC, "window")]["mean"] == pytest.approx(-8.5)
    assert (A, GLC, "exponential") not in r["cells"]
    assert r["cells"][(A, AC, "exponential")]["mean"] == pytest.approx(4.0)
    assert any("second time window, 0 h to the last sample" in w for w in r["warnings"])


def test_a_second_window_with_an_end_and_a_name_given_as_the_recorded_acid(client):
    # "acetic acid" is how the study recorded acetate; 0 to 6 h: 0 to 1 in both replicates
    r = run(client, second_window_metabolites="acetic acid", second_window_start=0.0, second_window_end=6.0)
    assert r["cells"][(A, AC, "window")]["mean"] == pytest.approx(1.0)


def test_a_second_window_name_that_matches_nothing_is_reported(client):
    r = run(client, second_window_metabolites="trehalose")
    assert any("trehalose, which no culture" in w for w in r["warnings"])


def _with_mucin_variant():
    """Study one gains Alpha alpha in WC plus mucin beads, where it takes up 2 mM of glucose more."""
    import conftest
    conftest.STUDIES["SMGDB00000001"]["experiments"].append({"id": "EMGDB000000006", "name": "A_MUCIN"})
    exp = conftest._experiment("EMGDB000000006", "SMGDB00000001", 1, "Wilkins-Chalgren Anaerobe Broth", ["a9", "a9b"])
    exp["description"] = "Alpha with WC plus mucin beads"
    conftest.EXPERIMENTS["EMGDB000000006"] = exp
    series = dict(conftest.SERIES)
    series["a9"] = {"growth": ("fc", conftest.FAST), "glucose": [(0, 10), (6, 6), (12, 0.1), (24, 0)]}
    series["a9b"] = {"growth": ("fc", conftest.FAST), "glucose": [(0, 10), (6, 6), (12, -0.1), (24, 0)]}
    return FakeClient(series)


def _drop_mucin_variant():
    import conftest
    conftest.STUDIES["SMGDB00000001"]["experiments"].pop()
    del conftest.EXPERIMENTS["EMGDB000000006"]


def test_a_medium_altered_in_its_description_gives_presence_not_values():
    # Karoline, 2026-10-06: "check descriptions that suggest something altered the medium ... and treat it as
    # another medium". "Wilkins" matches both; plain WC gives the values (A's glucose -8 as before)
    client = _with_mucin_variant()
    try:
        r = run(client, conditions="Wilkins")
        assert r["cells"][(A, GLC, "exponential")]["mean"] == pytest.approx(-8.0)
        assert r["value_rule"]["also_matched"] == ["Wilkins-Chalgren Anaerobe Broth (+mucin)"]
        assert any("also matched" in w for w in r["warnings"])
        # switched off, the mucin variant pools into the values: experiments at -8 (three replicates) and -10
        # (two), each counted once, mean -9
        r = run(client, conditions="Wilkins", strict_media=False)
        assert r["cells"][(A, GLC, "exponential")]["mean"] == pytest.approx(-9.0)
        assert r["cells"][(A, GLC, "exponential")]["n"] == 5
        assert r["cells"][(A, GLC, "exponential")]["n_experiments"] == 2
    finally:
        _drop_mucin_variant()


def test_experiments_can_be_excluded_by_id():
    client = _with_mucin_variant()
    try:
        r = run(client, conditions="Wilkins", strict_media=False, exclude_experiments="EMGDB000000006")
        assert r["cells"][(A, GLC, "exponential")]["mean"] == pytest.approx(-8.0)
        assert any("Exclude these experiments" in reason for _, reason in r["skipped"])
    finally:
        _drop_mucin_variant()
