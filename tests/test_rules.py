"""The rules a review of 0.1.0 changed (Karoline, 2026-10-06), each on a small synthetic case."""
import pytest
from conftest import run

from foodnet import derive, phase
from foodnet.reading import Culture, _growth_rank
from foodnet.search import compound_limits_of

A = "ncbi:1"
GLC = "chebi:17234"


def _culture(experiment, values, taxon="t1", medium="m", study="S1", times=(0, 4, 8, 12), growth=None):
    c = Culture(taxon={"id": taxon, "name": taxon}, study=study, experiment=experiment, experiment_name=experiment,
                description="", medium=medium, medium_key=medium, replicate=f"{experiment}{len(values)}")
    c.metabolites = {"x": {"name": "x", "chebi_id": "", "series": list(zip(times, values, strict=True))},
                     "y": {"name": "y", "chebi_id": "", "series": list(zip(times, values, strict=True))}}
    if growth:
        c.curves = growth
        c.growth = growth[0]
    return c


def _rows(changes):
    """One row per (culture index, change), for pool."""
    return [{"culture": i, "taxon": "t1", "metabolite": "x", "phase": "exponential", "change": v, "start": 0,
             "end": 12, "initial": 0, "cautions": []} for i, v in changes]


# ---- the spread, and the experiment as the unit -----------------------------------------------------

def test_the_mean_alone_and_one_value_decide_as_before():
    assert derive.classify([5.0, 5.2], 0.2) == 1
    assert derive.classify([-0.2, 1.9], 0.2) is None           # +0.84 +/- 1.04: inconclusive
    assert derive.classify([-0.2, 1.9], 0.2, agree=False) == 1
    assert derive.classify([0.9], 0.2) is None                 # one value: inconclusive
    assert derive.classify([0.9], 0.2, agree=False) == 1       # the mean alone, as 0.1.0


def test_each_experiment_counts_once_so_a_large_study_does_not_outvote_a_small_one():
    cultures = [_culture("E1", [0, 1, 2, 3]) for _ in range(10)] + [_culture("E2", [0, 1, 2, 3]) for _ in range(2)]
    cells = derive.pool(_rows([(i, 1.0) for i in range(10)] + [(10, -1.0), (11, -1.0)]), cultures, 0.2)
    cell = cells[("t1", "x", "exponential")]
    assert cell["mean"] == pytest.approx(0.0) and cell["n"] == 12 and cell["n_experiments"] == 2
    assert cell["state"] == "inconclusive" and "conflict" in cell["cautions"]


def test_identical_replicates_give_no_test():
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3])]
    cell = derive.pool(_rows([(0, 0.21), (1, 0.21)]), cultures, 0.2)[("t1", "x", "exponential")]
    assert cell["p_value"] is None and "no_variance" in cell["cautions"]


def test_a_compound_can_have_a_limit_of_its_own():
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3])]
    rows = _rows([(0, 0.05), (1, 0.06)])
    assert derive.pool(rows, cultures, 0.2)[("t1", "x", "exponential")]["state"] == "no_change"
    assert derive.pool(rows, cultures, 0.2, limits={"x": 0.01})[("t1", "x", "exponential")]["state"] == "produced"
    assert compound_limits_of({"compound_limits": "Thiamine=0.01; riboflavin: 0.005"}) == {
        "thiamine": 0.01, "riboflavin": 0.005}
    with pytest.raises(ValueError, match="name=mM"):
        compound_limits_of({"compound_limits": "thiamine"})


# ---- merging to genus --------------------------------------------------------------------------------

def _cell(mean, state, n=3):
    return {"mean": mean, "sd": None, "n": n, "n_experiments": 1, "state": state, "cautions": [], "notes": [],
            "experiments": ["E"], "studies": ["S"], "media": ["m"], "start": 0, "end": 1, "initial": 0,
            "exponential_h": 1, "direction": {"produced": "produced", "consumed": "consumed"}.get(state)}


def test_a_genus_whose_taxa_disagree_is_inconclusive_and_never_below_the_limit():
    taxa = {"a": {"name": "Bacteroides fragilis"}, "b": {"name": "Bacteroides thetaiotaomicron"}}
    cells = {("a", "x", "exponential"): _cell(0.30, "produced"), ("b", "x", "exponential"): _cell(-0.25, "consumed")}
    merged = derive.merge_genus_cells(cells, taxa, 0.2)[("genus:Bacteroides", "x", "exponential")]
    assert merged["direction"] is None and merged["state"] == "inconclusive" and "conflict" in merged["cautions"]
    cells = {("a", "x", "exponential"): _cell(0.30, "produced"), ("b", "x", "exponential"): _cell(0.1, "no_change")}
    merged = derive.merge_genus_cells(cells, taxa, 0.2)[("genus:Bacteroides", "x", "exponential")]
    assert merged["direction"] is None and "conflict" in merged["cautions"]


# ---- duplicate deposits ------------------------------------------------------------------------------

def test_two_experiments_that_only_sit_at_the_medium_level_are_not_one_deposit():
    flat = [_culture("E1", [30.0, 30.1, 29.9, 30.0]), _culture("E2", [30.05, 30.0, 30.0, 30.1], study="S2")]
    assert derive.duplicates(flat) == (set(), [])


def test_a_deposit_rounded_in_time_is_still_found():
    a = _culture("E1", [0.0, 2.0, 5.0, 9.0], times=(0, 8.333, 16.667, 25.0))
    b = _culture("E2", [0.0, 2.0, 5.0, 9.0], study="S2", times=(0, 8.33, 16.67, 25.0))
    exhausted = [_culture("E3", [10.0, 0.0, 0.0, 0.0]), _culture("E4", [10.0, 0.0, 0.0, 0.0], study="S2")]
    assert derive.duplicates(exhausted) == (set(), [])         # exhaustion alone identifies nothing
    dropped, found = derive.duplicates([a, b])
    assert dropped and found


# ---- growth phases -----------------------------------------------------------------------------------

def test_a_lone_late_jump_does_not_end_exponential_growth_at_the_last_point():
    times = [0, 4, 8, 12, 16, 20, 24]
    values = [1, 10, 100, 1000, 1000, 1000, 3000]
    assert phase.exponential_end(times, values)["end"] == 12


def test_a_boundary_on_too_few_samples_is_coarse():
    assert phase.exponential_end([0, 24, 48, 72], [0.7, 1.4, 1.4, 1.4])["coarse"] is True
    assert phase.exponential_end([0, 4, 8, 12, 24], [1, 10, 100, 1000, 1000])["coarse"] is False


def test_a_curve_that_gives_no_boundary_gives_way_to_the_next():
    short = {"times": [0, 8], "values": [1, 100], "technique": "fc", "level": "strain", "unit": ""}
    od = {"times": [0, 4, 8, 12, 24], "values": [0.01, 0.1, 0.4, 0.5, 0.5], "technique": "od",
          "level": "bioreplicate", "unit": ""}
    c = _culture("E1", [0, 1, 2, 3], growth=[short, od])
    found, _ = derive.boundaries([c])
    assert found[0]["end"] == 12 and c.growth is od


def test_relative_16s_is_no_growth_curve():
    context = {"subject": {"type": "strain"}, "techniqueType": "16S", "techniqueUnits": "% of reads"}
    assert _growth_rank(context) is None
    assert _growth_rank({**context, "techniqueUnits": "copies/mL"}) is not None


def test_replicates_that_end_growth_far_apart_are_marked():
    curve = {"times": [0, 4, 8, 12, 16], "technique": "fc", "level": "", "unit": ""}
    early = {**curve, "values": [1, 100, 1000, 1000, 1000]}
    late = {**curve, "values": [1, 10, 100, 300, 1000]}
    cultures = [_culture("E1", [0, 1, 2, 3], growth=[early]), _culture("E1", [0, 1, 2, 3], growth=[late])]
    found, _ = derive.boundaries(cultures)
    assert found[0]["differ"] and found[1]["differ"]


# ---- the value medium --------------------------------------------------------------------------------

def test_a_taxon_without_a_phase_gives_its_whole_run_change_and_does_not_vote(client):
    # C's cultures (mMCB) do not grow: its value is the change over the whole run (Karoline, 2026-10-06), marked
    # an OD read without its blank: from 0.5 to 0.65, 1.3-fold but a rise of 0.15, so it grew
    for rep in ("c1", "c1b"):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 0.5), (8, 0.55), (16, 0.65), (24, 0.6)]
    r = run(client, ignore_media=True)
    from foodnet import matrix
    assert matrix.entry(r, "ncbi:3", GLC, "exponential", "consumed") == (5.0, "whole_run")
    assert any("no end of exponential growth" in w for w in r["warnings"])
    arc = next(e for e in r["network"].edges if e.taxon == "ncbi:3" and e.metabolite == GLC)
    assert arc.phase == "whole_run" and "whole_run" in arc.cautions
    # with Both, the stationary column says no_phase
    r = run(client, ignore_media=True, phase="both")
    assert matrix.entry(r, "ncbi:3", GLC, "stationary", "consumed") == (None, "no_phase")
    # and it does not vote for the value medium: mMCB keeps two taxa with phase values only for B
    r = run(client)
    assert r["value_rule"]["media"] == ["Wilkins-Chalgren Anaerobe Broth (WC)"]


def test_a_culture_that_did_not_grow_gives_no_value(client):
    # flat OD (+0.05): grew by neither rule, so its glucose drift is no consumption
    for rep in ("c1", "c1b"):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 0.5), (8, 0.5), (16, 0.55), (24, 0.5)]
    r = run(client, ignore_media=True)
    from foodnet import matrix
    assert matrix.entry(r, "ncbi:3", GLC, "exponential", "consumed") == (None, "not_grown")
    assert not [e for e in r["network"].edges if e.taxon == "ncbi:3"]
    assert any("show no growth" in w and "Gamma" in w for w in r["warnings"])


def test_a_taxon_with_more_data_in_another_medium_is_named():
    import conftest
    from conftest import FakeClient
    series = dict(conftest.SERIES)
    series["b3"] = series["b2"]
    exp = dict(conftest.EXPERIMENTS["EMGDB000000003"], id="EMGDB000000007", bioreplicates=[{"id": "b3", "name": "b3"}])
    conftest.EXPERIMENTS["EMGDB000000007"] = exp
    conftest.STUDIES["SMGDB00000002"]["experiments"].append({"id": "EMGDB000000007", "name": "B_MCB2"})
    try:
        r = run(FakeClient(series))
        # mMCB now ties with WC on taxa and replicates and wins on its name; A has all its data in WC
        assert r["value_rule"]["media"] == ["mMCB"]
        assert any("Alpha alpha A1 has more data in Wilkins" in w for w in r["warnings"])
    finally:
        conftest.STUDIES["SMGDB00000002"]["experiments"].pop()
        del conftest.EXPERIMENTS["EMGDB000000007"]
    # in the plain world C has data only in mMCB, so it is named, and A and B (with data in WC) are not
    told = [w for w in run(FakeClient())["warnings"] if "has more data" in w]
    assert len(told) == 1 and "Gamma gamma C1 has more data in mMCB" in told[0] and "Alpha" not in told[0]


# ---- small samples and the interval (rounds 2 to 4) ---------------------------------------------------

def test_a_change_is_decided_on_a_confidence_interval_that_more_replicates_narrow():
    assert derive.classify([0.25, 0.3], 0.2) == 1                     # a pair: both beyond the limit
    assert derive.classify([-3.2, -1.52], 0.2) == -1                  # a clear pair is not inconclusive
    assert derive.classify([0.25, 0.3, 0.1], 0.2) is None              # three: the interval reaches inside
    assert derive.classify([-3.2, -1.52, -0.3], 0.2) == -1            # a confirming third never undoes a pair
    assert derive.classify([0.25, 0.3, 0.28, 0.27, 0.29], 0.2) == 1    # five pin it down beyond it
    assert derive.classify([0.05, -0.1, 0.12, 0.0], 0.2) == 0
    assert derive.classify([-0.4, 0.7], 0.2) is None
    # one failed sample does not erase what four replicates agree on (E. coli LF82 pyruvate)
    assert derive.classify([0.37, -8.56, -8.17, -8.06, -8.11], 0.2) == -1


def test_an_experiment_that_contradicts_the_others_is_a_conflict_even_when_noisy():
    cultures = [_culture("E1", [0, 1, 2, 3]) for _ in range(3)] + \
        [_culture("E2", [0, 1, 2, 3], study="S2") for _ in range(3)]
    rows = _rows([(0, 1.0), (1, 1.1), (2, 0.9), (3, -3.0), (4, -3.1), (5, 0.1)])
    cell = derive.pool(rows, cultures, 0.2)[("t1", "x", "exponential")]
    assert cell["state"] == "inconclusive" and "conflict" in cell["cautions"]


def test_experiments_that_agree_in_direction_but_not_in_size_stay_a_change():
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3]),
                _culture("E2", [0, 1, 2, 3], study="S2"), _culture("E2", [0, 1, 2, 3], study="S2")]
    cell = derive.pool(_rows([(0, 0.5), (1, 0.6), (2, 3.0), (3, 3.2)]), cultures, 0.2)[("t1", "x", "exponential")]
    assert cell["state"] == "produced" and "amounts_differ" in cell["cautions"]
    assert cell["mean"] == pytest.approx((0.55 + 3.1) / 2)


def test_an_inconclusive_experiment_does_not_veto_the_others_but_is_named():
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3]),
                _culture("E2", [0, 1, 2, 3], study="S2"), _culture("E2", [0, 1, 2, 3], study="S2")]
    cell = derive.pool(_rows([(0, 1.0), (1, 1.2), (2, -0.1), (3, 1.5)]), cultures, 0.2)[("t1", "x", "exponential")]
    # E1 (+1.0, +1.2) decides; E2 (-0.1, +1.5) is inconclusive and does not contradict. The amount is over both
    # experiments (+1.1 and +0.7), since leaving out the smaller effect would inflate it
    assert cell["state"] == "produced" and cell["mean"] == pytest.approx(0.9)
    assert any("left out as inconclusive" in n for n in cell["notes"]) and "experiment_left_out" in cell["cautions"]
    assert cell["n"] == 4 and cell["n_experiments"] == 2


def test_identical_replicates_count_as_one_measurement():
    assert derive.classify([0.21, 0.21], 0.2) is None        # one series twice: one value, inconclusive
    assert derive.classify([0.0, 0.0], 0.2) == 0              # but an exhausted compound is exactly 0 twice
    assert derive.classify([-8.23, -8.23, -8.23], 0.2) == -1   # and a substrate exhausted from one start
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3])]
    cell = derive.pool(_rows([(0, 0.21), (1, 0.21)]), cultures, 0.2)[("t1", "x", "exponential")]
    assert "no_variance" in cell["cautions"]


def test_repeats_of_one_protocol_that_differ_by_more_than_rounding_are_not_one_deposit():
    a = _culture("E1", [27.5, 14.0, 5.0, 0.0], times=(0, 8, 24, 48))
    b = _culture("E2", [27.5, 14.1, 5.02, 0.0], study="S2", times=(0, 8, 24, 48))
    assert derive.duplicates([a, b]) == (set(), [])


def test_a_zero_decided_on_a_pair_is_flagged_too():
    cultures = [_culture("E1", [0, 1, 2, 3]), _culture("E1", [0, 1, 2, 3])]
    cell = derive.pool(_rows([(0, 0.15), (1, 0.18)]), cultures, 0.2)[("t1", "x", "exponential")]
    assert cell["state"] == "no_change" and "pair_decided" in cell["cautions"]


def test_replicates_that_start_far_apart_are_flagged():
    cultures = [_culture("E1", [0, 1, 2, 3]) for _ in range(3)]
    rows = _rows([(0, 1.0), (1, 1.1), (2, 0.9)])
    rows[2]["start"] = 4                       # one replicate's first sample is at 4 h of a 0 to 12 h phase
    cell = derive.pool(rows, cultures, 0.2)[("t1", "x", "exponential")]
    assert "start_differs" in cell["cautions"]


def test_a_taxon_with_phase_cultures_takes_its_biomass_from_them_only():
    from foodnet import crm
    curve = {"times": [0, 4, 8, 12], "technique": "od", "level": "", "unit": "od"}
    grew = _culture("E1", [0, 1, 2, 3], growth=[{**curve, "values": [0.1, 0.4, 1.0, 1.0]}])
    flat = _culture("E2", [0, 1, 2, 3], growth=[{**curve, "values": [0.7, 0.7, 0.71, 0.7]}])
    rows = [{"culture": 0, "phase": "exponential", "change": 2.0, "start": 0, "end": 8, "cautions": []},
            {"culture": 1, "phase": "exponential", "change": 1.0, "start": 0, "end": 12, "cautions": ["whole_run"]}]
    found = crm.biomass_changes([(0, grew), (1, flat)], rows, "exponential")["t1"]
    assert found["n"] == 1 and found["change"] == pytest.approx(0.9) and found["unit"] == "od"


def test_a_culture_that_shows_no_growth_but_metabolizes_gives_its_whole_run_change():
    curve = {"times": [0, 4, 8, 12], "technique": "od", "level": "", "unit": "od", "values": [0.63, 0.62, 0.63, 0.62]}
    c = _culture("E1", [0, 1, 2, 3], growth=[curve])
    def met(name, values):
        return {"name": name, "chebi_id": "", "series": list(zip((0, 4, 8, 12), values, strict=True))}
    c.metabolites = {"lactate": met("lactate", (1.0, 0.6, 0.1, 0.0)), "butyrate": met("butyrate", (1.2, 2.0, 3.1, 3.8))}
    rows, _ = derive.changes([c], "exponential")
    assert all(r["cautions"][:2] == ["whole_run", "growth_unclear"] and r["change"] is not None for r in rows)
    assert c.not_grown


def test_a_fold_on_near_blank_od_is_no_growth_but_low_density_growth_is():
    with pytest.raises(phase.NoBoundary, match="did not grow"):
        phase.exponential_end([0, 4, 8, 12, 16, 24], [0.020, 0.024, 0.019, 0.031, 0.027, 0.022], smoothed=True)
    # Variovorax-like: from 0.002, rising over several samples
    assert phase.exponential_end([0, 4, 8, 12, 16], [0.002, 0.006, 0.02, 0.04, 0.045], smoothed=True)["end"] == 16


def test_evaporation_is_not_metabolism():
    curve = {"times": [0, 4, 8, 12], "technique": "od", "level": "", "unit": "od", "values": [0.5, 0.5, 0.5, 0.5]}
    c = _culture("E1", [0, 1, 2, 3], growth=[curve])
    # 10% volume loss: lactate concentrates steadily, volatile ethanol falls steadily
    def met(name, values):
        return {"name": name, "chebi_id": "", "series": list(zip((0, 4, 8, 12), values, strict=True))}
    c.metabolites = {"lactate": met("lactate", (20, 20.7, 21.5, 22.2)), "ethanol": met("ethanol", (5, 4.7, 4.3, 4.0))}
    assert not derive.active(c)


def test_a_volatile_falling_is_no_uptake_even_past_the_evaporation_share():
    curve = {"times": [0, 4, 8, 12], "technique": "od", "level": "", "unit": "od", "values": [0.5, 0.5, 0.5, 0.5]}
    c = _culture("E1", [0, 1, 2, 3], growth=[curve])

    def met(name, values):
        return {"name": name, "chebi_id": "", "series": list(zip((0, 4, 8, 12), values, strict=True))}
    # 12% volume loss: lactate concentrates past 10%, ethanol evaporates
    c.metabolites = {"lactate": met("lactate", (20, 20.8, 21.6, 22.4)), "ethanol": met("ethanol", (5, 4.4, 3.8, 3.0))}
    assert not derive.active(c)
