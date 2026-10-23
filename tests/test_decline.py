"""The decline after the growth peak, for a consumer-resource fit's biomass loss rate (Karoline, 2026-10-09)."""
import io
import math
import zipfile

import pytest
from conftest import run

from foodnet import crm, matrix, rates


def test_a_death_phase_gives_its_log_slope_from_the_peak_to_the_last_sample():
    times = [0, 4, 8, 12, 24, 36, 48]
    values = [1, 10, 100, 1000, 500, 250, 125]               # halves every 12 h after the peak
    d = rates.decline(times, values)
    assert d["rate"] == pytest.approx(math.log(2) / 12)
    assert (d["start"], d["end"], d["points"], d["death"]) == (12, 48, 4, True)


def test_a_plateau_gives_a_small_rate_that_only_bounds_death():
    times = [0, 4, 8, 12, 24, 36, 48]
    values = [1, 10, 100, 1000, 990, 980, 970]               # stays within the plateau band
    d = rates.decline(times, values)
    assert 0 < d["rate"] < 0.001 and d["death"] is False


def test_no_usable_decline_is_absent_not_zero():
    with pytest.raises(rates.RateUnavailable):
        rates.decline([0, 4, 8, 12, 24], [1, 10, 100, 1000, 1000])   # the peak one sample before the end
    with pytest.raises(rates.RateUnavailable):
        rates.decline([0, 4, 8], [1, 10, 100])                       # still rising at the end


def test_a_taxon_summary_carries_scatter_n_phase_and_window():
    found = [{"rate": r, "start": 10.0, "end": 48.0, "points": 4, "death": death, "study": "S"}
             for r, death in ((0.006, True), (0.007, True), (0.008, False))]
    s = crm.summarize_declines(found, crm.FROM_METABOLITE_REPLICATES)
    assert s["rate"] == pytest.approx(0.007) and s["n"] == 3
    assert s["sd"] == pytest.approx(0.001) and s["se"] == pytest.approx(0.001 / math.sqrt(3))
    assert (s["phase"], s["death_replicates"], s["start_h"], s["end_h"]) == ("death", 2, 10.0, 48.0)
    one = crm.summarize_declines(found[:1], crm.FROM_SAME_MEDIUM)
    assert one["sd"] is None and one["se"] is None


def test_the_crm_parameters_carry_decline_rates_or_say_why_not(client):
    # the test curves end on their peak's plateau one sample after it, so no taxon has a usable decline
    r = run(client, report_rates=True)
    assert r["declines"] == {} and r["without_a_decline"]
    p = matrix.crm_payload(r)
    assert all(v is None for v in p["decline_rates"]) and p["decline_rate_detail"] == {}
    table = zipfile.ZipFile(io.BytesIO(matrix.crm_package(r))).read("decline_rates.csv").decode()
    assert table.splitlines()[0].startswith("taxon,decline_rate,unit,se,sd,replicates,experiments,phase")
    assert all(row.split(",")[1] == "NA" for row in table.splitlines()[1:])


def test_a_declining_taxon_gets_its_rate_from_its_metabolite_replicates(client):
    A = "ncbi:1"
    for i, rep in enumerate(("a1", "a2", "a3")):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 1), (4, 10), (8, 100), (12, 1000), (24, 500 + 10 * i), (36, 250), (48, 125)]
    r = run(client, report_rates=True)
    d = r["declines"][A]
    assert d["n"] == 3 and d["phase"] == "death" and d["source"] == crm.FROM_METABOLITE_REPLICATES
    assert d["se"] is not None and (d["start_h"], d["end_h"]) == (12, 48)
    p = matrix.crm_payload(r)
    a = p["taxa"].index("Alpha alpha A1")
    assert p["decline_rates"][a] == pytest.approx(d["rate"]) and "replicates" not in p["decline_rate_detail"][
        "Alpha alpha A1"]


def test_a_plateau_runs_from_the_peak_to_its_last_sample_within_the_band():
    assert rates.plateau([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 990, 980, 100]) == (12, 36)
    with pytest.raises(rates.RateUnavailable):
        rates.plateau([0, 4, 8, 12, 24], [1, 10, 100, 1000, 100])    # falls right after the peak


def test_the_maintenance_coefficient_is_uptake_per_biomass_per_hour_on_the_plateau(client):
    # A's curves sit at 1000 from 12 h to 24 h; glucose falls by 1.9-1.9=0, 2.1-1.1=1 and 2-1.5=0.5 mM there
    r = run(client, report_rates=True)
    A, GLC = "ncbi:1", next(m for m, n in r["names"].items() if n == "glucose")
    m = r["maintenance"][A][GLC]
    assert m["coefficient"] == pytest.approx(0.5 / 1000 / 12) and m["n"] == 3
    assert (m["start_h"], m["end_h"]) == (12, 24) and m["unit"] == "mM h-1 per (fc)" and m["se"] is not None
    # glucose at the plateau's start: 1.9, 2.1 and 2 mM, so it was there to take up; two growth samples on it
    assert m["start_mM"] == pytest.approx(2.0) and m["below_limit_replicates"] == 0
    assert (m["plateau_samples"], m["n_experiments"]) == (2, 1)
    p = matrix.crm_payload(r)
    row = next(x for x in p["maintenance"] if x["taxon"] == "Alpha alpha A1" and x["resource"] == "glucose")
    assert row["coefficient"] == pytest.approx(m["coefficient"]) and "cultures" not in row
    table = zipfile.ZipFile(io.BytesIO(matrix.crm_package(r))).read("maintenance.csv").decode()
    assert table.splitlines()[0].startswith("taxon,resource,maintenance,unit,se,sd,replicates")
    # B's glucose does not change on its plateau (16 h to 24 h): a measured 0, not absent
    B = "ncbi:2"
    assert r["maintenance"][B][GLC]["coefficient"] == 0


def test_rates_name_the_growth_curve_technique_and_warn_when_taxa_differ(client):
    # Karoline, 2026-10-22: "the measurement technique, to prevent mismatches across techniques"; the fixture's
    # A is counted by flow cytometry, B and C read by OD (a window of 3 gives every taxon a rate)
    r = run(client, report_rates=True, rate_window=3)
    A, B = "ncbi:1", "ncbi:2"
    assert r["rates"][A]["technique"].startswith("fc") and r["rates"][B]["technique"].startswith("od")
    assert matrix.rates_csv(r).splitlines()[0].endswith(",technique,other_techniques_left_out")
    assert ",technique," in matrix.declines_csv(r).splitlines()[0]
    warning = matrix.technique_mismatch(r)
    assert warning.startswith("Growth curves of different techniques: fc") and "od for" in warning
    assert warning in matrix.readme(r, "crm")


def test_the_decline_starts_at_the_end_of_growth_not_at_the_highest_sample():
    # Karoline, 2026-10-22: "First sample >= 90% of max": 1 + 0.9 * 1299 = 1170 is first reached at 24 h, not at
    # the band's lower edge (12 h, still growing) nor on the noise that made 24 h the highest
    d = rates.decline([0, 4, 8, 12, 24, 36, 48, 60], [1, 10, 100, 800, 1300, 1000, 1000, 1000])
    assert d["start"] == 24 and d["points"] == 4
    # a curve that falls to zero is fitted up to the zero and counts as death
    z = rates.decline([0, 4, 8, 12, 24, 36, 48, 60], [1, 10, 100, 1000, 500, 250, 125, 0])
    assert z["to_zero"] and z["death"] and z["end"] == 48 and z["rate"] == pytest.approx(math.log(2) / 12)
    with pytest.raises(rates.RateUnavailable, match="to the zero"):
        rates.decline([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 100, 0, 0])


def test_of_a_duplicate_pair_the_richer_growth_curve_gives_the_rates():
    # Karoline, 2026-10-22: "Richer growth curve": the metabolite rule kept E1, whose curve lacks the 8 h sample
    import types

    from foodnet import search

    def culture(experiment, times):
        return types.SimpleNamespace(experiment=experiment, growth={"times": times,
                                                                    "values": [10 ** i for i in range(len(times))]})
    sparse = [culture("E1", [0, 4, 12, 16]) for _ in range(3)]
    rich = [culture("E2", [0, 4, 8, 12, 16]) for _ in range(3)]
    out, maint_out = search._growth_copies([(sparse, rich, ["acetate", "glucose"])])
    assert out == {"E1"} and maint_out == {("E1", "acetate"), ("E1", "glucose")}
    # a tie keeps the metabolite rule's copy
    assert search._growth_copies([(rich, [culture("E3", [0, 4, 8, 12, 16])], ["acetate"])])[0] == {"E3"}


def test_a_left_out_copy_gives_no_growth_rate_or_decline(client, monkeypatch):
    # pins the search's use of the copies: leaving A's only experiment out leaves A without either
    from foodnet import search
    monkeypatch.setattr(search, "_growth_copies", lambda pairs, fraction=0.9, factor=1.5: ({"EMGDB000000001"}, set()))
    r = run(client, report_rates=True, rate_window=3)
    assert "ncbi:1" not in r["rates"] and "ncbi:1" not in r["declines"]


def test_a_tie_between_units_always_goes_to_the_preferred_technique():
    # a review: a tie went by string hashing, which changes from run to run
    assert crm.majority(["od", "fc Cells/mL", "od", "fc Cells/mL"]) == "fc Cells/mL"
    assert crm.majority(["od", "od", "fc Cells/mL"]) == "od"


def test_declines_keep_the_technique_most_replicates_use():
    found = [{"rate": r, "start": 12.0, "end": 48.0, "points": 4, "death": True, "study": "S", "experiment": e,
              "technique": t} for r, e, t in ((0.05, "E1", "fc Cells/mL (strain)"),
                                              (0.06, "E1", "fc Cells/mL (strain)"), (0.0, "E2", "od (culture)"))]
    s = crm.summarize_declines(found, crm.FROM_METABOLITE_REPLICATES)
    assert s["rate"] == pytest.approx(0.055) and s["technique"] == "fc Cells/mL (strain)"
    assert (s["n"], s["n_experiments"], s["other_techniques_left_out"]) == (2, 1, 1)


def _plateau_culture(experiment, glucose, not_grown=False, technique="fc"):
    import types
    return types.SimpleNamespace(
        growth={"times": [0, 4, 8, 12, 24], "values": [1, 10, 100, 1000, 1000], "technique": technique,
                "unit": "", "level": "culture"},
        not_grown=not_grown, taxon={"id": "ncbi:9"}, experiment=experiment, study="S", label=experiment,
        metabolites={"glc": {"series": [(0, 10), (12, glucose[0]), (24, glucose[1])]}})


def test_maintenance_leaves_out_repeats_and_cultures_that_did_not_grow_and_flags_a_used_up_compound():
    cultures = [_plateau_culture("E1", (2.0, 1.0)), _plateau_culture("E1", (2.2, 1.0)),
                _plateau_culture("DUP", (2.0, 1.0)), _plateau_culture("E3", (2.0, 1.0), not_grown=True),
                _plateau_culture("E4", (2.0, 1.0), technique="od")]
    found, skipped = crm.maintenance(cultures, 0.2, dropped={("DUP", "glc")})
    m = found["ncbi:9"]["glc"]
    assert (m["n"], m["n_experiments"], m["biomass_unit"]) == (2, 1, "fc")      # not DUP, E3, nor E4's OD
    assert m["other_units_left_out"] == 1                                   # E4's OD culture, counted
    used_up, _ = crm.maintenance([_plateau_culture("E1", (0.05, 0.0)), _plateau_culture("E1", (0.1, 0.0))], 0.2)
    gone = used_up["ncbi:9"]["glc"]
    assert gone["coefficient"] == 0 and gone["below_limit_replicates"] == 2 and gone["start_mM"] == pytest.approx(0.075)


def test_fit_hours_span_every_replicates_fit():
    # a review: means of 2-6 h and 20-30 h named 11-18 h, which no fit used
    import types
    curves = [([0, 2, 4, 6, 8, 10], [1, 10, 100, 1000, 1000, 1000]),
              ([0, 10, 20, 25, 30, 40], [1, 1, 10, 100, 1000, 1000])]
    cultures = [types.SimpleNamespace(taxon={"id": "t"}, study="S", label=str(i), growth={
        "times": t, "values": v, "technique": "od", "unit": "", "level": "culture"}) for i, (t, v) in enumerate(curves)]
    found, _, _ = crm.growth_rates(["t"], cultures, [], set(), window=3)
    assert (found["t"]["fit_start_h"], found["t"]["fit_end_h"]) == (0, 30)


def test_merged_genera_keep_every_strains_studies(client):
    r = run(client, report_rates=True, rate_window=3, merge_genera=True)
    assert all(set(v["studies"]) for v in r["rates"].values())
    from foodnet import search
    merged, _ = search._genus_rates({"a": {"rate": 0.1, "n": 1, "studies": ["S1"], "source": "x", "fit_start_h": 0,
                                           "fit_end_h": 8, "technique": "fc (culture)"},
                                     "b": {"rate": 0.3, "n": 2, "studies": ["S2"], "source": "y", "fit_start_h": 2,
                                           "fit_end_h": 12, "technique": "fc (strain)"},
                                     "c": {"rate": 0.01, "n": 1, "studies": ["S3"], "source": "z", "fit_start_h": 0,
                                           "fit_end_h": 40, "technique": "od (culture)"}},
                                    {}, {t: {"genus": "G", "name": f"G {t}"} for t in "abc"})
    g = next(iter(merged.values()))
    assert g["studies"] == ["S1", "S2"] and g["source"] == "x; y" and (g["fit_start_h"], g["fit_end_h"]) == (0, 12)
    # the OD strain is of another technique than most of the genus's replicates (a review: a Bacteroides median
    # of flow cytometry and OD strains)
    assert g["rate"] == pytest.approx(0.2) and g["other_techniques_left_out"] == 1


def test_growth_curves_judged_unusable_give_no_maintenance(client, monkeypatch):
    # an OD read without its blank: C's curves give no end of exponential growth, so its values span the whole
    # run; as for its growth rate and decline, it gets no maintenance (a review)
    for rep in ("c1", "c1b"):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 0.5), (8, 0.55), (16, 0.65), (24, 0.6)]
    r = run(client, ignore_media=True, report_rates=True)
    assert "ncbi:3" not in r["maintenance"] and "ncbi:3" not in r["declines"]


def _curve_culture(experiment, values, technique="fc", level="strain", not_grown=False, times=None):
    import types
    return types.SimpleNamespace(
        growth={"times": times or [0, 4, 8, 12, 24, 36, 48][:len(values)], "values": values, "technique": technique,
                "unit": "", "level": level},
        not_grown=not_grown, taxon={"id": "t"}, experiment=experiment, study="S", label=experiment,
        metabolites={}, medium_key="m")


DYING = [1, 10, 100, 1000, 500, 250, 125]


def test_a_culture_that_did_not_grow_gives_no_decline():
    found, _, skipped = crm.decline_rates(["t"], [_curve_culture("E1", DYING, not_grown=True)], [], {"m"})
    assert "t" not in found and any("did not grow" in why for _, why in skipped)


def test_another_monoculture_only_where_the_value_replicates_have_no_growth_curve():
    flat = [1, 10, 100, 1000, 1000, 1000, 1000]
    short = _curve_culture("E1", [1, 10, 100, 1000, 500])         # a curve, but too short for a decline
    found, missing, _ = crm.decline_rates(["t"], [short], [_curve_culture("E9", flat)], {"m"})
    assert "t" not in found and "t" in missing
    import types
    curveless = types.SimpleNamespace(growth=None, taxon={"id": "t"}, experiment="E1", medium_key="m")
    found, _, _ = crm.decline_rates(["t"], [curveless], [_curve_culture("E9", flat)], {"m"})
    assert found["t"]["source"] == crm.FROM_SAME_MEDIUM


def test_the_fraction_setting_moves_the_start_of_the_decline():
    curve = [1, 10, 100, 800, 1000, 900, 800]
    early = crm.decline_rates(["t"], [_curve_culture("E1", curve)], [], {"m"}, fraction=0.7)
    late = crm.decline_rates(["t"], [_curve_culture("E1", curve)], [], {"m"}, fraction=0.9)
    assert early[0]["t"]["start_h"] == 12 and late[0]["t"]["start_h"] == 24


def test_a_tie_goes_to_the_preferred_technique_whatever_the_name():
    assert crm.majority(["od", "fc"]) == "fc" and crm.majority(["cfu", "plates"]) == "plates"


def test_a_genus_decline_keeps_strain_and_culture_level_counts_of_one_technique():
    # a review: B. fragilis' culture-level counts were left out of a Bacteroides decline as another technique
    found = [{"rate": r, "start": 12.0, "end": 48.0, "points": 4, "death": True, "study": "S", "experiment": e,
              "technique": t} for r, e, t in ((0.05, "E1", "fc Cells/mL (strain)"),
                                              (0.01, "E2", "fc Cells/mL (culture)"), (0.0, "E3", "od (culture)"))]
    s = crm.summarize_declines(found, crm.FROM_METABOLITE_REPLICATES)
    assert s["n"] == 2 and s["other_techniques_left_out"] == 1
    assert s["technique"] == "fc Cells/mL (culture); fc Cells/mL (strain)"


def test_initial_concentrations_count_a_duplicate_once():
    import types

    def c(experiment, start):
        return types.SimpleNamespace(experiment=experiment, metabolites={"glc": {"name": "glucose",
                                                                                 "series": [(0, start), (8, 1)]}})
    found = crm.initial_concentrations([c("E1", 10), c("DUP", 10), c("E2", 8)], dropped={("DUP", "glc")})
    assert found["glc"]["n"] == 2 and found["glc"]["mean"] == pytest.approx(9)


def test_biomass_takes_the_preferred_technique_on_a_tie_too():
    import types
    cultures = []
    for i, tech in enumerate(("od", "fc")):
        cultures.append((i, types.SimpleNamespace(
            growth={"times": [0, 12], "values": [1, 100], "technique": tech, "unit": ""}, taxon={"id": "t"},
            experiment=f"E{i}")))
    rows = [{"culture": i, "phase": "exponential", "change": -1.0, "start": 0, "end": 12, "cautions": []}
            for i in range(2)]
    assert crm.biomass_changes(cultures, rows, "exponential")["t"]["unit"] == "fc"


def test_a_single_zero_followed_by_counts_is_a_dropout_not_a_fall_to_zero():
    d = rates.decline([0, 4, 8, 12, 24, 36, 48, 60], [1, 10, 100, 1000, 1000, 0, 1000, 1000])
    assert not d["to_zero"] and not d["death"] and d["end"] == 60 and d["rate"] == pytest.approx(0, abs=1e-12)


def test_growth_rates_take_the_technique_most_replicates_use():
    # Karoline, 2026-10-22: "growth rates should follow the same rules"
    import types

    def c(i, technique, values):
        return types.SimpleNamespace(taxon={"id": "t"}, study=f"S{i}", label=str(i), growth={
            "times": [0, 2, 4, 6, 8, 10], "values": values, "technique": technique, "unit": "", "level": "culture"})
    fast, slow = [1, 10, 100, 1000, 1000, 1000], [1, 2, 4, 8, 16, 16]
    found, _, _ = crm.growth_rates(["t"], [c(0, "fc", fast), c(1, "fc", fast), c(2, "od", slow)], [], set(), window=3)
    assert found["t"]["n"] == 2 and found["t"]["other_techniques_left_out"] == 1
    assert found["t"]["technique"] == "fc (culture)" and found["t"]["studies"] == ["S0", "S1"]


def test_the_richer_copy_is_counted_to_the_end_of_growth_not_to_a_late_maximum():
    # a review: a copy sampled once during growth, whose plateau crept up, outranked one sampled every 2 h
    import types

    from foodnet import search

    def culture(experiment, times, values):
        return types.SimpleNamespace(experiment=experiment, growth={"times": times, "values": values})
    dense = [culture("A", [0, 2, 4, 6, 8, 10, 12, 24, 48], [1, 3, 9, 27, 81, 243, 700, 720, 730])]
    sparse = [culture("B", [0, 8, 14, 20, 26, 32, 38, 44, 50], [1, 300, 700, 710, 720, 730, 740, 750, 760])]
    assert search._growth_copies([(sparse, dense, ["acetate"])])[0] == {"B"}
    assert search._growth_copies([(dense, sparse, ["acetate"])])[0] == {"B"}


def test_a_copy_outside_the_search_never_takes_the_rates_of_one_inside_it(client, monkeypatch):
    # a review: with --outside-evidence, study 7's copy took B. thetaiotaomicron's rates from study 2's
    from foodnet import derive, search
    real = derive.duplicate_pairs

    def pairs(cultures):
        a = [c for c in cultures if c.experiment == "EMGDB000000001"]
        c4 = [c for c in cultures if c.experiment == "EMGDB000000004"]
        return real(cultures) + ([(a, c4, [])] if a and c4 else [])
    monkeypatch.setattr(search.d, "duplicate_pairs", pairs)
    monkeypatch.setattr(search, "_richness", lambda cultures, *a: (9.0, 9.0) if cultures[0].experiment ==
                        "EMGDB000000004" else (1.0, 1.0))
    r = run(client, report_rates=True, rate_window=3, conditions="SMGDB00000001", outside_evidence=True)
    assert "ncbi:1" in r["rates"] and r["rates"]["ncbi:1"]["source"] == crm.FROM_METABOLITE_REPLICATES


def test_the_readme_and_report_state_the_searchs_fraction(client):
    from foodnet.report import report_text
    r = run(client, report_rates=True, fraction=0.8)
    assert "80% in this search" in matrix.readme(r, "crm") and "80% of the rise" in report_text(r)


def test_a_flat_fit_reads_zero_not_minus_zero():
    d = rates.decline([0, 4, 8, 12, 24, 36], [1, 10, 100, 1000, 1000, 1000])
    assert str(d["rate"]) == "0.0"


def test_biomass_takes_the_preferred_technique_where_the_alphabet_disagrees():
    # plates before cfu by preference, cfu before plates by name
    import types
    cultures = [(i, types.SimpleNamespace(growth={"times": [0, 12], "values": [1, 100], "technique": tech, "unit": ""},
                                          taxon={"id": "t"}, experiment=f"E{i}"))
                for i, tech in enumerate(("cfu", "plates"))]
    rows = [{"culture": i, "phase": "exponential", "change": -1.0, "start": 0, "end": 12, "cautions": []}
            for i in range(2)]
    assert crm.biomass_changes(cultures, rows, "exponential")["t"]["unit"] == "plates"


def test_r_writes_the_same_columns_as_foodnet(client):
    # the R test pins crm_write's headers; these are foodnet's own, so the two stay the same (a review)
    import pathlib
    import re
    r = run(client, report_rates=True)
    rtest = (pathlib.Path(__file__).parents[1] / "r" / "tests" / "testthat" / "test-payload.R").read_text()
    written = ",".join(re.findall(r'"([A-Za-z0-9_]+)"', rtest[rtest.index("writes the download's columns"):]))
    for table in (matrix.rates_csv, matrix.declines_csv, matrix.maintenance_csv, matrix.biomass_csv):
        names = table(r).splitlines()[0].split(",")
        assert ",".join(names) in written, names


def test_a_curve_whose_highest_sample_is_its_last_still_has_a_plateau():
    # Karoline, 2026-10-22: the plateau starts where the decline starts, so a creeping plateau is not left out
    assert rates.plateau([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 990, 980, 1010]) == (12, 48)


def test_a_taxons_maintenance_unit_is_the_one_most_of_its_cultures_use_whatever_they_measured():
    # a review: two OD cultures measuring many compounds outvoted three qPCR cultures measuring one
    many = _plateau_culture("E1", (2.0, 1.0), technique="od")
    many.metabolites.update({f"m{i}": {"series": [(0, 10), (12, 2.0), (24, 1.0)]} for i in range(10)})
    few = [_plateau_culture(f"E{i}", (2.0, 1.0), technique="qpcr") for i in (2, 3)]
    for c in few:
        c.label = c.experiment
    found, skipped = crm.maintenance([many] + few, 0.2)
    assert found["ncbi:9"]["glc"]["biomass_unit"] == "qpcr" and "m0" not in found["ncbi:9"]
    assert any("in od" in why for _, why in skipped)


def test_the_fraction_setting_moves_the_maintenance_plateau():
    curve = _plateau_culture("E1", (2.0, 1.0))
    curve.growth = {"times": [0, 4, 8, 12, 24], "values": [1, 10, 100, 800, 1000], "technique": "fc", "unit": "",
                    "level": "culture"}
    early, _ = crm.maintenance([curve], 0.2, fraction=0.7)
    late, _ = crm.maintenance([curve], 0.2, fraction=0.9)
    assert early["ncbi:9"]["glc"]["start_h"] == 12 and "ncbi:9" not in late


def test_fit_r2_is_one_for_a_clean_exponential_decline():
    d = rates.decline([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 500, 250, 125])
    assert d["r2"] == pytest.approx(1.0)


def test_growth_rates_leave_out_cultures_that_did_not_grow(client):
    # Karoline, 2026-10-22: growth rates "should follow the same rules"; B's mMCB cultures read OD without
    # growth, beside its WC cultures that grew
    for rep in ("b2", "b2b"):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 0.5), (8, 0.52), (16, 0.55), (24, 0.56)]
    r = run(client, ignore_media=True, report_rates=True, rate_window=3)
    assert r["rates"]["ncbi:2"]["n"] == 2


def test_a_dead_cell_count_is_never_a_growth_curve_and_a_repeated_gate_is_named():
    from foodnet import reading
    dead = {"subject": {"type": "bioreplicate"}, "techniqueType": "fc", "techniqueCellType": "dead"}
    live = {"subject": {"type": "bioreplicate"}, "techniqueType": "fc", "techniqueCellType": "live"}
    assert reading._growth_rank(dead) is None and reading._growth_rank(live) is not None


def test_the_cell_gate_is_part_of_the_technique():
    g = {"technique": "fc", "unit": "Cells/mL", "level": "strain", "gate": "live"}
    assert crm.technique(g) == "fc Cells/mL live (strain)"
    found = [{"rate": r, "start": 12.0, "end": 48.0, "points": 4, "death": True, "study": "S", "experiment": e,
              "technique": t} for r, e, t in ((0.05, "E1", "fc Cells/mL live (strain)"),
                                              (0.05, "E2", "fc Cells/mL live (strain)"),
                                              (0.002, "E3", "fc Cells/mL (culture)"))]
    s = crm.summarize_declines(found, crm.FROM_METABOLITE_REPLICATES)
    assert s["n"] == 2 and s["other_techniques_left_out"] == 1


def test_maintenance_says_how_much_the_biomass_changed_and_when_a_compound_ran_out_in_the_window():
    c = _plateau_culture("E1", (2.0, 0.1))
    c.growth["values"] = [1, 10, 100, 1000, 600]
    c.growth["times"] = [0, 4, 8, 12, 24]
    found, _ = crm.maintenance([c], 0.2, fraction=0.9)
    m = found["ncbi:9"]["glc"]
    assert m["biomass_relative_change"] == pytest.approx(-0.4) and m["below_limit_by_end"] == 1
    assert m["replicate_values"] == [pytest.approx(m["coefficient"])]


def test_an_untagged_copy_of_a_live_count_is_named_live_and_one_of_a_dead_count_is_dropped():
    # a review: the naming compared only against the curves kept, so an untagged copy of a dead count stayed
    from foodnet import reading
    times = [0, 4, 8]
    live = {"times": times, "values": [1, 10, 100], "gate": "live"}
    copy_of_live = {"times": times, "values": [1, 10, 100], "gate": ""}
    copy_of_dead = {"times": times, "values": [0, 5, 50], "gate": ""}
    kept = reading.name_gates([copy_of_dead, copy_of_live, live], [list(zip(times, [0, 5, 50], strict=True))])
    assert kept == [copy_of_live, live] and copy_of_live["gate"] == "live"


def test_a_genus_whose_taxa_are_all_inconclusive_draws_and_downloads(client):
    # a review: a merged genus cell with no mean crashed the figure, the matrices zip and the page (as in 0.4.0)
    from foodnet.figure import matrices_svg
    r = run(client, merge_genera=True)
    key = next(k for k, c in r["cells"].items() if c.get("n"))
    r["cells"][key].update({"mean": None, "state": "inconclusive", "sd": None, "direction": None})
    assert "<svg" in matrices_svg(r)
    assert zipfile.ZipFile(io.BytesIO(matrix.pair_package(r))).namelist()


def test_technique_mismatch_keys_on_technique_unit_and_gate_not_level():
    r = {"names": {}, "rates": {"a": {"technique": "fc Cells/mL live (strain)"},
                                "b": {"technique": "fc Cells/mL live (culture)"}}, "declines": {}, "biomass": {}}
    assert matrix.technique_mismatch(r) == ""
    r["rates"]["b"]["technique"] = "fc Cells/mL (culture)"
    assert matrix.technique_mismatch(r).startswith("Growth curves of different techniques")


def test_maintenance_csv_writes_the_relative_change_and_the_compounds_gone_in_the_window(client):
    r = run(client, report_rates=True)
    header, *rows = matrix.maintenance_csv(r).splitlines()
    names = header.split(",")
    alpha = next(row.split(",") for row in rows if row.startswith("Alpha alpha A1,glucose"))
    assert alpha[names.index("biomass_relative_change")] == "0" and alpha[names.index("below_limit_by_end")] == "0"


def test_a_compound_already_below_the_limit_at_the_start_is_not_counted_as_gone_in_the_window():
    c = _plateau_culture("E1", (0.1, 0.0))
    found, _ = crm.maintenance([c], 0.2)
    assert found["ncbi:9"]["glc"]["below_limit_by_end"] == 0 and found["ncbi:9"]["glc"]["below_limit_replicates"] == 1


def test_crm_json_keeps_each_replicates_decline_rate(client):
    for i, rep in enumerate(("a1", "a2", "a3")):
        client.contexts[client.bioreplicates[rep]["measurementContexts"][0]["id"]] = \
            [(0, 1), (4, 10), (8, 100), (12, 1000), (24, 500 + 10 * i), (36, 250), (48, 125)]
    p = matrix.crm_payload(run(client, report_rates=True))
    rates_ = p["decline_rate_detail"]["Alpha alpha A1"]["replicate_rates"]
    assert len(rates_) == 3 and sum(rates_) / 3 == pytest.approx(p["decline_rates"][p["taxa"].index("Alpha alpha A1")])


def test_an_untagged_copy_of_a_dead_count_never_becomes_a_growth_curve_in_a_search(client):
    # pins the wiring of the dead-copy rule through reading a replicate (a review: the unit test alone let it go)
    rising = [(0, 1), (4, 50), (8, 300), (12, 2000), (24, 2400)]
    for n, rep in enumerate(("a1", "a2", "a3")):
        dead_id, copy_id = 9000 + 2 * n, 9001 + 2 * n
        client.contexts[dead_id] = client.contexts[copy_id] = rising
        client.bioreplicates[rep]["measurementContexts"] += [
            {"id": dead_id, "techniqueType": "fc", "techniqueUnits": "", "techniqueCellType": "Dead",
             "subject": {"type": "bioreplicate", "name": rep}},
            {"id": copy_id, "techniqueType": "fc", "techniqueUnits": "", "subject": {"type": "strain", "name": rep}}]
    r = run(client, report_rates=True, rate_window=3)
    assert r["rates"]["ncbi:1"]["technique"] == "fc (culture)"


def test_a_genus_of_inconclusive_taxa_says_so_once():
    from foodnet import derive
    taxa = {"t1": {"genus": "G", "name": "G one"}}
    cells = {("t1", "m", "exponential"): {"mean": 1.0, "sd": None, "n": 2, "n_experiments": 1, "values": [],
                                          "state": "inconclusive", "direction": None, "cautions": ["inconclusive"],
                                          "notes": [], "experiments": ["E1"], "studies": ["S1"], "media": ["m"],
                                          "start": 0.0, "end": 8.0, "initial": 5.0}}
    merged = derive.merge_genus_cells(cells, taxa)
    notes = next(iter(merged.values()))["notes"]
    assert notes == ["inconclusive, so not counted: G one"]


def test_a_studys_ungated_counts_take_the_gate_its_authors_named():
    # Karoline, 2026-10-23: "FC counts in study 2 are live counts"
    from foodnet import reading
    fc = {"times": [0, 4], "values": [1, 10], "gate": "", "technique": "fc"}
    od = {"times": [0, 4], "values": [0.1, 0.4], "gate": "", "technique": "od"}
    reading.name_gates([fc, od], [], reading.KNOWN_GATES["SMGDB00000002"])
    assert fc["gate"] == "live" and od["gate"] == ""


def test_a_known_gate_reaches_a_searchs_technique(client, monkeypatch):
    # pins the wiring of KNOWN_GATES through reading a replicate (a review: only name_gates was tested)
    from foodnet import reading
    monkeypatch.setitem(reading.KNOWN_GATES, "SMGDB00000001", {"fc": "live"})
    r = run(client, report_rates=True, rate_window=3)
    assert r["rates"]["ncbi:1"]["technique"] == "fc live (culture)"


def test_an_identical_gated_copy_outranks_a_known_gate():
    from foodnet import reading
    total = {"times": [0, 4], "values": [1, 10], "gate": "total", "technique": "fc"}
    copy = {"times": [0, 4], "values": [1, 10], "gate": "", "technique": "fc"}
    reading.name_gates([copy, total], [], {"fc": "live"})
    assert copy["gate"] == "total"


def test_biomass_counts_the_cultures_of_another_unit_it_leaves_out():
    import types
    cultures = [(i, types.SimpleNamespace(growth={"times": [0, 12], "values": [1, 100], "technique": tech, "unit": ""},
                                          taxon={"id": "t"}, experiment=f"E{i}"))
                for i, tech in enumerate(("fc", "fc", "od"))]
    rows = [{"culture": i, "phase": "exponential", "change": -1.0, "start": 0, "end": 12, "cautions": []}
            for i in range(3)]
    b = crm.biomass_changes(cultures, rows, "exponential")["t"]
    assert (b["unit"], b["n"], b["other_units_left_out"]) == ("fc", 2, 1)


def test_merged_genera_name_the_cultures_their_maintenance_leaves_out():
    many = {"glc": [{"unit": "fc live", "culture": "A1"}, {"unit": "fc live", "culture": "A2"},
                    {"unit": "fc", "culture": "B1"}]}
    assert crm.unit_skips({"G": many}) == [("B1", "no maintenance coefficient: its growth curve is in fc, and most "
                                                   "of the cultures of its taxon (or genus) are in fc live")]


def test_crm_json_carries_the_biomass_counts_of_the_crms_own_phase(client):
    # a review: they reached only other_phases, so R read NA for every taxon
    p = matrix.crm_payload(run(client, report_rates=True, phase="both"))
    a = p["taxa"].index("Alpha alpha A1")
    assert p["biomass_replicates"][a] == 3 and p["biomass_other_units_left_out"][a] == 0
    assert p["other_phases"]["stationary"]["biomass_replicates"][a] == 3


def test_merged_genera_pool_every_strains_cultures_before_choosing_the_unit(client, monkeypatch):
    # a review: a strain's culture whose compound had no entry in the strain's own unit never reached the genus
    from foodnet import search
    seen = {}
    real = search.crm.maintenance

    def spy(*a, **k):
        seen["raw"] = k.get("raw")
        return real(*a, **k)
    monkeypatch.setattr(search.crm, "maintenance", spy)
    run(client, report_rates=True, merge_genera=True)
    assert seen["raw"] is True


def test_biomass_start_is_every_culture_at_one_hour_averaged_as_logs():
    # Karoline, 2026-10-23: "fix biomass_start first": cultures whose windows start at 0 h and at 4 h are all read
    # at 0 h (the median start), not averaged as they stood at their own starts
    import types
    curve = {"times": [0, 4, 8, 12], "values": [3e5, 3e7, 3e8, 4e8], "technique": "fc", "unit": ""}
    cultures = [(i, types.SimpleNamespace(growth=dict(curve), taxon={"id": "t"}, experiment="E1")) for i in range(3)]
    rows = [{"culture": i, "phase": "exponential", "change": -1.0, "start": s, "end": 12, "cautions": []}
            for i, s in enumerate((0, 0, 4))]
    b = crm.biomass_changes(cultures, rows, "exponential")["t"]
    assert b["start_h"] == 0 and b["start"] == pytest.approx(3e5)
