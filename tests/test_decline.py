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
    assert m["start_mM"] == pytest.approx(2.0) and not m["depleted_at_start"]
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
    assert matrix.rates_csv(r).splitlines()[0].endswith(",technique")
    assert matrix.declines_csv(r).splitlines()[0].endswith(",technique")
    warning = matrix.technique_mismatch(r)
    assert warning.startswith("Growth curves of different techniques: fc") and "od for" in warning
    assert warning in matrix.readme(r, "crm")


def test_the_decline_starts_where_the_curve_enters_its_plateau_not_at_its_highest_sample():
    # Karoline, 2026-10-22: "Anchor at plateau entry"; a noisy high point at 24 h no longer starts the slope
    d = rates.decline([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 1300, 1000, 1000])
    assert d["start"] == 12 and d["points"] == 4
    with pytest.raises(rates.RateUnavailable, match="falls to zero at 36 h"):
        rates.decline([0, 4, 8, 12, 24, 36, 48], [1, 10, 100, 1000, 100, 0, 0])


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
    assert any("in od" in why for _, why in skipped)
    used_up, _ = crm.maintenance([_plateau_culture("E1", (0.05, 0.0)), _plateau_culture("E1", (0.1, 0.0))], 0.2)
    gone = used_up["ncbi:9"]["glc"]
    assert gone["coefficient"] == 0 and gone["depleted_at_start"] and gone["start_mM"] == pytest.approx(0.075)
