"""Where exponential growth ends, and the change of a metabolite over a phase (foodnet.phase)."""
import pytest

from foodnet import phase


def test_the_boundary_is_the_first_sample_at_ninety_percent_of_the_maximum():
    # the docstring's example: start 1, maximum 1000, threshold 1 + 0.9 * 999 = 900.1; 900 at 16 h is just
    # below it, so the first sample at or above is 24 h, which is also the last point
    b = phase.exponential_end([0, 4, 8, 12, 16, 24], [1, 10, 100, 500, 900, 1000])
    assert b == {"end": 24.0, "index": 5, "last": True, "coarse": False}


def test_the_boundary_ignores_the_decline_after_the_maximum():
    # the threshold is 1 + 0.9 * 999 = 900.1; 950 at 12 h is the first sample above it
    b = phase.exponential_end([0, 4, 8, 12, 16, 24], [1, 10, 100, 950, 1000, 200])
    assert b["end"] == 12.0 and not b["last"]


def test_the_fraction_is_a_setting():
    # with 0.5: threshold 1 + 0.5 * 999 = 500.5, first reached at 12 h (950)
    assert phase.exponential_end([0, 4, 8, 12], [1, 10, 100, 950], fraction=0.5)["end"] == 12.0
    # with 0.05: threshold 1 + 0.05 * 949 = 48.45, first reached at 8 h (100)
    assert phase.exponential_end([0, 4, 8, 12], [1, 10, 100, 950], fraction=0.05)["end"] == 8.0


def test_a_culture_that_did_not_grow_has_no_boundary():
    # 1.4 times its start, below the no-growth factor 1.5
    with pytest.raises(phase.NoBoundary, match="did not grow"):
        phase.exponential_end([0, 4, 8], [1.0, 1.2, 1.4])


def test_too_short_a_curve_has_no_boundary():
    with pytest.raises(phase.NoBoundary, match="at least 3"):
        phase.exponential_end([0, 4], [1, 100])


def test_a_concentration_between_samples_is_interpolated():
    series = [(0, 10.0), (8, 6.0), (16, 2.0)]
    assert phase.value_at(series, 4) == (8.0, False)       # halfway from 10 to 6
    assert phase.value_at(series, 12) == (4.0, False)
    assert phase.value_at(series, 16) == (2.0, False)


def test_a_time_outside_the_samples_takes_the_nearest_sample_and_says_so():
    series = [(2, 10.0), (8, 6.0)]
    assert phase.value_at(series, 20) == (6.0, True)
    assert phase.value_at(series, 0) == (10.0, True)


def test_the_change_is_end_minus_start():
    series = [(0, 10.0), (8, 6.0), (16, 2.0)]
    d = phase.change(series, 0, 12)
    assert d["change"] == -6.0 and not d["beyond"] and d["initial"] == 10.0
    assert phase.change(series, 0, 30)["beyond"]


def test_phases_split_at_the_boundary():
    series = [(0, 10.0), (12, 2.0), (24, 1.0)]
    windows = phase.phase_windows(series, {"end": 12.0, "last": False}, "both")
    assert windows == {"exponential": (0.0, 12.0, []), "stationary": (12.0, 24.0, [])}


def test_no_stationary_phase_when_the_culture_was_still_growing():
    series = [(0, 10.0), (12, 2.0), (24, 1.0)]
    windows = phase.phase_windows(series, {"end": 24.0, "last": True}, "stationary")
    assert windows == {"stationary": (None, None, ["stationary_not_reached"])}


def test_a_window_replaces_the_phases_and_flags_an_end_after_the_data():
    series = [(0, 10.0), (12, 2.0), (24, 1.0)]
    assert phase.phase_windows(series, None, "exponential", (0, 48)) == {"window": (0, 48, ["window_beyond_data"])}
    assert phase.phase_windows(series, None, "exponential", (6, 24)) == {"window": (6, 24, [])}


def test_a_record_shorter_than_24_hours_is_short():
    assert phase.short_record([(0, 1), (20, 2)])
    assert not phase.short_record([(0, 1), (24, 2)])


def test_hours_from_the_units_mgrowthdb_records():
    assert phase.hours("h") == 1.0 and phase.hours("min") == pytest.approx(1 / 60) and phase.hours("d") == 24.0
    assert phase.hours("fortnights") is None


def test_a_spiked_growth_curve_is_found():
    # 1e5 between neighbors of 10 and 30: 3333 times the larger neighbor
    assert phase.spike([1, 10, 1e5, 30, 40]) == pytest.approx(1e5 / 30)
    assert phase.spike([1, 10, 100, 1000]) is None
