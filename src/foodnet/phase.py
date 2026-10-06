"""Where exponential growth ends, and what a metabolite did before and after.

Karoline, 2026-10-04: "metabolite production and consumption during exponential phase can differ from
stationary phase. So let's do this differently: estimate when the exponential phase ends and take that as
the threshold", with a choice of "Exponential phase" (the default), "Stationary phase" or "Both", and "This
will not treat diauxic shifts well, but we are also not able to identify them always clearly, so OK for now."

**The boundary.** The end of exponential growth is the first sampled time at which the culture has risen
`fraction` (0.9 by default, an advanced setting) of the way from its starting abundance to its maximum, on
the linear scale. Only measured time points are used; nothing is interpolated on the growth curve, so the
boundary is always a time at which the culture was sampled. The linear scale is the choice for a reason
found on real data (study SMGDB00000007, R. intestinalis): on a log scale the last doubling before the
plateau is a tenth of the rise, so a 90% rule on log abundance ended the phase at 12 h, while glucose was
still being taken up until 16 h, when the cells reached their maximum.

**The curve is smoothed first** (`smooth`, a running median of three with Tukey's end-point rule), and the
maximum and the crossing are read from the smoothed curve (Karoline, 2026-10-06, after a review found the
raw maximum on a noisy plateau: E. coli LF82 in study SMGDB00000009 stops growing by about 10 h, but counting
noise of 10% on its plateau put the boundary of single replicates at 84 to 108 h; and a late rise of one
point, B. hydrogenotrophica at 48 h in study SMGDB00000004, made the whole run "exponential"). A median of
three removes one stray point and leaves a real rise, which lasts more than one sample. The boundary also
says when it rests on fewer than three samples (`coarse`: a culture sampled every 24 h cannot place it more
finely), which the arcs carry as the caution `coarse_sampling`.

**A culture that did not grow** has no exponential phase. A rise below `NO_GROWTH_FACTOR` (1.5, grownet's
default for the same question) gives no boundary, and its metabolites give no value in the phases: the
matrices mark them `no_phase`, and a time window (which needs no boundary) gives them values.

**The change in a phase** is the metabolite concentration at the end of the phase minus the concentration
at its start, from the metabolite's own series. The phases are:

  * exponential: from the first metabolite sample to the boundary;
  * stationary: from the boundary to the last metabolite sample;
  * window (an advanced setting that overrides both): from a start to an end the user gives, in hours.

A metabolite is often sampled at other times than the cells, so its concentration at a phase boundary is
interpolated linearly between the two samples around it. A boundary outside the metabolite's samples is
never extrapolated: the nearest sample is used and the arc is marked `window_beyond_data`.
"""
from __future__ import annotations

import statistics

FRACTION = 0.9               # the share of the maximal abundance that ends exponential growth
NO_GROWTH_FACTOR = 1.5       # below this rise (maximum over start) the culture did not grow
SHORT_RECORD_H = 24.0        # a metabolite series shorter than this is flagged (Karoline, 2026-10-04)

# time units to hours
TIME_UNITS = {"h": 1.0, "hr": 1.0, "hrs": 1.0, "hour": 1.0, "hours": 1.0,
              "min": 1 / 60, "mins": 1 / 60, "minute": 1 / 60, "minutes": 1 / 60,
              "s": 1 / 3600, "sec": 1 / 3600, "seconds": 1 / 3600,
              "d": 24.0, "day": 24.0, "days": 24.0}


class NoBoundary(ValueError):
    """A growth curve gives no end of exponential growth, with the reason."""


def hours(unit: str) -> float | None:
    """The factor that turns a time in `unit` into hours, or None for a unit foodnet does not know."""
    return TIME_UNITS.get((unit or "").strip().lower())


def smooth(values) -> list:
    """A running median of three, with Tukey's end-point rule at the end: the last point becomes the median of
    itself, its smoothed neighbor and the value that neighbor's trend extrapolates to, so a curve still rising
    at its end still ends at its maximum and a lone jump at the end is taken down to its neighbors. The first
    point is the culture's start and stays as measured (the end-point rule there would lift a culture that
    had already grown by its second sample to its plateau). Curves of fewer than four points are returned as
    they are."""
    v = [float(x) for x in values]
    if len(v) < 4:
        return v
    out = [statistics.median(v[i - 1:i + 2]) for i in range(1, len(v) - 1)]
    last = statistics.median([v[-1], out[-1], 3 * out[-1] - 2 * out[-2]])
    return [v[0], *out, last]


def exponential_end(times, values, fraction: float = FRACTION, factor: float = NO_GROWTH_FACTOR) -> dict:
    """{"end": time, "index": i, "last": bool, "coarse": bool}: where exponential growth ends on one curve.

    `last` is True when the boundary is the curve's last point: the culture had not stopped growing when
    sampling ended, so there is no stationary phase in the data. `coarse` is True when fewer than three
    samples lead up to the boundary. Raises NoBoundary when the curve is too short, holds no positive value,
    or did not grow by `factor`.

    Worked example (the docstring of the test repeats it): times 0, 4, 8, 12, 16, 24 and abundances
    1, 10, 100, 500, 900, 1000. Smoothed, they are 1, 10, 100, 500, 900, 1000 (a rising curve keeps its
    shape). The start is 1 and the maximum 1000, so the threshold is 1 + 0.9 * (1000 - 1) = 900.1. The first
    sample at or above it is 24 h (900 at 16 h is just below), so the boundary is 24 h, and it is the last
    point.
    """
    pairs = [(float(t), float(v)) for t, v in zip(times, values, strict=True)]
    if len(pairs) < 3:
        raise NoBoundary(f"{len(pairs)} growth time point(s); the end of growth needs at least 3")
    start = next((v for _, v in pairs if v > 0), None)
    if start is None:
        raise NoBoundary("no positive abundance in the growth curve")
    smoothed = smooth([v for _, v in pairs])
    top = max(smoothed)
    if top < factor * start:
        raise NoBoundary(f"did not grow: the maximum is {top / start:.2g} times the start, below {factor:g}")
    threshold = start + fraction * (top - start)
    for i, ((t, _), v) in enumerate(zip(pairs, smoothed, strict=True)):
        if v >= threshold:
            return {"end": t, "index": i, "last": i == len(pairs) - 1, "coarse": i < 2}
    raise NoBoundary("no point reaches the threshold")       # not reached: the maximum always does


def value_at(series, t: float) -> tuple:
    """(concentration at time t, whether t lay outside the samples) from a sorted [(time, value)] series.

    Inside the samples the value is interpolated linearly between the two around t; outside, the nearest
    sample is returned and the second element is True."""
    if not series:
        raise ValueError("an empty series")
    if t <= series[0][0]:
        return series[0][1], t < series[0][0]
    if t >= series[-1][0]:
        return series[-1][1], t > series[-1][0]
    for (t0, v0), (t1, v1) in zip(series, series[1:], strict=False):
        if t0 <= t <= t1:
            if t1 == t0:
                return v1, False
            return v0 + (v1 - v0) * (t - t0) / (t1 - t0), False
    return series[-1][1], False                               # not reached


def change(series, start: float, end: float) -> dict:
    """{"change", "start", "end", "beyond"}: the net change of a metabolite from `start` to `end` (hours).

    `beyond` is True when `end` lies after the last sample (or `start` before the first), so the value at
    that end is the nearest sample's and not an interpolation."""
    first, before = value_at(series, start)
    last, after = value_at(series, end)
    return {"change": last - first, "start": start, "end": end, "beyond": after or before, "initial": first}


def phase_windows(series, boundary: dict | None, phase: str, window: tuple | None = None) -> dict:
    """{phase: (start, end, cautions)} for the phases asked for, from one metabolite series.

    `phase` is "exponential", "stationary" or "both"; `window`, when given, replaces them with one
    "window" phase. `boundary` is what `exponential_end` returned for the replicate's growth curve.
    """
    first, last = series[0][0], series[-1][0]
    out = {}
    if window is not None:
        start, end = window
        cautions = ["window_beyond_data"] if end > last or start < first else []
        out["window"] = (max(start, first), end, cautions)
        return out
    if boundary is None:
        return out
    end = boundary["end"]
    if phase in ("exponential", "both"):
        out["exponential"] = (first, end, ["window_beyond_data"] if end > last else [])
    if phase in ("stationary", "both"):
        if boundary.get("last"):
            out["stationary"] = (None, None, ["stationary_not_reached"])
        elif end >= last:
            out["stationary"] = (None, None, ["window_beyond_data"])
        else:
            out["stationary"] = (end, last, [])
    return out


SPIKE_FACTOR = 100.0   # grownet's default; an advanced setting (0 switches the check off)
SPIKE_RUN = 2


def spike(values, factor: float = SPIKE_FACTOR) -> float | None:
    """The ratio of an implausible spike in a growth curve, or None: grownet's rule, one or two consecutive
    interior points above both neighbors by more than `factor` (in study SMGDB00000004 a qPCR trace holds
    two points near 5e13 cells/mL between neighbors of 1e8 and 5e8). A spike would make its point the
    maximum and so move the end of exponential growth, which is why a spiked curve sets no boundary."""
    if not factor:
        return None
    found = None
    for i in range(1, len(values) - 1):
        for length in range(1, SPIKE_RUN + 1):
            j = i + length
            if j > len(values) - 1:
                break
            reference = max(values[i - 1], values[j])
            if reference <= 0:
                continue
            ratio = min(values[i:j]) / reference
            if ratio > factor and (found is None or ratio > found):
                found = ratio
    return found


def short_record(series) -> bool:
    """Whether a metabolite series covers less than SHORT_RECORD_H hours."""
    return (series[-1][0] - series[0][0]) < SHORT_RECORD_H if series else True
