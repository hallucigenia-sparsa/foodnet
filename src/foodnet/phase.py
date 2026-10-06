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

**Or earlier, where growth stops first** (Karoline, 2026-10-06, after two rounds of review). The 90% rule
ends growth late on two kinds of real curves: a plateau that keeps creeping up (E. coli LF82 in study
SMGDB00000009 stops growing by about 10 h, but rises another 20% by 168 h, so 90% of the final maximum is
reached at 84 to 108 h) and a late rise that is not growth (B. hydrogenotrophica in study SMGDB00000004, whose
qPCR count rises again at 30 and 48 h while OD falls: DNA from lysing cells). So the boundary is the earlier of
the 90% rule and the growth-rate rule (`rate_end`): the first sample after the fastest growth from which the
specific growth rate stays below a tenth of its maximum over two consecutive intervals. Two intervals, so one
noisy interval neither sets the maximum nor ends growth. Each rule fails by being late, in different ways, so
the earlier one was right on every curve checked; on a curve sampled too coarsely for two slow intervals the
rate rule finds nothing and the 90% rule decides. A smoothed 90% rule was tried first and dropped: it moved
E. coli later still and clipped a real one-sample peak. A boundary on fewer than three samples is `coarse`,
which the arcs carry as the caution `coarse_sampling`.

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

import math

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


RATE_FRACTION = 0.1          # growth has stopped when the specific rate stays below this share of its maximum


def rate_end(times, values, fraction: float = RATE_FRACTION) -> float | None:
    """The time growth stops by the growth rate, or None when the curve cannot say (fewer than four positive
    samples, or no two consecutive slow intervals after the fastest growth).

    The specific rate of each interval is ln(x1 / x0) / (t1 - t0). The maximum is the best rate held over two
    consecutive intervals, and growth stops at the start of the first two consecutive intervals after it whose
    rates are both below `fraction` of that maximum."""
    pts = [(float(t), float(v)) for t, v in zip(times, values, strict=True) if float(v) > 0]
    if len(pts) < 4:
        return None
    rates = [math.log(v1 / v0) / (t1 - t0) for (t0, v0), (t1, v1) in zip(pts, pts[1:], strict=False) if t1 > t0]
    if len(rates) < 3:
        return None
    held = [min(a, b) for a, b in zip(rates, rates[1:], strict=False)]
    top = max(held)
    if top <= 0:
        return None
    for i in range(held.index(top), len(rates) - 1):
        if rates[i] < fraction * top and rates[i + 1] < fraction * top:
            return pts[i][0]
    return None


def exponential_end(times, values, fraction: float = FRACTION, factor: float = NO_GROWTH_FACTOR) -> dict:
    """{"end": time, "index": i, "last": bool, "coarse": bool, "by": "90%" or "rate"}: where exponential
    growth ends on one growth curve: the earlier of the 90% rule and `rate_end`.

    `last` is True when the boundary is the curve's last point: the culture had not stopped growing when
    sampling ended, so there is no stationary phase in the data. `coarse` is True when fewer than three
    samples lead up to the boundary. Raises NoBoundary when the curve is too short, holds no positive value,
    or did not grow by `factor`.

    Worked example (the docstring of the test repeats it): times 0, 4, 8, 12, 16, 24 and abundances
    1, 10, 100, 500, 900, 1000. The start is 1 and the maximum 1000, so the 90% threshold is
    1 + 0.9 * (1000 - 1) = 900.1, first reached at 24 h (900 at 16 h is just below). The rates fall from
    0.58 per hour to 0.013 over the last interval, but only one interval is slow, so the rate rule finds no
    end, and the boundary is 24 h, the last point.
    """
    pairs = [(float(t), float(v)) for t, v in zip(times, values, strict=True)]
    if len(pairs) < 3:
        raise NoBoundary(f"{len(pairs)} growth time point(s); the end of growth needs at least 3")
    start = next((v for _, v in pairs if v > 0), None)
    if start is None:
        raise NoBoundary("no positive abundance in the growth curve")
    top = max(v for _, v in pairs)
    if top < factor * start:
        raise NoBoundary(f"did not grow: the maximum is {top / start:.2g} times the start, below {factor:g}")
    threshold = start + fraction * (top - start)
    end, by = next(t for t, v in pairs if v >= threshold), "90%"
    times = [t for t, _ in pairs]
    late = times.index(end)
    by_rate = rate_end(times, [v for _, v in pairs])
    if by_rate is not None and by_rate < end:
        end, by = by_rate, "rate"
    i = times.index(end)
    # "moved": the 90% rule alone would have ended growth more than one sample later (E. coli LF82: 8 h
    # instead of 84 h), which the arcs flag; one sample earlier is within the sampling's own resolution
    return {"end": end, "index": i, "last": i == len(pairs) - 1, "coarse": i < 2, "by": by,
            "moved": by == "rate" and late - i > 1}


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
