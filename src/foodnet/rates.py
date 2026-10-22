"""Growth rates from a growth curve, two ways, standard library only (grownet #41).

Karoline decided on 2026-09-19 (on grownet #41) that the growth rate ships as a metric with two implementations,
and on 2026-09-27 that the implementation is its own option:

* `easylinear` (the default): the maximum specific growth rate by the method of Hall et al. (2014) as the
  growthrates R package implements it in `fit_easylinear`, which mGrowthDB uses for the rates it reports.
  Fit a straight line to log abundance over every window of `window` consecutive points, take the steepest,
  widen it to span every window whose slope reaches `QUOTA` of that maximum, and fit once more over the
  widened range. With the default window of 5 it reproduced mGrowthDB's reported rates on SMGDB00000004
  (median ratio 1.00, 43 of 46 within 10%), so a network on growth rate agrees with the database.
* `baranyi`: the Baranyi-Roberts model fitted by Levenberg-Marquardt, with a multi-start and two guards,
  because an unguarded fit converged on every curve, some to nonsense: mu must stay within 5 times the
  steepest observed slope, and the fit must explain at least R2 = 0.9 of log abundance. A rejected fit is
  reported, never replaced by another number. It is a different quantity from easylinear (median 1.16
  times the reported rate on SMGDB00000004). It is fitted to the part of the curve the model describes: up
  to the maximum and on along the plateau after it, cut where the curve declines below the plateau
  (Karoline, 2026-09-28: "Fit up to the maximum", then "End of the plateau"). Fitted to whole curves it was
  rejected on 152 curves of mGrowthDB that decline after their peak or grow in two phases.

A rate is per time unit of the curve; the comparison takes log2 of rates with and without the partner,
like any other metric. `RateUnavailable` says why a curve has no rate, and the caller reports it.
"""
from __future__ import annotations

import math

METHODS = ("easylinear", "baranyi")
DEFAULT_METHOD = "easylinear"
DEFAULT_WINDOW = 5          # growthrates' default h, which matched mGrowthDB's reported rates
QUOTA = 0.95                # growthrates' default quota
BARANYI_MIN_POINTS = 6
BARANYI_MIN_R2 = 0.9
PLATEAU = 0.10              # the plateau after the maximum: within this share of the rise below it (log)


def _until_decline(xs, ys):
    """The points up to the end of the plateau after the curve's maximum: the plateau goes on while the log
    abundance stays within PLATEAU times the rise (maximum minus first value) below the maximum, and ends at
    the first point below that, where the decline begins."""
    if not ys:
        return xs, ys
    top = max(ys)
    band = top - PLATEAU * (top - ys[0])
    end = ys.index(top)
    while end + 1 < len(ys) and ys[end + 1] >= band:
        end += 1
    return xs[:end + 1], ys[:end + 1]


_MAX_EXP = 700.0


class RateUnavailable(ValueError):
    """A curve has no growth rate by this method, with the reason."""


def _positive_logs(times, values):
    pairs = [(t, math.log(v)) for t, v in zip(times, values, strict=True) if v > 0]
    return [t for t, _ in pairs], [y for _, y in pairs]


def _line(xs, ys):
    """Least-squares slope and intercept, or None when the x values do not vary."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / sxx
    return slope, my - slope * mx


def easylinear(times, values, window: int = DEFAULT_WINDOW, quota: float = QUOTA) -> float:
    """The maximum specific growth rate by growthrates' fit_easylinear (see the module docstring)."""
    return easylinear_fit(times, values, window, quota)[0]


def easylinear_fit(times, values, window: int = DEFAULT_WINDOW, quota: float = QUOTA) -> tuple:
    """(rate, first hour, last hour of the stretch the final line is fitted over), by fit_easylinear; the hours
    are None when no stretch rises."""
    if window < 2:
        raise ValueError(f"the growth rate window must hold at least 2 points, not {window}")
    xs, ys = _positive_logs(times, values)
    if len(xs) < window + 1:
        raise RateUnavailable(f"{len(xs)} positive time point(s); easylinear with a window of {window} needs "
                              f"{window + 1}")
    # growthrates fits the windows starting at points 1 to N - h, so the last window it tries starts one
    # point before the last one that fits; kept, so the numbers agree with mGrowthDB's
    slopes = []
    for i in range(len(xs) - window):
        fit = _line(xs[i:i + window], ys[i:i + window])
        slopes.append(fit[0] if fit else float("-inf"))
    best = max(slopes)
    if best == float("-inf"):
        raise RateUnavailable("no window with varying times")
    if best <= 0:
        # no window rises: the steepest slope is itself not positive, and the quota below would select no
        # window at all (0.95 of a negative slope lies above it). Returned as it is, so the caller treats it
        # as any non-positive property (code review of 2026-09-28: it raised a bare ValueError that dropped
        # the whole comparison, SMGDB00000014)
        return best, None, None
    candidates = [i for i, s in enumerate(slopes) if s >= quota * best]
    first, last = min(candidates), max(candidates) + window
    fit = _line(xs[first:last], ys[first:last])
    return fit[0], xs[first], xs[last - 1]


# ---- Baranyi-Roberts --------------------------------------------------------------------------------

def _baranyi(t, y0, mu, d, h0):
    """ln abundance at time t; d = ymax - y0 (positive), h0 = mu * lag. Rearranged so nothing overflows."""
    inner = 1.0 - (1.0 - math.exp(-min(mu * t, _MAX_EXP))) * (1.0 - math.exp(-min(h0, _MAX_EXP)))
    a = t + math.log(max(inner, 1e-300)) / mu
    q = mu * a
    u = q - d
    if u > _MAX_EXP:
        term = u + math.log1p((1.0 - math.exp(-d)) * math.exp(-u))
    else:
        term = math.log(max(1.0 - math.exp(-d) + math.exp(u), 1e-300))
    return y0 + q - term


def _predict(times, params):
    y0, log_mu, log_d, log_h0 = params
    mu = max(math.exp(min(log_mu, _MAX_EXP)), 1e-9)
    d = max(math.exp(min(log_d, _MAX_EXP)), 1e-9)
    h0 = max(math.exp(min(log_h0, _MAX_EXP)), 1e-12)
    return [_baranyi(t, y0, mu, d, h0) for t in times]


def _steepest(times, logs, span=2):
    best = (0.0, times[0], logs[0])
    for i in range(len(times) - span):
        dt = times[i + span] - times[i]
        if dt > 0 and (logs[i + span] - logs[i]) / dt > best[0]:
            best = ((logs[i + span] - logs[i]) / dt, times[i], logs[i])
    return best


def _solve(matrix, rhs):
    """Gaussian elimination with partial pivoting; None when singular."""
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-14:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        for row in range(col + 1, n):
            f = a[row][col] / a[col][col]
            for k in range(col, n + 1):
                a[row][k] -= f * a[col][k]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        x[row] = (a[row][n] - sum(a[row][k] * x[k] for k in range(row + 1, n))) / a[row][row]
    return x


def _fit_once(times, logs, params, max_iter=200):
    def rss(p):
        try:
            return sum((m - o) ** 2 for m, o in zip(_predict(times, p), logs, strict=True))
        except (ValueError, OverflowError):
            return float("inf")

    params, lam = params[:], 1e-3
    current = rss(params)
    for _ in range(max_iter):
        residuals = [m - o for m, o in zip(_predict(times, params), logs, strict=True)]
        jac = [[0.0] * 4 for _ in times]
        for j in range(4):
            step = 1e-5 * max(abs(params[j]), 1.0)
            up, down = params[:], params[:]
            up[j] += step
            down[j] -= step
            try:
                f_up, f_down = _predict(times, up), _predict(times, down)
            except (ValueError, OverflowError):
                return None
            for i in range(len(times)):
                jac[i][j] = (f_up[i] - f_down[i]) / (2 * step)
        jtj = [[sum(r[a] * r[b] for r in jac) for b in range(4)] for a in range(4)]
        jtr = [sum(jac[i][a] * residuals[i] for i in range(len(times))) for a in range(4)]
        for _ in range(30):
            damped = [[jtj[a][b] + (lam * jtj[a][a] if a == b else 0.0) for b in range(4)] for a in range(4)]
            delta = _solve(damped, [-v for v in jtr])
            if delta is None:
                lam *= 10
                continue
            candidate = [params[j] + delta[j] for j in range(4)]
            trial = rss(candidate)
            if trial < current:
                improvement = current - trial
                params, current, lam = candidate, trial, max(lam / 10, 1e-12)
                if improvement < 1e-12 * max(current, 1e-12):
                    return params, current
                break
            lam *= 10
        else:
            return params, current
    return params, current


def baranyi(times, values) -> float:
    """The Baranyi-Roberts maximum specific growth rate, guarded, fitted up to the end of the plateau after the
    curve's maximum (see the module docstring)."""
    return baranyi_fit(times, values)[0]


def baranyi_fit(times, values) -> tuple:
    """(rate, first hour, last hour the model was fitted over: the first sample to the end of the plateau)."""
    xs, ys = _positive_logs(times, values)
    xs, ys = _until_decline(xs, ys)
    if len(xs) < BARANYI_MIN_POINTS:
        raise RateUnavailable(f"{len(xs)} positive time point(s) up to the end of the plateau; a Baranyi fit needs "
                              f"{BARANYI_MIN_POINTS}")
    slope, t_at, y_at = _steepest(xs, ys)
    y0 = sum(ys[:2]) / 2
    mu0 = max(slope, 1e-3)
    lag = max(t_at - (y_at - y0) / mu0, 0.0)
    base = [y0, math.log(mu0), math.log(max(max(ys) - y0, 0.1)), math.log(max(mu0 * lag, 1e-3))]
    starts = [base]
    for mu_factor in (0.5, 2.0):
        for lag_factor in (0.25, 4.0):
            start = base[:]
            start[1] = math.log(max(math.exp(base[1]) * mu_factor, 1e-3))
            start[3] = math.log(max(math.exp(base[3]) * lag_factor, 1e-3))
            starts.append(start)
    best = None
    for start in starts:
        result = _fit_once(xs, ys, start)
        if result is None:
            continue
        params, residual = result
        mu = math.exp(params[1])
        if not 1e-3 < mu < 5 * max(slope, 1e-3):
            continue                        # an optimum far from the steepest observed slope is not credible
        if best is None or residual < best[1]:
            best = (params, residual)
    if best is None:
        raise RateUnavailable("Baranyi fit rejected: no start converged to a rate near the observed slope")
    mean = sum(ys) / len(ys)
    total = sum((y - mean) ** 2 for y in ys)
    r2 = 1 - best[1] / total if total > 0 else float("nan")
    if not r2 >= BARANYI_MIN_R2:
        raise RateUnavailable(f"Baranyi fit rejected: R2 {r2:.2f} below {BARANYI_MIN_R2}; the model does not "
                              "describe this curve up to the end of its plateau (for example two growth phases)")
    return math.exp(best[0][1]), xs[0], xs[-1]


DECLINE_MIN_POINTS = 3      # the peak and two samples after it: a slope, not a step


def decline(times, values, fraction: float = 0.9) -> dict:
    """The sustained decline once a growth curve has stopped growing (Karoline, 2026-10-09, for the biomass loss
    rate of a consumer-resource fit): the least-squares slope of log abundance from the first sample at or above
    `fraction` of the way from the first sample to the maximum (foodnet's 90% rule for the end of exponential
    growth, the `fraction` setting, without its growth-rate rule) to the last sample, as a rate of loss per hour
    (positive: the curve falls; negative: it still creeps up), with that window (Karoline, 2026-10-22: "First
    sample >= 90% of max"; reviews: a slope
    from the highest sample started on noise, and one from the plateau band's lower edge began inside
    exponential growth). `death` says whether the curve leaves its plateau band (PLATEAU, the cut the Baranyi
    fit uses) before the last sample. A curve that falls to zero is fitted over its positive samples up to the
    zero, marked death and `to_zero`, when that leaves DECLINE_MIN_POINTS of them: the fall to zero itself is not
    in the slope, so it bounds the decline from below. What a curve counts decides
    what the rate is: optical density and total cell counts keep dead cells until they lyse, so for anything
    but a viable count it is a loss of counted biomass, which bounds death from below. RateUnavailable when
    there is no usable decline."""
    pairs = sorted(zip(times, values, strict=True))
    positive = [(t, v) for t, v in pairs if v > 0]
    if not positive:
        raise RateUnavailable("no positive values")
    first, top = positive[0][1], max(v for _, v in positive)
    start = next(t for t, v in positive if v >= first + fraction * (top - first))
    # a fall to zero is one the curve never recovers from; a single zero followed by counts again is a sample
    # below detection, left out, not an end (a review: one dropout cut a plateau's fit short and read as death)
    zero = next((t for t, v in pairs if t > start and v <= 0 and all(w <= 0 for u, w in pairs if u >= t)), None)
    tail = [(t, math.log(v)) for t, v in positive if t >= start and (zero is None or t < zero)]
    if len(tail) < DECLINE_MIN_POINTS:
        raise RateUnavailable(f"fewer than {DECLINE_MIN_POINTS} samples from the end of growth to the "
                              + ("zero" if zero is not None else "last sample"))
    fit = _line([t for t, _ in tail], [y for _, y in tail])
    if fit is None:
        raise RateUnavailable("the samples after the end of growth share one time")
    xs, ys = _positive_logs(times, values)
    plateau_x, _ = _until_decline(xs, ys)
    return {"rate": -fit[0] + 0.0, "start": tail[0][0], "end": tail[-1][0], "points": len(tail),
            "death": zero is not None or len(plateau_x) < len(xs), "to_zero": zero is not None}


def plateau(times, values) -> tuple:
    """(start, end) of a growth curve's plateau: from its peak (the highest sample) to the last sample still
    within the plateau band (PLATEAU, the rule `decline` and the Baranyi fit use). RateUnavailable when the
    curve has no sample on its plateau after the peak."""
    xs, ys = _positive_logs(times, values)
    if not ys:
        raise RateUnavailable("no positive values")
    peak = ys.index(max(ys))
    kept_x, _ = _until_decline(xs, ys)
    if len(kept_x) - 1 <= peak:
        raise RateUnavailable("no sample on the plateau after the peak")
    return xs[peak], kept_x[-1]


def method_name(rate_method: str = DEFAULT_METHOD, window: int = DEFAULT_WINDOW) -> str:
    """The metric name an edge records: "growth_rate:easylinear:5" or "growth_rate:baranyi"."""
    if rate_method not in METHODS:
        raise ValueError(f"unknown growth rate method {rate_method!r}; choose one of {METHODS}")
    return f"growth_rate:easylinear:{int(window)}" if rate_method == "easylinear" else "growth_rate:baranyi"


def fitted(name: str):
    """The (times, values) -> (rate, first hour, last hour fitted over) function a metric name stands for."""
    parts = name.split(":")
    kind = parts[1] if len(parts) > 1 else DEFAULT_METHOD
    if parts[0] == "growth_rate" and kind == "baranyi" and len(parts) == 2:
        return baranyi_fit
    if parts[0] == "growth_rate" and kind == "easylinear" and len(parts) <= 3:
        window = int(parts[2]) if len(parts) == 3 else DEFAULT_WINDOW
        return lambda times, values: easylinear_fit(times, values, window)
    raise ValueError(f"unknown growth rate metric {name!r}")


def feature(name: str):
    """The (times, values) -> rate function a metric name stands for, or None if it is not a rate."""
    parts = name.split(":")
    if parts[0] != "growth_rate":
        return None
    kind = parts[1] if len(parts) > 1 else DEFAULT_METHOD
    if kind == "baranyi" and len(parts) == 2:
        return baranyi
    if kind == "easylinear" and len(parts) <= 3:
        window = int(parts[2]) if len(parts) == 3 else DEFAULT_WINDOW
        return lambda times, values: easylinear(times, values, window)
    raise ValueError(f"unknown growth rate metric {name!r}")
