"""
trendfit.py -- fit a straight line to the measured Amp:Kan trend.

WHAT IS BEING FITTED
--------------------
The ratio panel is a log axis, so a straight line there is an EXPONENTIAL
decay in the ratio itself:

    log10(Amp:Kan) = intercept + slope * round

The reported "fold change per round" is 10**slope. A slope of -1 means the
ratio falls 10x per round. This is a description of the observed trend, not a
mechanism -- the model is what carries the mechanism, and the model's own
curve is NOT what is being fitted here.

Numpy only, like the rest of the model. Weighted least squares via lstsq, and
a small table of t critical values in place of scipy.stats.t.

WHAT THE DATA CAN AND CANNOT SUPPORT -- READ BEFORE QUOTING A SLOPE
-------------------------------------------------------------------
Four limits, in descending order of how much they should worry you:

1. THE POINTS ARE NOT REPLICATES. They span seven constructs across two
   experiments. A single line through them treats construct-to-construct and
   run-to-run differences as random scatter about one common decay. They are
   not: at Phase 2 rd 1 the values run 0.70 (HP17KO7, run A) down to 0.068
   (Strep3, run B), a 10x spread WITHIN one round. The fitted slope is
   therefore an average over heterogeneous things, and its confidence
   interval understates the real uncertainty because it treats that spread as
   noise. This is the dominant caveat and no amount of extra statistics fixes
   it -- only more rounds per lineage would.

2. THREE ROUNDS, AND THE LAST ONE RESTS ON ONE POINT. Uncensored points sit
   at rounds 0, 1, 2 with counts 2, 4, 1. A two-parameter line over three
   x-levels has little leverage, and the slope is set almost entirely by the
   Jumpstart-to-round-1 drop. Drop the single round-2 point and the fit
   barely moves, which is a sign of weak constraint, not robustness.

3. SIX OF THIRTEEN POINTS ARE CENSORED. One-sided bounds ("> 63", "< 0.016")
   carry real information but cannot enter a least-squares fit as if they
   were measurements. Treating a "< 0.016" as the number 0.016 would bias the
   slope, so they are EXCLUDED from the fit and then checked against it
   (`censored_consistency`). A fit that contradicts a bound is a failed fit
   however good its R^2 looks.

4. THE JUMPSTART IS A DIFFERENT EXPERIMENT. Round 0 is the Jumpstart, which
   uses a different host and protocol from the Phase 2 rounds. Including it
   assumes one continuous process across a protocol change. `fit(skip_jumpstart
   =True)` drops it -- but then only two x-levels remain and the slope has no
   residual degrees of freedom at all.

So: the fit below is a legitimate summary of the observed downward trend. It
is not an estimate of any single lineage's decay rate.
"""

import math

import numpy as np

from labdata import markers

# Two-sided 95% t critical values by residual degrees of freedom. scipy is not
# installed beside this model, and these are the only quantiles needed.
_T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
        7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179,
        15: 2.131, 20: 2.086, 30: 2.042}


def _t95(df):
    """Two-sided 95% t value, nearest tabulated df, falling back to normal."""
    if df <= 0:
        return float("nan")
    if df in _T95:
        return _T95[df]
    smaller = [k for k in _T95 if k <= df]
    return _T95[max(smaller)] if smaller else 1.96


def _sigma_log10(m):
    """
    Per-point standard deviation in log10 units, from the sheet's interval.

    The intervals are roughly symmetric in LOG space (they come from a ratio),
    so half-width / 1.96 is a reasonable sd there. Returns None when the sheet
    gives no interval, in which case the point is fitted unweighted.
    """
    if m["lo"] is None or m["hi"] is None or m["lo"] <= 0:
        return None
    half = (np.log10(m["hi"]) - np.log10(m["lo"])) / 2.0
    return half / 1.96 if half > 0 else None


def fit(weighted=True, skip_jumpstart=False):
    """
    Weighted log-linear fit to the UNCENSORED points.

    weighted        use the sheet's intervals as 1/variance weights. Points
                    with no interval get the median weight rather than being
                    dropped or treated as infinitely precise.
    skip_jumpstart  exclude round 0 (see caveat 4 in the module docstring).

    Returns a dict, or None when fewer than three points remain (a line
    through two points is exact and says nothing about uncertainty).
    """
    pts = [m for m in markers() if not m["bound"] and m["value"] > 0]
    if skip_jumpstart:
        pts = [m for m in pts if m["index"] > 0]
    if len(pts) < 3:
        return None

    x = np.array([m["index"] for m in pts], float)
    y = np.log10(np.array([m["value"] for m in pts], float))

    if weighted:
        sig = [_sigma_log10(m) for m in pts]
        known = [s for s in sig if s is not None]
        fill = float(np.median(known)) if known else 1.0
        sig = np.array([fill if s is None else s for s in sig], float)
        w = 1.0 / sig                      # lstsq works on sqrt-weights
    else:
        w = np.ones_like(y)

    # Design matrix [round, 1]; solve the weighted normal equations.
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
    slope, intercept = float(coef[0]), float(coef[1])

    resid = y - (slope * x + intercept)
    df = len(pts) - 2
    # Weighted residual variance, scaled so the covariance is in y units.
    chi2 = float(np.sum((resid * w) ** 2))
    s2 = chi2 / df if df > 0 else float("nan")
    cov = np.linalg.inv((A * w[:, None]).T @ (A * w[:, None])) * s2
    se_slope = float(np.sqrt(cov[0, 0]))
    t = _t95(df)

    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_res = float(np.sum(resid ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    return {
        "slope": slope,                      # log10 ratio per round
        "intercept": intercept,
        "se_slope": se_slope,
        "slope_lo": slope - t * se_slope,
        "slope_hi": slope + t * se_slope,
        "fold_per_round": 10.0 ** slope,     # multiply the ratio each round
        "fold_lo": 10.0 ** (slope - t * se_slope),
        "fold_hi": 10.0 ** (slope + t * se_slope),
        "r2": r2,
        "n": len(pts),
        "df": df,
        "x_levels": sorted(set(int(v) for v in x)),
        "weighted": weighted,
        "skip_jumpstart": skip_jumpstart,
        "points": pts,
    }


def predict(f, index):
    """The fitted ratio at a round index."""
    return 10.0 ** (f["slope"] * index + f["intercept"])


def censored_consistency(f):
    """
    Check each censored bound against the fitted line.

    A "< v" row says the true ratio is BELOW v, so a fit predicting more than
    v contradicts it. A "> v" row says the opposite. Returns a list of
    (construct, run, index, kind, bound, fitted, ok).

    These rows were excluded from the fit, so this is a genuine out-of-sample
    check -- the one real test available here.
    """
    out = []
    for m in markers():
        if not m["bound"]:
            continue
        if f["skip_jumpstart"] and m["index"] == 0:
            continue
        yhat = predict(f, m["index"])
        ok = (yhat >= m["value"]) if m["bound"] == "lower" else (yhat <= m["value"])
        out.append((m["construct"], m["run"], m["index"], m["bound"],
                    m["value"], yhat, ok))
    return out


# ---------------------------------------------------------------------------
# Censored (Tobit) fit -- the statistically correct one here.
#
# WHY THE PLAIN FIT ABOVE IS NOT ENOUGH: the censoring is INFORMATIVE. A plate
# is reported as "< v" precisely because its Amp count collapsed to zero --
# that is, because the ratio was LOW. Dropping those rows therefore removes
# the smallest values systematically, not at random, so a fit to the survivors
# is biased upward by construction. `censored_consistency` on the plain fit
# shows this plainly: every held-out bound is contradicted.
#
# The fix is to keep the bounds in the likelihood as inequalities. A "< v" row
# contributes the probability that the true value lies below log10(v); a "> v"
# row the probability it lies above. Uncensored rows contribute the usual
# normal density. Maximising that total is a Tobit regression.
#
# Fitted by coarse-to-fine grid search over (slope, intercept, log10 sigma).
# Three parameters make this cheap and, unlike a hand-rolled Nelder-Mead, it
# cannot quietly converge to a local optimum or fail to converge at all.
# ---------------------------------------------------------------------------

def _norm_cdf(z):
    """Standard normal CDF via erf -- scipy is not available."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _log_lik(slope, intercept, log_sigma, rows):
    """Total log-likelihood of the censored model. -inf on impossible params."""
    sigma = 10.0 ** log_sigma
    if sigma <= 0:
        return -np.inf
    total = 0.0
    tiny = 1e-300
    for idx, yval, bound in rows:
        mu = slope * idx + intercept
        z = (yval - mu) / sigma
        if bound is None:
            # Normal density in log10 space.
            total += (-0.5 * z * z - math.log(sigma)
                      - 0.5 * math.log(2.0 * math.pi))
        elif bound == "upper":          # true value is BELOW yval
            total += math.log(max(_norm_cdf(z), tiny))
        else:                            # 'lower': true value is ABOVE yval
            total += math.log(max(1.0 - _norm_cdf(z), tiny))
    return total


def fit_censored(skip_jumpstart=False, refine=4):
    """
    Tobit fit using ALL 13 points -- bounds included as inequalities.

    Returns a dict shaped like fit()'s, plus 'log_lik'. The slope CI comes
    from a likelihood-ratio walk (the set of slopes within 1.92 log-likelihood
    units of the maximum, which is the standard 95% profile interval), not
    from a standard error -- the likelihood here is not symmetric.
    """
    rows = []
    for m in markers():
        if skip_jumpstart and m["index"] == 0:
            continue
        if m["value"] <= 0:
            continue
        rows.append((m["index"], math.log10(m["value"]), m["bound"]))
    if len(rows) < 3:
        return None

    lo = [-3.0, -1.5, -2.0]          # slope, intercept, log10 sigma
    hi = [1.0, 3.0, 1.0]
    best = None
    for step in range(refine):
        grids = [np.linspace(lo[i], hi[i], 41) for i in range(3)]
        best = None
        for s in grids[0]:
            for b in grids[1]:
                for ls in grids[2]:
                    ll = _log_lik(s, b, ls, rows)
                    if best is None or ll > best[0]:
                        best = (ll, s, b, ls)
        # Shrink the box around the current optimum and search again.
        span = [(hi[i] - lo[i]) / 8.0 for i in range(3)]
        centre = [best[1], best[2], best[3]]
        lo = [centre[i] - span[i] for i in range(3)]
        hi = [centre[i] + span[i] for i in range(3)]

    ll_max, slope, intercept, log_sigma = best

    # Profile the slope: widen outward until the likelihood drops by 1.92,
    # re-optimising the other parameters at each trial slope.
    def profile_at(s_try):
        sub = None
        for b in np.linspace(intercept - 2.0, intercept + 2.0, 81):
            for ls in np.linspace(log_sigma - 1.0, log_sigma + 1.0, 41):
                ll = _log_lik(s_try, b, ls, rows)
                if sub is None or ll > sub:
                    sub = ll
        return sub

    def find_edge(direction):
        s_edge = slope
        for _ in range(60):
            s_try = s_edge + direction * 0.02
            if profile_at(s_try) < ll_max - 1.92:
                return s_try
            s_edge = s_try
        return s_edge

    slope_lo, slope_hi = find_edge(-1), find_edge(+1)

    return {
        "slope": slope,
        "intercept": intercept,
        "sigma": 10.0 ** log_sigma,
        "se_slope": float("nan"),            # profile CI instead, see above
        "slope_lo": slope_lo,
        "slope_hi": slope_hi,
        "fold_per_round": 10.0 ** slope,
        "fold_lo": 10.0 ** slope_lo,
        "fold_hi": 10.0 ** slope_hi,
        "r2": float("nan"),                  # undefined with censored terms
        "n": len(rows),
        "n_censored": sum(1 for r in rows if r[2] is not None),
        "df": len(rows) - 3,
        "x_levels": sorted(set(int(r[0]) for r in rows)),
        "log_lik": ll_max,
        "weighted": False,
        "skip_jumpstart": skip_jumpstart,
        "censored_fit": True,
        "points": [m for m in markers()
                   if not (skip_jumpstart and m["index"] == 0)],
    }


def lineage_coverage():
    """
    Rounds measured per (construct, run) -- the only unit a decay rate would
    actually mean something for.

    Returns a list of (construct, run, n_rounds, n_uncensored, rounds), worst
    coverage last. A line needs at least three rounds from ONE lineage before
    its slope describes that lineage rather than an average over several.
    """
    from collections import defaultdict
    g = defaultdict(list)
    for m in markers():
        g[(m["construct"], m["run"])].append(m)
    out = []
    for (con, run), ms in g.items():
        rounds = sorted(set(m["index"] for m in ms))
        n_unc = sum(1 for m in ms if not m["bound"])
        out.append((con, run, len(rounds), n_unc, rounds))
    return sorted(out, key=lambda t: (-t[2], -t[3]))


def report_coverage():
    """Print per-lineage coverage and say plainly what it does not support."""
    print("  rounds measured per construct+run (the unit a decay rate would "
          "describe):")
    print(f"    {'construct':<16} {'run':<4} {'rounds':>6} {'uncensored':>11}"
          f"   which rounds")
    print("    " + "-" * 58)
    for con, run, n_r, n_u, rounds in lineage_coverage():
        print(f"    {con:<16} {run:<4} {n_r:>6} {n_u:>11}   {rounds}")
    best = max((t[3] for t in lineage_coverage()), default=0)
    print(f"\n  Most uncensored points from any single lineage: {best}.")
    print("  A two-parameter line needs 3+ to have any residual degrees of")
    print("  freedom, so NO single lineage supports a fitted decay rate. Every")
    print("  fit above is therefore pooled across constructs and runs, which")
    print("  is why the residual scatter is so large.")


def report(f):
    """Print the fit and everything needed to read it honestly."""
    if f is None:
        print("  not enough uncensored points to fit a line.")
        return
    print(f"  log10(Amp:Kan) = {f['slope']:+.3f} * round {f['intercept']:+.3f}")
    print(f"  fold change per round : {f['fold_per_round']:.3f}x "
          f"(95% CI {f['fold_lo']:.3f} - {f['fold_hi']:.3f})")
    if f["fold_per_round"] < 1:
        print(f"  i.e. the ratio falls about "
              f"{1.0 / f['fold_per_round']:.1f}x per round")
    if f.get("censored_fit"):
        print(f"  n={f['n']} points ({f['n_censored']} of them bounds, kept in "
              f"as inequalities), rounds {f['x_levels']}")
        print(f"  residual scatter sigma = {f['sigma']:.2f} in log10 units "
              f"(~{10 ** f['sigma']:.1f}x)")
        print("  (CI is a likelihood-ratio profile; R^2 is undefined when "
              "some terms are bounds)")
    else:
        print(f"  R^2 = {f['r2']:.3f} on n={f['n']} uncensored points, "
              f"df={f['df']}, rounds {f['x_levels']}")
        if f["weighted"]:
            print("  (weighted by the sheet's 95% intervals)")

    cc = censored_consistency(f)
    if cc:
        if f.get("censored_fit"):
            print("\n  Bounds (these were IN the fit; shown for inspection):")
        else:
            print("\n  Censored bounds, held out of the fit, checked "
                  "against it:")
        for con, run, idx, kind, bound, yhat, ok in cc:
            sym = ">" if kind == "lower" else "<"
            if f.get("censored_fit"):
                # Not a pass/fail here: the bound was IN the fit and is
                # modelled as a probability, so calling it "contradicted"
                # would overstate what a central line above it means.
                flag = "on expected side" if ok else "line sits past bound"
            else:
                flag = "consistent" if ok else "** CONTRADICTED **"
            print(f"    round {idx}  {con:<14} run {run}  "
                  f"measured {sym} {bound:<7.3g} fitted {yhat:<9.3g} {flag}")
        bad = [c for c in cc if not c[6]]
        if f.get("censored_fit"):
            # The Tobit fit does not have to pass BELOW every "<" bound: it
            # models each as a probability, so a central line sitting a little
            # above one is expected, not a failure. Only a gross violation
            # matters, and sigma already reports the scatter.
            print(f"\n  {len(cc) - len(bad)} of {len(cc)} bounds fall on the "
                  "expected side of the fitted line.")
            print("  (A censored fit models bounds probabilistically, so the "
                  "central line\n  need not clear every one -- unlike the "
                  "uncensored fit above, where each\n  contradiction is a "
                  "real out-of-sample failure.)")
        elif bad:
            print(f"\n  {len(bad)} of {len(cc)} bounds contradict the fit. A "
                  "line that predicts a value a plate\n  has ruled out is "
                  "wrong there, whatever its R^2.")
        else:
            print("\n  All bounds consistent.")


if __name__ == "__main__":
    print(__doc__)
    print("=" * 68)
    print("WHAT ONE LINEAGE CAN SUPPORT")
    print("=" * 68)
    report_coverage()
    print()
    print("=" * 68)
    print("FIT, all uncensored points (weighted) -- THE NAIVE FIT, BIASED")
    print("=" * 68)
    report(fit())
    print()
    print("=" * 68)
    print("FIT, unweighted -- does the weighting drive the answer?")
    print("=" * 68)
    report(fit(weighted=False))
    print()
    print("=" * 68)
    print("FIT, Phase 2 rounds only (Jumpstart excluded)")
    print("=" * 68)
    report(fit(skip_jumpstart=True))
    print()
    print("=" * 68)
    print("CENSORED (TOBIT) FIT -- every plotted point, bounds kept in")
    print("=" * 68)
    report(fit_censored())
