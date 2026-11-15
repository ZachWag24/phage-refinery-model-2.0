"""
Mechanistic Jumpstart -- the proposed extension from MODEL_REPORT.md section 9,
built fresh (no prior-session prototype was available to port).

WHAT THIS FILE IS FOR
------------------------
jumpstart.py is honest that its Amp:Kan output is phi/(1-phi), and phi is a
MEASURED input -- the real packaging competition between pHK002.2's intact f1
origin and M13KO7's KanR-interrupted origin is not modelled at all (see
jumpstart.py, section 3, and MODEL_REPORT.md section 7.1). This file is an
attempt to model that competition, so the Jumpstart ratio becomes an OUTPUT
instead of an input.

THE MECHANISM
----------------
Inside a dual-infected (I) cell, both genomes replicate toward a shared,
per-cell copy-number ceiling (they compete for the same replication
machinery). The phagemid enters the cell already established as a resident
plasmid, at Nres copies; the helper enters as a single injected genome. If
the two replicated at the SAME rate, this would just be a fixed multiplier
on whatever ratio walked in the door (see IDENTIFIABILITY below). The reason
to expect otherwise: M13KO7's origin is interrupted by the KanR insertion,
so its origin-driven replication may simply be slower than the phagemid's
intact origin. That gives the phagemid a growing head start over the course
of the incubation, not a fixed one -- and THAT is falsifiable in a way a
fixed phi never can be (see PREDICTIONS below).

Per dual-infected cell, tracking replicative template copies gp (phagemid)
and gh (helper):

    dgp/dt = rp * gp * (1 - (gp+gh)/Ncap)
    dgh/dt = rh * gh * (1 - (gp+gh)/Ncap)
    gp(infection) = Nres,   gh(infection) = 1

This is standard two-species logistic competition for one shared capacity
Ncap. Export (packaging) draws on these templates without consuming them
(same assumption as jumpstart.py's shared k2, and as the fixed-phi model):

    phi(t) = alpha*gp(t) / (alpha*gp(t) + gh(t))
    dPp/dt = phi(t)   * k2 * I
    dPh/dt = (1-phi(t)) * k2 * I

alpha is a residual packaging-affinity multiplier for anything that is NOT
copy number. See "WHY ALPHA IS CALIBRATED NOW, NOT Nres" below for what it
represents and why it is the calibrated parameter, not a fixed 1.

Nres IS NOW A REAL, MEASURED VALUE -- READ THIS BEFORE CHANGING IT
------------------------------------------------------------------------
The Benchling map for pHK002.2 shows an "ori" and a "rop" feature -- that is
the pBR322/ColE1-family signature, and rop is specifically the gene that
REPRESSES that origin's copy number. rop+ pBR322-family plasmids run about
15-20 copies/cell (Lee et al. 2005, cited in the Genspace iGEM 2016
measurement writeup; SnapGene's own pBR322 reference map independently
agrees, "~20 copies per cell"; a third source gives the same "15-20 copies
per cell" range). Nres defaults to 18.0 here -- the middle of that range,
not a free parameter anymore. Change it only if better data on pHK002.2
itself (not just "a pBR322-type plasmid" in general) comes in.

WHY ALPHA IS CALIBRATED NOW, NOT Nres
--------------------------------------------
The file's first version calibrated Nres (holding alpha=1, "copy number
alone") and got Nres~3.6 -- and flagged, honestly, that this was far below
typical pBR322-type copy numbers and would need revisiting once a real
number was available. It now is: Nres=18. Plugged in directly, with
replication_ratio=1.05 (phagemid still favoured) and alpha=1, the Jumpstart
would produce an Amp:Kan WELL ABOVE the measured 4.88 -- copy number alone
overshoots the data by roughly 4x, and the replication-rate mechanism only
ever pushes the ratio UP from wherever Nres starts it, never down. So alpha=1
is no longer viable at a real Nres. Rather than invent a specific biological
story for why phagemid packages less efficiently per template than its copy
number would predict (nothing here identifies one -- see the module's
conversation history for what was considered and ruled out for lack of a
mechanism), alpha now absorbs that gap explicitly, as an admitted unknown:
it is CALIBRATED to reproduce jump_ratio, holding the real Nres fixed. This
is the same move the file always made, just at the parameter that no longer
has a measurement to fall back on. Nres is real; alpha is now where "what we
don't know" lives, exactly as phi was where it lived before this file
existed.

At current defaults this calibration gives alpha ~ 0.21 -- roughly a 5x
packaging discount for the phagemid relative to what its raw copy-number
advantage alone would produce. That is a real, sizeable number, not a
rounding correction, and it is worth remembering it is doing exactly the
kind of unexplained work the file's first version was built to avoid
needing (see MODEL_REPORT.md section 7.8, where a similar-looking
helper-favoring packaging bias was rejected for lacking a mechanism). The
difference here is that this one is forced by a real Nres measurement, not
chosen to fit a trend -- but "forced" is not the same as "explained."

MEAN-FIELD APPROXIMATION, NAMED
-----------------------------------
Real cells get infected at different times and so sit at different points
along their own replication trajectory. This file does NOT track that
per-cell age structure -- it tracks population TOTALS Gp, Gh and uses the
POPULATION AVERAGE (Gp/I, Gh/I) as if every cell were at that average state,
including for the shared capacity term. That is an approximation, not a
measurement: a fully age-structured version would track a distribution over
infection times rather than one average. It is the same kind of
simplification jumpstart.py and phase2.py already make everywhere (no
eclipse delay, no cell-to-cell variation) -- named here because it is new to
THIS file, not inherited from them.

IDENTIFIABILITY (why this file exists, and its limit)
----------------------------------------------------------
Derived and checked (test_model.py):
if rp == rh, this model's output is EXACTLY alpha * Nres, constant --
independent of Ncap, T, k1, moi, or anything else. That is the fixed-phi
model wearing two extra parameters that cancel. The whole point of this file
is the case rp != rh, where Gp/Gh no longer stays at Nres: while B cells are
still being infected, fresh cells keep entering the dual-infected pool at
ratio Nres, holding the population average near Nres; once B is exhausted (in
this protocol, once kanamycin has killed every uninfected cell), no more
Nres-ratio cells arrive, and the replication-rate gap is free to compound
for the rest of the incubation. That is where the time- and Ncap-dependence
below comes from -- it is NOT built in as a rule, it falls out of the
equations once rp != rh.

PREDICTIONS THIS MODEL MAKES THAT A FIXED PHI CANNOT
----------------------------------------------------------
At the defaults below (Nres=18 real, replication_ratio=1.05, Ncap=300,
alpha calibrated to 0.2121 so the 18 h output matches the measured 4.88):

    T =  6 h : ratio ~ 4.64
    T = 18 h : ratio ~ 4.88   (the calibration point)
    T = 36 h : ratio ~ 4.91

    Ncap =   50 : ratio ~ 4.48
    Ncap =  300 : ratio ~ 4.88   (the calibration point)
    Ncap = 2000 : ratio ~ 5.30

Both rise, both weakly -- more capacity or more time gives the faster
replicator (phagemid) more room to pull ahead, on top of whatever alpha is
discounting it by. These numbers will shift if Nres is later revised from a
pHK002.2-specific measurement rather than "pBR322-type in general"; run
`--sensitivity` for the current numbers rather than trusting this snapshot
indefinitely.

WHAT THIS FILE DELIBERATELY DOES NOT TOUCH
-----------------------------------------------
Everything else is identical to jumpstart.py: same Holling-II infection
saturation, same growth penalty, same one-time kanamycin kill, same absence
of an eclipse delay. Only the phi mechanism changes. This file does not feed
into phase2.py directly -- it stands alone, an alternative account of the
Jumpstart stage only, but recursion.py's --mechanistic flag can swap it in
for the whole chain (see recursion.py).

USAGE
-----
    python3 jumpstart_mechanistic.py                    # real Nres, calibrated alpha
    python3 jumpstart_mechanistic.py --sensitivity       # the two prediction tables
    python3 jumpstart_mechanistic.py --alpha 0.3          # PREDICT instead of calibrate
    python3 jumpstart_mechanistic.py --Nres 20            # override Nres if better data comes in
    python3 jumpstart_mechanistic.py --replication-ratio 1.0   # collapse check: ratio == alpha*Nres
    python3 jumpstart_mechanistic.py --plot out.png --csv trace.csv
"""
import argparse
import csv as csvlib

import numpy as np


# ==============================================================================
# THE EQUATIONS.
# ==============================================================================
def derivatives(y, k1, h, k2, k3, k4, K, rp, rh, Ncap, Nres, alpha):
    """
    y = [B, I, Pp, Ph, Gp, Gh]. Gp, Gh are POPULATION TOTALS of replicative
    template copies across all dual-infected cells (not per-cell counts --
    see module docstring, MEAN-FIELD APPROXIMATION).
    """
    B, I, Pp, Ph, Gp, Gh = np.maximum(y, 0.0)

    # Infection and growth: identical to jumpstart.py.
    sat = 1.0 / (1.0 + k1 * h * Ph)
    infection = k1 * Ph * B * sat
    room = 1.0 - (B + I) / K
    dB = -infection + k3 * B * room
    dI = infection + k4 * I * room

    # Per-cell averages, guarded against I == 0 (true at t=0, before any
    # infection: Gp = Gh = 0 there too, so the guard never hides real state).
    if I > 1e-12:
        gp_avg, gh_avg = Gp / I, Gh / I
    else:
        gp_avg = gh_avg = 0.0

    # Both genomes replicate toward the SAME shared capacity Ncap -- this is
    # the "one set of replication machinery, split between two genomes"
    # assumption, the same shape as jumpstart.py's shared export rate k2.
    cap_room = max(0.0, 1.0 - (gp_avg + gh_avg) / Ncap)

    # dGp/dt = rp*gp_avg*cap_room*I + Nres*infection simplifies to the line
    # below because gp_avg*I == Gp (verified in test_model.py -- this
    # identity is also what makes the rp==rh collapse to alpha*Nres exact
    # rather than approximate; see IDENTIFIABILITY in the module docstring
    # for the derivation).
    dGp = rp * Gp * cap_room + Nres * infection
    dGh = rh * Gh * cap_room + 1.0 * infection

    denom = alpha * gp_avg + gh_avg
    phi = (alpha * gp_avg / denom) if denom > 0 else 0.5   # I=0 here too, so phi is moot
    dPp = phi * k2 * I
    dPh = (1.0 - phi) * k2 * I

    return np.array([dB, dI, dPp, dPh, dGp, dGh])


def _simulate(B0, moi, T, steps, k1, h, k2, k3, K, kan_at, growth_penalty,
              rp, rh, Ncap, Nres, alpha):
    k4 = growth_penalty * k3
    y0 = np.array([B0, 0.0, 0.0, moi * B0, 0.0, 0.0])
    dt = T / steps
    Y = np.zeros((steps + 1, 6))
    Y[0] = y0
    kan_step = None if kan_at is None else int(round(kan_at / dt))
    args = (k1, h, k2, k3, k4, K, rp, rh, Ncap, Nres, alpha)
    for i in range(steps):
        y = Y[i]
        s1 = derivatives(y, *args)
        s2 = derivatives(y + 0.5 * dt * s1, *args)
        s3 = derivatives(y + 0.5 * dt * s2, *args)
        s4 = derivatives(y + dt * s3, *args)
        Y[i + 1] = np.maximum(y + (dt / 6.0) * (s1 + 2 * s2 + 2 * s3 + s4), 0.0)
        if kan_step is not None and (i + 1) == kan_step:
            Y[i + 1, 0] = 0.0
    return np.arange(steps + 1) * dt, Y


def resolve_alpha(jump_ratio, lo=1e-4, hi=1e4, iters=40, **run_kwargs):
    """
    Find the alpha that makes this model's output ratio equal jump_ratio,
    holding every other parameter fixed via run_kwargs (Nres included --
    Nres is now a real measured value, not something this function touches).
    Bisection: the output ratio is monotonically increasing in alpha (more
    packaging affinity for the phagemid never lowers its share), so this
    always has exactly one solution in a wide enough bracket.

    This plays the role resolve_nres used to play in the file's first
    version -- see module docstring, "WHY ALPHA IS CALIBRATED NOW, NOT
    Nres" for why the calibrated parameter changed.
    """
    def ratio_for(a):
        _, Y = _simulate(alpha=a, **run_kwargs)
        Pp, Ph = Y[-1, 2], Y[-1, 3]
        return Pp / Ph if Ph > 0 else float("inf")

    r_lo, r_hi = ratio_for(lo), ratio_for(hi)
    if not np.isfinite(r_lo) and not np.isfinite(r_hi):
        raise ValueError(
            "alpha cannot be calibrated: the ratio is undefined (both bounds "
            "give inf) regardless of alpha -- this happens when no helper "
            "ever enters the system at all (e.g. moi=0), so Ph stays exactly "
            "0 and nothing depends on alpha. Pass --alpha explicitly if you "
            "want to run this regime; calibration has nothing to solve for here.")
    if not (r_lo <= jump_ratio <= r_hi):
        raise ValueError(
            f"jump_ratio={jump_ratio} is not reachable by varying alpha in "
            f"[{lo}, {hi}] (gives [{r_lo:.4g}, {r_hi:.4g}]); widen the "
            f"bracket, or reconsider Nres/replication_ratio/Ncap.")
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if ratio_for(mid) < jump_ratio:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run(B0=4.0e7, moi=2.5, T=1080.0, steps=4000, k1=3.0e-11, h=1.0, k2=3.0,
        k3=np.log(2) / 25.0, K=1.0e9, kan_at=75.0, growth_penalty=0.6,
        rp=None, replication_ratio=1.05, Ncap=300.0, Nres=18.0,
        alpha=None, jump_ratio=4.88):
    """
    Simulate the Jumpstart with the mechanistic (replication-competition)
    phi, from t=0 to t=T.

    rp                : phagemid's replication rate, /min. None (default)
                         uses k3 -- an arbitrary shared clock; only its
                         RATIO to rh (set by replication_ratio) matters for
                         the output ratio, not its absolute value.
    replication_ratio : rp/rh -- how much faster the phagemid's intact
                         origin replicates than M13KO7's interrupted one.
                         UNMEASURED. Default 1.05.
    Ncap              : shared per-cell copy-number ceiling both genomes
                         replicate toward. UNMEASURED. Default 300.
    Nres              : phagemid's resident copy number at infection.
                         Default 18.0 -- a REAL measured-backbone value now,
                         not a free parameter. See module docstring, "Nres
                         IS NOW A REAL, MEASURED VALUE".
    alpha             : residual packaging-affinity multiplier beyond copy
                         number. None (default): CALIBRATED so the output
                         ratio equals jump_ratio, holding Nres fixed at its
                         real value (see resolve_alpha). A number: used
                         directly, and the output ratio becomes a genuine
                         PREDICTION that need not equal jump_ratio.
    jump_ratio        : the measured Jumpstart Amp:Kan (particle-level, same
                         convention as jumpstart.py). Only used to calibrate
                         alpha when alpha=None.

    Returns (t, Y, alpha_used): Y has columns [B, I, Pp, Ph, Gp, Gh].
    """
    if rp is None:
        rp = k3
    rh = rp / replication_ratio
    common = dict(B0=B0, moi=moi, T=T, steps=steps, k1=k1, h=h, k2=k2, k3=k3,
                  K=K, kan_at=kan_at, growth_penalty=growth_penalty,
                  rp=rp, rh=rh, Ncap=Ncap, Nres=Nres)
    alpha_used = alpha if alpha is not None else resolve_alpha(jump_ratio, **common)
    t, Y = _simulate(alpha=alpha_used, **common)
    return t, Y, alpha_used


# ==============================================================================
# SENSITIVITY TABLES -- the two falsifiable predictions from the module
# docstring, reproduced on demand rather than hard-coded, so they track
# whatever parameters you actually pass.
# ==============================================================================
def sensitivity_tables(**kw):
    """
    Returns (t_table, ncap_table, alpha_used): (value, ratio) pairs, plus
    the alpha this was computed at.

    alpha is resolved ONCE, at the baseline T/Ncap (calibrated to jump_ratio
    if not given explicitly, holding Nres fixed at its real value), and then
    held FIXED while T or Ncap is swept. Re-calibrating alpha at every sweep
    point would force every single point back to jump_ratio by construction
    and make the table meaningless -- the whole point is to see what a FIXED
    mechanism predicts elsewhere.
    """
    T0 = kw.pop("T", 1080.0)
    Ncap0 = kw.pop("Ncap", 300.0)
    jump_ratio = kw.pop("jump_ratio", 4.88)
    alpha = kw.pop("alpha", None)
    if alpha is None:
        _, _, alpha = run(T=T0, Ncap=Ncap0, jump_ratio=jump_ratio, **kw)

    t_table = []
    for hours in (6, 18, 36):
        t, Y, _ = run(T=hours * 60.0, Ncap=Ncap0, alpha=alpha, **kw)
        t_table.append((hours, Y[-1, 2] / Y[-1, 3]))

    ncap_table = []
    for ncap in (50, 300, 2000):
        t, Y, _ = run(T=T0, Ncap=float(ncap), alpha=alpha, **kw)
        ncap_table.append((ncap, Y[-1, 2] / Y[-1, 3]))

    return t_table, ncap_table, alpha


def sensitivity_curves(n=15, **kw):
    """
    Denser version of sensitivity_tables, for plotting: n log-spaced points
    across 6-36 h and across Ncap=50-2000, at the SAME fixed (once-calibrated)
    alpha as sensitivity_tables uses. The 3-point table is for reading; this
    is for seeing that the underlying relationship is a smooth curve, not
    two straight-line segments joined at 18 h / Ncap=300.

    Returns (t_curve, ncap_curve, alpha_used).
    """
    T0 = kw.pop("T", 1080.0)
    Ncap0 = kw.pop("Ncap", 300.0)
    jump_ratio = kw.pop("jump_ratio", 4.88)
    alpha = kw.pop("alpha", None)
    if alpha is None:
        _, _, alpha = run(T=T0, Ncap=Ncap0, jump_ratio=jump_ratio, **kw)

    hours = np.linspace(6, 36, n)
    t_curve = []
    for hh in hours:
        _, Y, _ = run(T=hh * 60.0, Ncap=Ncap0, alpha=alpha, **kw)
        t_curve.append((hh, Y[-1, 2] / Y[-1, 3]))

    ncaps = np.geomspace(50, 2000, n)
    ncap_curve = []
    for nc in ncaps:
        _, Y, _ = run(T=T0, Ncap=nc, alpha=alpha, **kw)
        ncap_curve.append((nc, Y[-1, 2] / Y[-1, 3]))

    return t_curve, ncap_curve, alpha


def build_parser():
    p = argparse.ArgumentParser(
        description="Mechanistic Jumpstart: phi from replication competition, not a fixed input.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--jump-ratio", type=float, default=4.88, dest="jump_ratio",
        help="measured Jumpstart Amp:Kan, used to calibrate alpha when --alpha is not given")
    p.add_argument("--Nres", type=float, default=18.0,
        help="phagemid resident copy number at infection. REAL value now: pHK002.2 is "
             "pBR322-type with an intact rop gene, which runs ~15-20 copies/cell -- "
             "override only if better data on pHK002.2 specifically comes in")
    p.add_argument("--alpha", type=float, default=None,
        help="residual packaging-affinity multiplier beyond copy number. Unset = "
             "calibrated to --jump-ratio, holding Nres fixed (see module docstring, "
             "\"WHY ALPHA IS CALIBRATED NOW, NOT Nres\"); given = this becomes a "
             "PREDICTION, not a fit")
    p.add_argument("--replication-ratio", type=float, default=1.05, dest="replication_ratio",
        help="rp/rh: how much faster the phagemid's intact origin replicates "
             "than M13KO7's interrupted one. UNMEASURED. 1.0 = the collapse "
             "case (output becomes exactly alpha*Nres, constant)")
    p.add_argument("--Ncap", type=float, default=300.0,
        help="shared per-cell copy-number ceiling. UNMEASURED")
    p.add_argument("--rp", type=float, default=None,
        help="phagemid replication rate, /min (default: same value as k3 -- "
             "see run()'s docstring for why the absolute value barely matters)")
    p.add_argument("-T", "--time", type=float, default=1080.0,
        help="total simulated duration, minutes")
    p.add_argument("--sensitivity", action="store_true",
        help="print the two falsifiable prediction tables (ratio vs T, ratio vs Ncap) and exit")
    p.add_argument("--k1", type=float, default=3.0e-11)
    p.add_argument("--h", type=float, default=1.0)
    p.add_argument("-k", "--k2", type=float, default=3.0)
    p.add_argument("--k3", type=float, default=np.log(2) / 25.0)
    p.add_argument("--K", type=float, default=1.0e9)
    p.add_argument("--B0", type=float, default=4.0e7)
    p.add_argument("-M", "--moi", type=float, default=2.5)
    p.add_argument("--kan-at", type=float, default=75.0, dest="kan_at")
    p.add_argument("--no-kan", action="store_true")
    p.add_argument("--growth-penalty", type=float, default=0.6, dest="growth_penalty")
    p.add_argument("--steps", type=int, default=4000)
    p.add_argument("--csv", metavar="FILE", default=None)
    p.add_argument("--plot", metavar="FILE", default=None)
    p.add_argument("-q", "--quiet", action="store_true")
    return p


def main():
    a = build_parser().parse_args()
    kan_at = None if a.no_kan else a.kan_at
    common = dict(B0=a.B0, moi=a.moi, steps=a.steps, k1=a.k1, h=a.h, k2=a.k2,
                  k3=a.k3, K=a.K, kan_at=kan_at, growth_penalty=a.growth_penalty,
                  rp=a.rp, replication_ratio=a.replication_ratio, Ncap=a.Ncap,
                  Nres=a.Nres, alpha=a.alpha, jump_ratio=a.jump_ratio)

    if a.sensitivity:
        t_table, ncap_table, alpha_used = sensitivity_tables(T=a.time, **common)
        print("=" * 60 + "\nSENSITIVITY (the two falsifiable predictions)\n" + "=" * 60)
        print(f"  Nres = {a.Nres:g} (real), alpha = {alpha_used:.4f} "
              f"({'given' if a.alpha is not None else 'calibrated once at baseline, held fixed'})")
        print(f"  ratio vs culture length (Ncap={a.Ncap:g}):")
        for hours, r in t_table:
            print(f"    {hours:3d} h : {r:.4f}")
        print(f"  ratio vs per-cell copy cap (T={a.time:g} min):")
        for ncap, r in ncap_table:
            print(f"    Ncap={ncap:6.0f} : {r:.4f}")

        if a.plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            try:
                import plotstyle as ps
                ps.apply_style()
                c = ps.DUAL
            except ImportError:
                ps = None
                c = "#0f172a"
            t_curve, ncap_curve, _ = sensitivity_curves(T=a.time, **common)
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
            hrs, r_t = zip(*t_curve)
            ax1.plot(hrs, r_t, "-", color=c, lw=2.2)
            th, tr = zip(*t_table)
            ax1.plot(th, tr, "o", color=c, markersize=8, zorder=3)
            ax1.set_xlabel("culture length (h)"); ax1.set_ylabel("Amp:Kan")
            ax1.set_title(f"Ratio vs. culture length (Ncap={a.Ncap:g})", loc="left", weight="bold")
            ncaps, r_n = zip(*ncap_curve)
            ax2.plot(ncaps, r_n, "-", color=c, lw=2.2)
            nh, nr = zip(*ncap_table)
            ax2.plot(nh, nr, "o", color=c, markersize=8, zorder=3)
            ax2.set_xscale("log"); ax2.set_xlabel("per-cell copy cap Ncap"); ax2.set_ylabel("Amp:Kan")
            ax2.set_title(f"Ratio vs. copy cap (T={a.time:g} min)", loc="left", weight="bold")
            for ax in (ax1, ax2):
                if ps:
                    ps.finish_axes(ax, legend=False)
                else:
                    ax.grid(alpha=0.2)
                    for s_ in ("top", "right"):
                        ax.spines[s_].set_visible(False)
            fig.suptitle("Two falsifiable predictions a fixed phi cannot make "
                         f"(Nres={a.Nres:g} real, replication_ratio={a.replication_ratio:g}, "
                         f"alpha={alpha_used:.3f})",
                         fontsize=10, color="#555")
            fig.tight_layout()
            fig.savefig(a.plot, dpi=200, facecolor="white")
            print(f"\n  wrote {a.plot}")
        return

    t, Y, alpha_used = run(T=a.time, **common)
    B, I, Pp, Ph, Gp, Gh = Y[-1]
    ratio = Pp / Ph if Ph > 0 else float("inf")

    if not a.quiet:
        print("=" * 60 + "\nMECHANISTIC JUMPSTART\n" + "=" * 60)
        mode = "GIVEN (a prediction)" if a.alpha is not None else "CALIBRATED to jump_ratio"
        print(f"  Nres  = {a.Nres:g}  (REAL: pBR322+rop backbone, ~15-20 copies/cell)")
        print(f"  alpha = {alpha_used:.4f}  ({mode})")
        print(f"  replication_ratio (rp/rh) = {a.replication_ratio:g}   "
              f"({'COLLAPSE CASE: ratio == alpha*Nres exactly' if a.replication_ratio == 1.0 else 'unmeasured'})")
        print(f"  Ncap = {a.Ncap:g}")
        print("-" * 60)
        print(f"  final Pp = {Pp:.3e}   Ph = {Ph:.3e}   Gp = {Gp:.3e}   Gh = {Gh:.3e}")
        print(f"  Amp:Kan  = {ratio:.4f}" + ("" if a.alpha is None else
              f"   (measured jump_ratio was {a.jump_ratio:g} -- "
              f"{'matches' if abs(ratio - a.jump_ratio) < 0.05 * a.jump_ratio else 'DIFFERS from'} it)"))

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csvlib.writer(f)
            w.writerow(["t_min", "B", "I", "Pp", "Ph", "Gp", "Gh"])
            for row_t, row_y in zip(t, Y):
                w.writerow([row_t] + list(row_y))
        print(f"\n  wrote {a.csv}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            import plotstyle as ps
            ps.apply_style()
            c_Pp, c_Ph, c_gp, c_gh = ps.PHAGEMID, ps.HELPER, ps.PHAGEMID, ps.HELPER
        except ImportError:
            ps = None
            c_Pp, c_Ph, c_gp, c_gh = "#1f6fb4", "#d1495b", "#1f6fb4", "#d1495b"

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
        ax1.plot(t, Y[:, 2], "-", color=c_Pp, lw=2.2, label="Pp (Amp)")
        ax1.plot(t, Y[:, 3], "-", color=c_Ph, lw=2.2, label="Ph (Kan)")
        ax1.set_yscale("log"); ax1.set_xlabel("time (min)"); ax1.set_ylabel("particles / mL")
        ax1.set_title("Free particles", loc="left", weight="bold"); ax1.legend(frameon=False)

        I = Y[:, 1]
        gp_avg = np.divide(Y[:, 4], I, out=np.zeros_like(I), where=I > 1e-9)
        gh_avg = np.divide(Y[:, 5], I, out=np.zeros_like(I), where=I > 1e-9)
        ax2.plot(t, gp_avg, "--", color=c_gp, lw=1.8, label="gp (avg. phagemid templates/cell)")
        ax2.plot(t, gh_avg, "--", color=c_gh, lw=1.8, label="gh (avg. helper templates/cell)")
        ax2.set_xlabel("time (min)"); ax2.set_ylabel("templates / dual-infected cell")
        ax2.set_title("Per-cell template competition", loc="left", weight="bold")
        ax2.legend(frameon=False)

        for ax in (ax1, ax2):
            if ps:
                ps.finish_axes(ax, legend=False)
            else:
                ax.grid(alpha=0.2)
                for s_ in ("top", "right"):
                    ax.spines[s_].set_visible(False)
        fig.suptitle(f"Nres={a.Nres:g} (real), alpha={alpha_used:.4f}, "
                     f"replication_ratio={a.replication_ratio:g}, Ncap={a.Ncap:g}  ->  "
                     f"Amp:Kan={ratio:.3f}",
                     fontsize=10, color="#555")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200, facecolor="white")
        print(f"\n  wrote {a.plot}")


if __name__ == "__main__":
    main()
