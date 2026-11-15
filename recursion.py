"""
The full recursion: Jumpstart once, then repeated Phase 2 rounds, joined by
passage() between every pair of stages.

THIS FILE HAS NO MODEL LOGIC OF ITS OWN
------------------------------------------
Every equation lives in jumpstart.py or phase2.py; passage() (its own file)
does the 1:11 dilution, with genuine extinction via Poisson sampling. This file
calls them in order and hands each stage's output to the next:

    Jumpstart.run() --> passage() --> Phase2.run() --> passage() --> Phase2.run() --> ...
        (once)                          (round 1)                       (round 2)

The Jumpstart runs once (it is the kick-start: colony, helper, kanamycin).
Every later round starts from fresh HP17KO7 cells and whatever phage the last
round produced -- that is Phase 2, and it is what repeats.

TWO LINEAGES, FOUR POOLS
---------------------------
pHK002.2 particles (Pp) and M13KO7 particles (Ph) are the two lineages. In
Phase 2 a pHK002.2 particle makes a cell that produces fully modified phage and
an M13KO7 particle makes a cell that produces hybrid phage; the hybrid coat
infects better (see phase2.py). Each lineage is tracked in two coat states:
the WILD-TYPE particles that came out of the Jumpstart (no HP17KO7 there), which
only ever carry over, and the particles Phase 2 makes. All four pools are
Poisson-sampled independently by passage().

TREND PRESET (read this before quoting a number)
---------------------------------------------------
The defaults are chosen to reproduce the general TREND of the plates --
Amp:Kan falling round after round while the helper persists -- NOT to fit them:

    jump_ratio=4.88 the Jumpstart's PLATE Amp:Kan (measured: about 5, range roughly
                   3-16). phi is derived from it; you rarely need --phi.
    helper_edge=1.0 an M13KO7 particle infecting better than a pHK002.2 particle
                   wearing the SAME coat -- NO LONGER AN ACTIVE DEFAULT (was 1.5).
                   Never had a real mechanism behind it: infection is coat-mediated,
                   so two particles with an identical coat should attach and enter
                   equally well, whatever genome is inside. Kept as a knob (see
                   phase2.py, HELPER EDGE, and MODEL_REPORT.md section 7.10) but no
                   longer stacked on top of the coat gap by default.
    eps=0.05       fully modified coat infects at 5% of wild-type; hybrid a_h
                   comes from coat gene dosage (see phase2.py), ~73%, so the
                   helper out-infects the phagemid ~15x. UNMEASURED.
    m13ko7_copies=20  M13KO7's own p15A-maintained copy number, NOT wild-type
                   M13's ~130-300 replicative form -- its f1 origin is
                   attenuated by the KanR/p15A insertion. This was 200 until
                   2026-11-14. The change is not cosmetic: a_h 0.963 -> 0.729,
                   coat advantage 19.3x -> 14.6x, round 1 0.718 -> 0.914.
                   UNMEASURED for this construct; plausible range 12-200, and
                   at the low end the chain stops crossing 1:1 at round 1.
    origin_penalty=1/1.05 (~0.952)  M13KO7's KanR-interrupted origin also
                   costs it export speed in Phase 2 (see phase2.py, ORIGIN
                   PENALTY). The one remaining coat-independent,
                   M13KO7-specific effect -- unlike the old helper_edge, this
                   one has a stated mechanism. Its SIZE is deliberately TIED
                   to the mechanistic Jumpstart's replication_ratio (1.05):
                   both describe the same origin defect, so the model makes
                   ONE guess about how big that defect is, not two
                   independent ones. It was an untied 0.9 until 2026-11-14.
    size_effect=0  no production advantage for the phagemid (it may well have
                   one; set --size-effect 1 to include it, and eps must then be
                   smaller for the ratio to keep falling)
    sigma=0.05     infected cells adsorb at 5% of a naive cell. UNMEASURED.

WHEN THE RATIO FLIPS: in Phase 2 round 1 the particles are still mostly
wild-type carry-over from the Jumpstart, so the coat gap has not yet had a
full round to act -- the made particles only start to dominate from round 2
onward (see phase2.py's plot, "What a plate would show", for round 1 itself).
Even so, the coat gap alone (a_h~0.73 vs eps~0.05, a ~15x difference) is still
large enough that round 1 flips below 1 at the measured Jumpstart ratio:
colony ratio 4.88 (Jumpstart) -> 0.91 (round 1) -> 0.096 (round 2) ->
0.012 (round 3). Note how much less margin that leaves than the old 200-copy
default did (0.72 at round 1): the flip now depends on the copy number being
nearer 20 than 12. helper_edge is no longer needed for this to happen -- it
happens with helper_edge=1 (removed) exactly as shown above.

The plate ratio ("COLONIES") is what to compare with the plates: it weights
each pool by its infectivity. The tube ratio ("particles") is what is
physically there and falls more slowly. With eps=1 and size_effect=0 the two
lineages are interchangeable and the ratio stays exactly where the Jumpstart
left it.
"""
import argparse
import csv as csvlib

import numpy as np

import jumpstart
import jumpstart_mechanistic
import phase2
from phase2 import plate_counts, hybrid_infectivity
from passage import passage


def run_recursion(rounds=6, seed=1,
                  Lp=3.8, Lh=8.7, k2=3.0,
                  k1=3.0e-11, h=1.0, k3=np.log(2) / 25.0, growth_penalty=0.6,
                  K=1.0e9, B0_jump=4.0e7, B0_phase2=3.0e8,
                  moi=2.5, kan_at=75.0, T1=1080.0, T2=900.0,
                  vol_transfer=0.5, vol_new=5.5, steps=4000,
                  eps=0.05, m13ko7_copies=20.0, hp17_copies=8.0, sigma=0.05,
                  phi=None, size_effect=0.0, helper_edge=1.0, jump_ratio=4.88,
                  origin_penalty=1/1.05,
                  mechanistic=False, mech_replication_ratio=1.05, mech_Ncap=300.0,
                  mech_alpha=None, mech_Nres=18.0):
    """
    Run the Jumpstart once, then `rounds` Phase 2 rounds, each preceded by a
    passage() dilution of the previous stage's output.

    jump_ratio  : the Jumpstart's PLATE Amp:Kan -- the quantity the plates
                  measured (about 5, range roughly 3-16). jumpstart.py's own
                  default is the same measurement (its PHI_MEASURED); the only
                  difference is that here it is corrected for the helper edge
                  first, since a plate counts infectious units. Plates count infectious
                  units, so the Jumpstart's PARTICLE ratio is jump_ratio *
                  helper_edge, and phi is derived from that. Default 4.88.
    phi         : override. Phagemid's share of a dual cell's export (0-1).
                  None (default) derives it from jump_ratio and helper_edge;
                  "length" derives it from genome length, Lh/(Lh+Lp) = 0.696.
                  Ignored if mechanistic=True.
    eps, m13ko7_copies, hp17_copies, sigma, size_effect, helper_edge : Phase 2,
                  see phase2.py -- m13ko7_copies/hp17_copies derive the hybrid
                  coat's infectivity from real gene dosage (hybrid_infectivity()),
                  replacing the old fitted f_hybrid exponent.
                  helper_edge also weights the Jumpstart's own plate reading:
                  its plate ratio is the particle ratio divided by helper_edge.
    origin_penalty : multiplier on an Ih cell's export rate (k2h), for
                  M13KO7's interrupted origin costing it export speed even
                  without direct competition (see phase2.py's run() docstring
                  and MODEL_REPORT.md section 7.13). UNMEASURED; default 0.9
                  here is the "current results" assumption -- phase2.py's own
                  standalone default is neutral (1.0), same convention as eps
                  and helper_edge.
    mechanistic : use jumpstart_mechanistic.py's replication-competition phi
                  instead of jumpstart.py's fixed input (see that module's
                  docstring). Calibrated to the SAME particle-level target
                  (jump_ratio * helper_edge) that the fixed-phi path uses, so
                  the Jumpstart's own output stays in the same ballpark --
                  what changes is that phi is no longer a bare number handed
                  in, it is the output of a mechanism, calibrated at the
                  parameter (alpha) that does NOT have a real measurement
                  behind it (mech_Nres=18 now IS real -- see that module's
                  docstring). mech_alpha given explicitly instead of
                  calibrated turns this into a genuine prediction that need
                  not match jump_ratio.

    Returns a list of dicts, one per stage (the Jumpstart plus each round):
    stage, Pp, Ph, ratio (particles), colony_ratio (plate), and 'extinct'
    ("M13KO7" or "pHK002.2") if a passage() draw ever leaves a species with
    exactly zero particles -- the run stops there, and that IS the outcome.
    """
    rng = np.random.default_rng(seed)
    results = []

    if mechanistic:
        _, Ymech, mech_alpha_used = jumpstart_mechanistic.run(
            B0=B0_jump, moi=moi, T=T1, steps=steps, k1=k1, h=h, k2=k2, k3=k3, K=K,
            kan_at=kan_at, growth_penalty=growth_penalty,
            replication_ratio=mech_replication_ratio, Ncap=mech_Ncap, alpha=mech_alpha,
            Nres=mech_Nres, jump_ratio=jump_ratio * helper_edge,
        )
        B, I, Pp, Ph = Ymech[-1, 0], Ymech[-1, 1], Ymech[-1, 2], Ymech[-1, 3]
    else:
        if phi is None:                       # hold the measured PLATE ratio fixed
            pr = jump_ratio * helper_edge
            phi_js = pr / (1.0 + pr)
        elif phi == "length":
            phi_js = -1.0          # jumpstart.resolve_phi treats any negative number as
                                    # "derive from genome length" -- None means "measured
                                    # default" instead, which silently discarded this option
        else:
            phi_js = phi

        t, Y = jumpstart.run(
            B0=B0_jump, moi=moi, T=T1, steps=steps,
            k1=k1, h=h, Lp=Lp, Lh=Lh, k2=k2, k3=k3, K=K,
            kan_at=kan_at, growth_penalty=growth_penalty, phi=phi_js,
        )
        B, I, Pp, Ph = Y[-1]
    # Jumpstart particles are all WILD-TYPE coat (no HP17KO7 there).
    Pw, Pt, Hw, Hh = Pp, 0.0, Ph, 0.0
    r0 = Pp / Ph if Ph > 0 else float("inf")
    results.append({"stage": "Jumpstart", "Pp": Pp, "Ph": Ph, "ratio": r0,
                    "colony_ratio": r0 / helper_edge,
                    **({"mech_alpha": mech_alpha_used} if mechanistic else {})})

    for r in range(1, rounds + 1):
        Pw_in, Pt_in, Hw_in, Hh_in = passage(Pw, Pt, Hw, Hh, vol_transfer=vol_transfer,
                                             vol_new=vol_new, rng=rng)
        Pp_in, Ph_in = Pw_in + Pt_in, Hw_in + Hh_in
        if Ph_in <= 0:
            results.append({"stage": f"round {r}", "Pp": Pp_in, "Ph": 0.0, "ratio": float("inf"),
                            "colony_ratio": float("inf"), "extinct": "M13KO7"})
            break
        if Pp_in <= 0:
            results.append({"stage": f"round {r}", "Pp": 0.0, "Ph": Ph_in, "ratio": 0.0,
                            "colony_ratio": 0.0, "extinct": "pHK002.2"})
            break
        t, Y = phase2.run(
            Pw_in, Hw_in, Pt_in=Pt_in, Hh_in=Hh_in, B0=B0_phase2, T=T2, steps=steps,
            k1=k1, h=h, Lp=Lp, Lh=Lh, k2=k2, k3=k3, K=K,
            growth_penalty=growth_penalty, eps=eps,
            m13ko7_copies=m13ko7_copies, hp17_copies=hp17_copies,
            sigma=sigma, size_effect=size_effect, helper_edge=helper_edge,
            origin_penalty=origin_penalty,
        )
        N, Ip, Ih, Pw, Pt, Hw, Hh = Y[-1]
        Pp, Ph = Pw + Pt, Hw + Hh
        a_h = hybrid_infectivity(eps, m13ko7_copies, hp17_copies)
        amp, kan = plate_counts(Pw, Pt, Hw, Hh, eps, a_h, helper_edge)
        results.append({"stage": f"round {r}", "Pp": Pp, "Ph": Ph, "Pw": Pw, "Pt": Pt, "Hw": Hw, "Hh": Hh,
                        "ratio": Pp / Ph if Ph > 0 else float("inf"),
                        "colony_ratio": amp / kan if kan > 0 else float("inf")})
    return results


def run_band(low=3.0, high=16.0, center=None, rounds=6, seed=1, **kwargs):
    """
    Run the SAME chain three times -- same seed, same every other parameter
    -- varying only jump_ratio across the measured Jumpstart range (default
    3-16, per MODEL_REPORT.md section 6/10). Reusing the seed across all
    three means any difference between them comes ONLY from jump_ratio, not
    from different passage() draws landing differently by chance.

    This exists because a single line (the measured jump_ratio's central
    value, 4.88) understates how much this one unmeasured-but-bounded input
    moves every downstream number -- see MODEL_REPORT.md section 5's
    sensitivity table, where jump_ratio alone gives a 4.7x spread in the
    round-2 plate ratio. Quoting only the center invites treating 0.306 as
    more precise than it is.

    center=None (default) uses jump_ratio's own default, 4.88 -- the
    measured central value, not the midpoint of [low, high].

    Returns a list of dicts, one per stage, each with:
        stage, colony_low, colony_center, colony_high,
        ratio_low, ratio_center, ratio_high
    A scenario that goes extinct before another stops contributes NaN to
    every later stage in ITS OWN columns only, rather than truncating the
    whole table -- the other two scenarios may still have more to say.
    """
    if center is None:
        center = kwargs.pop("jump_ratio", 4.88)
    runs = {
        "low": run_recursion(rounds=rounds, seed=seed, jump_ratio=low, **kwargs),
        "center": run_recursion(rounds=rounds, seed=seed, jump_ratio=center, **kwargs),
        "high": run_recursion(rounds=rounds, seed=seed, jump_ratio=high, **kwargs),
    }
    n_stages = max(len(r) for r in runs.values())
    band = []
    for i in range(n_stages):
        row = {"stage": None}
        for key, res in runs.items():
            if i < len(res):
                row["stage"] = row["stage"] or res[i]["stage"]
                row[f"colony_{key}"] = res[i]["colony_ratio"]
                row[f"ratio_{key}"] = res[i]["ratio"]
            else:                                    # this scenario already went extinct
                row[f"colony_{key}"] = float("nan")
                row[f"ratio_{key}"] = float("nan")
        band.append(row)
    return band


def _fmt_ratio(x):
    """Readable across six orders of magnitude (two decimals printed every small ratio as 0.00).
    NaN (no data here -- e.g. this scenario already stopped in run_band()) is rendered
    distinctly from a genuine infinite ratio (M13KO7 extinct), which np.isfinite alone
    cannot tell apart."""
    if np.isnan(x):
        return "--"
    if not np.isfinite(x):
        return "inf"
    if x >= 100:
        return f"{x:,.0f}"
    if x >= 1:
        return f"{x:.2f}"
    return f"{x:.3g}"


def build_parser():
    p = argparse.ArgumentParser(
        description="Jumpstart once, then repeated Phase 2 rounds (two lineages).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    c = p.add_argument_group("commonly-changed parameters")
    c.add_argument("-e", "--eps", type=float, default=0.05,
                   help="infectivity of the fully modified coat relative to wild-type "
                        "(1 = no difference between lineages). Default chosen to show the "
                        "decline, not measured. Hybrid a_h comes from coat gene dosage, see "
                        "--m13ko7-copies/--hp17-copies")
    c.add_argument("--jump-ratio", type=float, default=4.88, dest="jump_ratio",
                   help="the Jumpstart's PLATE Amp:Kan (measured about 5, range roughly 3-16); "
                        "the Jumpstart's phagemid share phi is derived from it and helper_edge")
    c.add_argument("--phi", type=float, default=None,
                   help="override the derived Jumpstart phagemid share of export (0-1). "
                        "Negative = derive from genome length, Lh/(Lh+Lp)")
    c.add_argument("--mechanistic", action="store_true",
                   help="use jumpstart_mechanistic.py's replication-competition phi instead "
                        "of a fixed input, calibrated at alpha (see jumpstart_mechanistic.py's "
                        "module docstring; still unvalidated against real data -- see "
                        "MODEL_REPORT.md open question 1)")
    c.add_argument("--mech-replication-ratio", type=float, default=1.05, dest="mech_replication_ratio",
                   help="mechanistic Jumpstart only: rp/rh, UNMEASURED (see jumpstart_mechanistic.py)")
    c.add_argument("--mech-ncap", type=float, default=300.0, dest="mech_Ncap",
                   help="mechanistic Jumpstart only: shared per-cell copy cap, UNMEASURED")
    c.add_argument("--mech-nres", type=float, default=18.0, dest="mech_Nres",
                   help="mechanistic Jumpstart only: phagemid resident copy number. REAL value "
                        "now (pHK002.2 is pBR322-type with an intact rop gene, ~15-20 copies/cell) "
                        "-- override only if better pHK002.2-specific data comes in")
    c.add_argument("--mech-alpha", type=float, default=None, dest="mech_alpha",
                   help="mechanistic Jumpstart only: residual packaging-affinity multiplier. "
                        "Unset = calibrated to jump_ratio, holding mech-nres fixed; given = "
                        "this becomes a genuine prediction instead of a calibration")
    c.add_argument("--size-effect", type=float, default=0.0, dest="size_effect",
                   help="exponent on (Lh/Lp) for the phagemid's export rate in Phase 2 "
                        "(0 = none, 1 = rate proportional to 1/length)")
    c.add_argument("--helper-edge", type=float, default=1.0, dest="helper_edge",
                   help="how many times better an M13KO7 particle infects than a pHK002.2 "
                        "particle with the SAME coat (1 = no difference, and no longer "
                        "an active default -- see origin_penalty for the mechanism with "
                        "actual backing)")
    c.add_argument("--sigma", type=float, default=0.05,
                   help="adsorption rate of an infected cell relative to a naive cell")
    c.add_argument("-r", "--rounds", type=int, default=6, help="Phase 2 rounds after the Jumpstart")
    c.add_argument("-s", "--seed", type=int, default=1, help="random seed for passage()")
    c.add_argument("--band", action="store_true",
                   help="report the plate ratio as a band across the measured Jumpstart "
                        "range (--jump-ratio-low to --jump-ratio-high) instead of a single "
                        "line at --jump-ratio. Same seed used for all three, so the spread "
                        "reflects jump_ratio alone, not passage() noise")
    c.add_argument("--jump-ratio-low", type=float, default=3.0, dest="jump_ratio_low",
                   help="low end of the measured Jumpstart range, used only with --band")
    c.add_argument("--jump-ratio-high", type=float, default=16.0, dest="jump_ratio_high",
                   help="high end of the measured Jumpstart range, used only with --band")
    f = p.add_argument_group("every other model parameter")
    f.add_argument("--m13ko7-copies", type=float, default=20.0, dest="m13ko7_copies",
                   help="M13 replicative-form DNA copy number in an infected cell "
                        "(typical published range ~130-300; not construct-specific)")
    f.add_argument("--hp17-copies", type=float, default=8.0, dest="hp17_copies",
                   help="HP17's own copy number (Rep101/pSC101 origin, typical published "
                        "range ~5-12; not construct-specific)")
    f.add_argument("--origin-penalty", type=float, default=1/1.05, dest="origin_penalty",
                   help="multiplier on an Ih cell's export rate, for M13KO7's interrupted "
                        "origin costing it export speed even without competition "
                        "(1 = no penalty). UNMEASURED")
    f.add_argument("-L", "--Lp", type=float, default=3.8, help="phagemid genome length, kb")
    f.add_argument("--Lh", type=float, default=8.7, help="M13KO7 genome length, kb")
    f.add_argument("-k", "--k2", type=float, default=3.0, help="export rate of an infected cell, particles/min")
    f.add_argument("--k1", type=float, default=3.0e-11)
    f.add_argument("--h", type=float, default=1.0)
    f.add_argument("--k3", type=float, default=np.log(2) / 25.0)
    f.add_argument("--growth-penalty", type=float, default=0.6, dest="growth_penalty")
    f.add_argument("--K", type=float, default=1.0e9)
    f.add_argument("--B0-jump", type=float, default=4.0e7, dest="B0_jump")
    f.add_argument("--B0-phase2", type=float, default=3.0e8, dest="B0_phase2")
    f.add_argument("--moi", type=float, default=2.5)
    f.add_argument("--kan-at", type=float, default=75.0, dest="kan_at")
    f.add_argument("--T1", type=float, default=1080.0)
    f.add_argument("--T2", type=float, default=900.0)
    f.add_argument("--vol-transfer", type=float, default=0.5, dest="vol_transfer")
    f.add_argument("--vol-new", type=float, default=5.5, dest="vol_new")
    f.add_argument("--steps", type=int, default=4000)
    o = p.add_argument_group("output")
    o.add_argument("--csv", metavar="FILE", default=None)
    o.add_argument("--plot", metavar="FILE", default=None)
    o.add_argument("--overlay-data", action="store_true", dest="overlay_data",
        help="mark the measured plate ratios on the ratio panel, aligned so "
             "the Jumpstart sits on lab R1. Reference context for where the "
             "preset came from, not a fit test. Every point is drawn the "
             "same, so the Jumpstart value jump_ratio was calibrated from "
             "is not marked -- say so in any caption. See labdata.py.")
    o.add_argument("--fit", action="store_true", dest="fit_line",
        help="fit a straight line (exponential decay) through the measured "
             "points and draw it on the ratio panel. Uses a CENSORED (Tobit) "
             "fit so the one-sided bounds are kept in as inequalities -- "
             "dropping them biases the slope badly, because a plate is "
             "censored precisely when the ratio collapses. Implies "
             "--overlay-data. See trendfit.py for what the fit can and "
             "cannot support.")
    o.add_argument("--compare", metavar="FILE", default=None,
        help="write the PRESENTATION figure: measured plate ratios on the "
             "left, the model's predicted ratio on the right, same axes. "
             "One message -- Amp:Kan falls round after round, and the model "
             "predicts it. Drops the titer panel and the particle curve, "
             "which do not serve that comparison. Use --plot for the fuller "
             "diagnostic figure.")
    o.add_argument("-q", "--quiet", action="store_true")
    return p


def _draw_measured(ax, mk, colour=None, spread=0.09, max_span=0.52):
    """
    Draw the measured points on `ax`. ONE marker shape for every point.

    EVERY point is drawn identically: same filled circle, same size, same
    colour as the rest of its lineage. Nothing on a marker distinguishes:

      * the ANCHOR -- the Jumpstart value jump_ratio was calibrated from. It
        now looks like a prediction the model matched, when in fact the model
        was set from it. Whatever uses this figure has to say so in words.
      * the six one-sided BOUNDS. A ">63" plotted as a plain point reads as
        "63 was measured" when it means "63 or anything above"; a "<0.016"
        reads as a measurement when the Amp count was simply zero.

    Both facts are still carried per point (`is_anchor`, `bound`) for any
    caller that wants them back. The only surviving visual hint is that a
    bound has no interval bar, so a bare circle with no bar is usually
    censored.

    Points sharing a round are spread sideways, with the total width capped so
    a crowded round cannot bleed into its neighbours.

    COLOUR identifies the LINEAGE -- one construct in one run -- so the two
    samples that span all three rounds can be followed across the plot, and
    the spread within a round is visibly several constructs rather than
    scatter on one measurement. Points of the same lineage are joined by a
    faint line for the same reason: without it a colour is just decoration.
    Pass `colour` to force a single colour for every point instead.

    Returns the list of y values actually drawn, for callers setting limits.
    """
    import labdata as _ld
    cmap = _ld.lineage_colours()

    def _col(m):
        return colour if colour else cmap.get(m["lineage"], "#ea580c")
    from collections import defaultdict
    by_round = defaultdict(list)
    for m in mk:
        by_round[m["index"]].append(m)

    drawn = []
    for m in mk:
        grp = by_round[m["index"]]
        n, k = len(grp), grp.index(m)
        span = min(spread * (n - 1), max_span) if n > 1 else 0.0
        x = m["index"] + ((k / (n - 1) - 0.5) * span if n > 1 else 0.0)
        y = m["value"]
        if not (y > 0):
            continue
        drawn.append(y)

        c = _col(m)
        ax.plot([x], [y], marker="o", markersize=7.0, linestyle="none",
                markerfacecolor=c, markeredgecolor=c, markeredgewidth=1.6,
                zorder=6)

        # 95% interval, where the sheet gives one.
        lo, hi = m["lo"], m["hi"]
        if lo is not None and hi is not None and lo > 0:
            ax.plot([x, x], [lo, hi], color=c, lw=1.3, alpha=0.55,
                    zorder=5)
            drawn += [lo, hi]

        # Limit arrow. Drawn in axes-fraction length so it looks the same on
        # any y range, rather than a fixed number of decades that would be a
        # stub on one figure and a streak on another.
    # Join the points of each lineage that appears in more than one round.
    # Only two lineages qualify, which is itself the point: everything else
    # is a single round and cannot show a trend on its own.
    tracks = defaultdict(list)
    for m in mk:
        if m["value"] > 0:
            grp = by_round[m["index"]]
            n, k = len(grp), grp.index(m)
            span = min(spread * (n - 1), max_span) if n > 1 else 0.0
            xx = m["index"] + ((k / (n - 1) - 0.5) * span if n > 1 else 0.0)
            tracks[m["lineage"]].append((xx, m["value"], m))
    for key, pts in tracks.items():
        if len(pts) < 2:
            continue
        pts.sort(key=lambda t: t[0])
        ax.plot([t[0] for t in pts], [t[1] for t in pts],
                color=(colour if colour else cmap.get(key, "#ea580c")),
                lw=1.2, alpha=0.45, zorder=4)
    return drawn


def _plot_compare(res, a, path):
    """
    The presentation figure: measured beside model, same axes, one message.

    MEASURED (left) and MODEL (right) share a y axis and an x axis so the eye
    compares them directly -- the point being that Amp:Kan falls round after
    round on the plates, and the model predicts the same fall. Everything not
    serving that comparison is left out: no titer panel, no particle curve, no
    parameter banner.

    Both panels span every simulated round, so the measured side visibly stops
    after round 2. That is deliberate: it shows how far the data actually goes
    rather than implying agreement over rounds nobody plated. A light band
    marks the measured rounds on the model panel for the same reason.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    try:
        import plotstyle as ps
        ps.apply_style()
        c_D = ps.DUAL
    except ImportError:
        ps = None
        c_D = "#0f172a"
    try:
        import labdata
    except ImportError:
        labdata = None

    MEAS, FITCOL, GUIDE = "#ea580c", "#475569", "#94a3b8"
    # This figure plots the COLONY ratio only -- the quantity plate counts are
    # comparable with -- so it keeps the dark line. Green is reserved for the
    # particle curve, which appears on the diagnostic figure.
    MODEL_LINE = c_D
    xs = list(range(len(res)))
    # Short tick labels: the 45-degree rotated "round 5" strings were the
    # noisiest thing on the old figure, and the axis title already says round.
    labels = ["JS" if r["stage"] == "Jumpstart"
              else r["stage"].replace("round ", "") for r in res]

    fig, (axm, axd) = plt.subplots(1, 2, figsize=(10.4, 5.0), sharey=True)

    mk = []
    if labdata is not None:
        mk = [m for m in labdata.markers() if m["index"] < len(res)]

    # ---- left: what the plates showed ------------------------------------
    _draw_measured(axm, mk)      # colour per lineage

    # A fitted trend through the measured points, if one can be had.
    fit_note = None
    if mk:
        try:
            import trendfit
            ft = trendfit.fit_censored()
        except Exception:
            ft = None
        if ft is not None:
            fx = np.array(ft["x_levels"], float)
            axm.plot(fx, 10.0 ** (ft["slope"] * fx + ft["intercept"]),
                     color=FITCOL, lw=2.0, dashes=(7, 3), zorder=4)
            fit_note = f"fitted: {1/ft['fold_per_round']:.0f}x down per round"

    # ---- right: what the model predicts ----------------------------------
    nz = lambda v: v if (np.isfinite(v) and v > 0) else np.nan
    axd.plot(xs, [nz(r["colony_ratio"]) for r in res], "o-",
             color=MODEL_LINE, lw=2.4, markersize=5.5, zorder=5)
    # How far to draw. Running 8 rounds but plating 3 leaves most of the
    # measured panel empty, which reads as missing rather than informative.
    # Crop both panels to just past the measured rounds and state the model's
    # continuation in words -- the number carries better than 5 blank rounds.
    last_meas = max((m["index"] for m in mk), default=None)
    if last_meas is not None and last_meas + 2 < len(res):
        xmax = last_meas + 1
        tail = res[-1]["colony_ratio"]
        if np.isfinite(tail) and tail > 0:
            axd.annotate(f"keeps falling to {tail:.0e}\nby round "
                         f"{len(res) - 1}",
                         xy=(xmax, 10.0 ** (np.log10(
                             max(res[xmax]["colony_ratio"], 1e-300)) - 0.3)),
                         fontsize=8.5, color=MODEL_LINE, alpha=0.75,
                         ha="right", va="top")
    else:
        xmax = len(res) - 1

    # ---- shared framing ---------------------------------------------------
    for ax, title in ((axm, "MEASURED on plates"), (axd, "MODEL prediction")):
        ax.set_yscale("log")
        ax.axhline(1.0, color=GUIDE, linestyle=":", lw=1.2, zorder=1)
        ax.set_xticks(xs[:xmax + 1])
        ax.set_xticklabels(labels[:xmax + 1])
        ax.set_xlabel("round")
        ax.set_xlim(-0.6, xmax + 0.4)
        ax.set_title(title, loc="left", weight="bold", fontsize=11)
        if ps:
            ps.finish_axes(ax, legend=False)
        else:
            ax.grid(alpha=0.2)
            for s_ in ("top", "right"):
                ax.spines[s_].set_visible(False)
    axm.set_ylabel("Amp : Kan")
    # Scale y to what is actually VISIBLE. The model's later rounds are still
    # in the plotted series even though the x axis is cropped, and matplotlib
    # autoscales to them -- which left most of both panels blank.
    vis = [r["colony_ratio"] for r in res[:xmax + 1]
           if np.isfinite(r["colony_ratio"]) and r["colony_ratio"] > 0]
    for m in mk:
        vis += [v for v in (m["value"], m["lo"], m["hi"])
                if v is not None and v > 0]
    if vis:
        axm.set_ylim(min(vis) / 6.0, max(vis) * 6.0)
    # One shared note for the 1:1 line rather than one per panel.
    axm.text(xmax + 0.35, 1.0, "1:1 ", color=GUIDE, fontsize=8.5,
             va="bottom", ha="right")

    # One shape throughout; COLOUR names the lineage. The key therefore lists
    # the constructs, which is the information the single-colour version threw
    # away -- then two short notes for the only non-colour variations.
    handles, lbls = [], []
    if mk:
        cmap = labdata.lineage_colours()
        seen = []
        for m in mk:
            if m["lineage"] in seen:
                continue
            seen.append(m["lineage"])
            handles.append(Line2D([], [], marker="o", linestyle="-", lw=1.2,
                                  markersize=7.0,
                                  color=cmap[m["lineage"]],
                                  markeredgewidth=1.6))
            lbls.append(m["lineage_label"])

    if fit_note:
        handles.append(Line2D([], [], color=FITCOL, lw=2.0, dashes=(7, 3)))
        lbls.append(fit_note)
    # The key goes BELOW the figure, not inside a panel. In-panel it covered
    # the round-1 cluster -- and on a figure whose whole job is showing those
    # points, hiding any of them to make room for their own legend is the one
    # thing it must not do.
    fig.legend(handles, lbls, loc="lower center", ncol=4, frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, -0.01),
               columnspacing=1.5, handletextpad=0.6)

    fig.suptitle("Amp : Kan falls with every round — "
                 "on the plates, and in the model",
                 fontsize=12.5, weight="bold", color="#1e293b")
    # Leave room at the bottom for the shared key and at the top for the title.
    fig.tight_layout(rect=(0, 0.11, 1, 0.94))
    fig.savefig(path, dpi=200, facecolor="white")
    print(f"\n  wrote {path}")


def _main_band(a):
    """CLI path for --band: reports the plate ratio as a range across the measured
    Jumpstart split instead of a single line. See run_band()'s docstring."""
    common = dict(Lp=a.Lp, Lh=a.Lh, k2=a.k2, k1=a.k1, h=a.h, k3=a.k3,
                  growth_penalty=a.growth_penalty, K=a.K, B0_jump=a.B0_jump, B0_phase2=a.B0_phase2,
                  moi=a.moi, kan_at=a.kan_at, T1=a.T1, T2=a.T2, vol_transfer=a.vol_transfer,
                  vol_new=a.vol_new, steps=a.steps, eps=a.eps,
                  m13ko7_copies=a.m13ko7_copies, hp17_copies=a.hp17_copies, sigma=a.sigma,
                  size_effect=a.size_effect, helper_edge=a.helper_edge, origin_penalty=a.origin_penalty,
                  mechanistic=a.mechanistic, mech_replication_ratio=a.mech_replication_ratio,
                  mech_Ncap=a.mech_Ncap, mech_alpha=a.mech_alpha, mech_Nres=a.mech_Nres)
    band = run_band(low=a.jump_ratio_low, high=a.jump_ratio_high, center=a.jump_ratio,
                    rounds=a.rounds, seed=a.seed, **common)
    a_h = hybrid_infectivity(a.eps, a.m13ko7_copies, a.hp17_copies)

    if not a.quiet:
        print("=" * 64 + "\nCONFIG (BAND MODE)\n" + "=" * 64)
        print(f"  Jumpstart plate ratio band = {a.jump_ratio_low:g} - {a.jump_ratio_high:g}  "
              f"(center {a.jump_ratio:g}, the measured value)")
        print(f"  fully modified coat eps = {a.eps:g}; hybrid a_h = {a_h:.3g}"
              f"  -> helper out-infects phagemid {a_h / a.eps:.3g}x")
        print(f"  helper edge = {a.helper_edge:g}x  (M13KO7 vs pHK002.2, same coat)")
        print(f"  sigma = {a.sigma:g}, Phase 2 size effect = {a.size_effect:g}, "
              f"origin penalty = {a.origin_penalty:g}x on Ih export (UNMEASURED)")
        print(f"  rounds, seed    = {a.rounds}, {a.seed}  (same seed for all three -- "
              "spread below is jump_ratio alone, not passage() noise)")
        print("\n" + "=" * 64 + "\nRESULT -- PLATE (COLONY) RATIO BAND\n" + "=" * 64)
        print(f"  {'stage':>10} {'low':>10} {'center':>10} {'high':>10}")
        for r in band:
            print(f"  {r['stage']:>10} {_fmt_ratio(r['colony_low']):>10} "
                  f"{_fmt_ratio(r['colony_center']):>10} {_fmt_ratio(r['colony_high']):>10}")
        print(f"\n  band = jump_ratio in [{a.jump_ratio_low:g}, {a.jump_ratio_high:g}], the measured "
              "Jumpstart range (MODEL_REPORT.md section 6); the center is not the midpoint of the "
              "band, it is the separately measured central value.")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csvlib.writer(f)
            w.writerow(["stage", "colony_low", "colony_center", "colony_high",
                        "particle_low", "particle_center", "particle_high"])
            for r in band:
                w.writerow([r["stage"], r["colony_low"], r["colony_center"], r["colony_high"],
                            r["ratio_low"], r["ratio_center"], r["ratio_high"]])
        print(f"\n  wrote {a.csv}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            import plotstyle as ps
            ps.apply_style()
            c_D = ps.DUAL
        except ImportError:
            ps = None
            c_D = "#0f172a"
        nz = lambda v: v if (np.isfinite(v) and v > 0) else np.nan
        xs = list(range(len(band)))
        lo = [nz(r["colony_low"]) for r in band]
        ct = [nz(r["colony_center"]) for r in band]
        hi = [nz(r["colony_high"]) for r in band]
        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        ax.fill_between(xs, lo, hi, color=c_D, alpha=0.15,
                        label=f"measured range ({a.jump_ratio_low:g}-{a.jump_ratio_high:g})")
        ax.plot(xs, ct, "o-", color=c_D, lw=2.4, label=f"center ({a.jump_ratio:g})")
        ax.axhline(1.0, color="#cbd5e1", lw=1)
        ax.set_yscale("log"); ax.set_ylabel("Amp:Kan (plate)")
        ax.set_title("Plate ratio, banded across the measured Jumpstart range", loc="left", weight="bold")
        ax.set_xticks(xs); ax.set_xticklabels([r["stage"] for r in band], rotation=45, ha="right")
        ax.legend(frameon=False)
        if ps:
            ps.finish_axes(ax, legend=False)
        else:
            ax.grid(alpha=0.2)
            for s_ in ("top", "right"):
                ax.spines[s_].set_visible(False)
        fig.suptitle(f"eps={a.eps:g} (hybrid a_h={a_h:.2g}), sigma={a.sigma:g}, "
                     f"size effect={a.size_effect:g}", fontsize=10, color="#555")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200, facecolor="white")
        print(f"\n  wrote {a.plot}")


def main():
    a = build_parser().parse_args()
    # Fitting a line to the data is meaningless without the data on the plot,
    # so --fit turns the overlay on rather than silently drawing a line
    # through points that are not shown.
    if a.fit_line:
        a.overlay_data = True
    if a.band:
        _main_band(a)
        return
    phi = None if a.phi is None else ("length" if a.phi < 0 else a.phi)
    res = run_recursion(
        rounds=a.rounds, seed=a.seed, Lp=a.Lp, Lh=a.Lh, k2=a.k2, k1=a.k1, h=a.h, k3=a.k3,
        growth_penalty=a.growth_penalty, K=a.K, B0_jump=a.B0_jump, B0_phase2=a.B0_phase2,
        moi=a.moi, kan_at=a.kan_at, T1=a.T1, T2=a.T2, vol_transfer=a.vol_transfer,
        vol_new=a.vol_new, steps=a.steps, eps=a.eps,
        m13ko7_copies=a.m13ko7_copies, hp17_copies=a.hp17_copies, sigma=a.sigma,
        phi=phi, size_effect=a.size_effect, helper_edge=a.helper_edge, jump_ratio=a.jump_ratio,
        origin_penalty=a.origin_penalty,
        mechanistic=a.mechanistic, mech_replication_ratio=a.mech_replication_ratio,
        mech_Ncap=a.mech_Ncap, mech_alpha=a.mech_alpha, mech_Nres=a.mech_Nres)
    a_h = hybrid_infectivity(a.eps, a.m13ko7_copies, a.hp17_copies)

    if a.mechanistic:
        mech_alpha_used = res[0].get("mech_alpha")
        mode = "GIVEN (a prediction)" if a.mech_alpha is not None else "calibrated to jump_ratio"
        jumpstart_desc = f"Jumpstart MECHANISTIC (alpha {mode})"
    else:
        if phi is None:
            _pr = a.jump_ratio * a.helper_edge
            phi_used = _pr / (1.0 + _pr)
        elif phi == "length":
            phi_used = a.Lh / (a.Lh + a.Lp)
        else:
            phi_used = phi
        jumpstart_desc = f"Jumpstart phi={phi_used:.2f}"

    if not a.quiet:
        print("=" * 64 + "\nCONFIG\n" + "=" * 64)
        if a.mechanistic:
            print(f"  Jumpstart phi   = MECHANISTIC (jumpstart_mechanistic.py; alpha {mode})")
            print(f"    Nres={a.mech_Nres:g} (REAL, pBR322+rop), "
                  f"replication_ratio={a.mech_replication_ratio:g} (unmeasured), "
                  f"Ncap={a.mech_Ncap:g} (unmeasured), alpha={mech_alpha_used:.4f}")
            print(f"    target jump_ratio (plate) = {a.jump_ratio:g}; Jumpstart's own colony_ratio "
                  f"below shows what the mechanism actually produced")
        else:
            _pp = phi_used / (1 - phi_used) if phi_used < 1 else float('inf')
            print(f"  Jumpstart phi   = {phi_used:.3f}  (particle ratio {_pp:.2f}:1; plate ratio "
                  f"{_pp / a.helper_edge:.2f}:1)")
        print(f"  fully modified coat eps = {a.eps:g}; hybrid a_h = {a_h:.3g}"
              f"  -> helper out-infects phagemid {a_h / a.eps:.3g}x")
        print(f"  helper edge = {a.helper_edge:g}x  (M13KO7 vs pHK002.2, same coat)")
        print(f"  sigma = {a.sigma:g}, Phase 2 size effect = {a.size_effect:g}, "
              f"origin penalty = {a.origin_penalty:g}x on Ih export (UNMEASURED)")
        print(f"  rounds, seed    = {a.rounds}, {a.seed}")
        print("\n" + "=" * 64 + "\nRESULT\n" + "=" * 64)
        print(f"  {'stage':>10} {'Pp':>10} {'Ph':>10} {'particles':>11} {'COLONIES':>11}")
        for r in res:
            note = f"   <- {r['extinct']} EXTINCT" if r.get("extinct") else ""
            print(f"  {r['stage']:>10} {r['Pp']:10.3e} {r['Ph']:10.3e} "
                  f"{_fmt_ratio(r['ratio']):>11} {_fmt_ratio(r['colony_ratio']):>11}{note}")
        print("\n  particles = Amp:Kan physically in the tube;  COLONIES = what a plate shows"
              " (weights each pool by its infectivity).")

    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csvlib.writer(f)
            w.writerow(["stage", "Pp", "Ph", "particle_ratio", "colony_ratio", "extinct"])
            for r in res:
                w.writerow([r["stage"], r["Pp"], r["Ph"], r["ratio"], r["colony_ratio"], r.get("extinct", "")])
        print(f"\n  wrote {a.csv}")

    if a.compare:
        _plot_compare(res, a, a.compare)

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            import plotstyle as ps          # shared palette; figures still work without it
            ps.apply_style()
            c_P, c_H, c_D, c_N = ps.PHAGEMID, ps.HELPER, ps.DUAL, ps.NAIVE
        except ImportError:
            ps = None
            c_P, c_H, c_D, c_N = "#1f6fb4", "#d1495b", "#0f172a", "#94a3b8"
        xs = list(range(len(res)))
        nz = lambda v: v if (np.isfinite(v) and v > 0) else np.nan
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.1))
        ax1.plot(xs, [nz(r["Pp"]) for r in res], "o-", color=c_P, lw=2.2, label="Pp (Amp)")
        ax1.plot(xs, [nz(r["Ph"]) for r in res], "s-", color=c_H, lw=2.2, label="Ph (Kan)")
        ax1.set_yscale("log"); ax1.set_ylabel("particles / mL"); ax1.set_title("Titer across rounds", loc="left", weight="bold")
        # Declared before the overlay so the legend code below can rely on
        # them whether or not the overlay ran at all.
        overlay_drawn = False
        MODEL_PARTICLES = "#059669"   # tube (particles), in clear green
        has_bounds = False
        fit_short = None
        MEAS = "#ea580c"
        FITCOL = "#475569"
        # Smaller markers than the titer panel, and a lighter tube curve: the
        # plate curve is the one the plate counts are comparable with, so it
        # leads, and the measured points need to stay the most prominent
        # things on the panel rather than competing with two marker trails.
        ax2.plot(xs, [nz(r["colony_ratio"]) for r in res], "o-", color=c_D,
                 lw=2.4, markersize=4.5, label="plate (COLONIES)")
        ax2.plot(xs, [nz(r["ratio"]) for r in res], "o--",
                 color=MODEL_PARTICLES, lw=2.0, markersize=4.0,
                 label="tube (particles)")
        ax2.axhline(1.0, color="#cbd5e1", lw=1)

        # ---- measured plate ratios, as reference markers ------------------
        # Context for where the preset came from, NOT a goodness-of-fit test.
        # They go on the COLONIES curve: a plate counts transducing units, so
        # the particle curve is a different quantity.
        # Alignment: the Jumpstart IS the first plated stage, so lab round N
        # sits at index N-1. labdata.py carries that map; do not re-derive it.
        # All constructs and both runs are drawn as ONE series -- see
        # labdata.py before reading the vertical spread as error.
        if a.overlay_data:
            try:
                import labdata
            except ImportError:
                labdata = None
                print("  (--overlay-data: labdata.py not importable, skipped)")
            if labdata is not None:
                mk = [m for m in labdata.markers() if m["index"] < len(res)]
                if not mk:
                    print("  (--overlay-data: no measured points to place)")
                # One colour, distinct from both model curves: grey would
                # collide with "tube (particles)", near-black with "plate".
                MEAS = "#ea580c"
                FITCOL = "#475569"
                overlay_drawn = True
                # Same drawing routine as the presentation figure, so the
                # two never drift apart in marker shape or spacing.
                _draw_measured(ax2, mk, max_span=0.60)
                has_bounds = any(m["bound"] for m in mk)
                # Nothing is dropped any more -- censored rows are drawn as
                # bounds -- but report anything labdata does hold back, so a
                # missing point never reads as data we never collected.
                if not a.quiet:
                    for item in labdata.unplottable():
                        print(f"  (--overlay-data: not plotted -- {item})")

                # ---- fitted trend line ----------------------------------
                # Censored fit: the bounds stay in as inequalities. See
                # trendfit.py for why dropping them biases the slope.
                if a.fit_line:
                    try:
                        import trendfit
                    except ImportError:
                        trendfit = None
                        print("  (--fit: trendfit.py not importable, skipped)")
                    if trendfit is not None:
                        ft = trendfit.fit_censored()
                        if ft is None:
                            print("  (--fit: too few points to fit)")
                        else:
                            # Draw only across the rounds that were MEASURED.
                            # Extrapolating a 3-round fit out to round 8 would
                            # invite reading it as a prediction, which with
                            # this much scatter it cannot be.
                            fx = np.array(ft["x_levels"], float)
                            fy = 10.0 ** (ft["slope"] * fx + ft["intercept"])
                            ax2.plot(fx, fy, color=FITCOL, lw=2.0,
                                     dashes=(7, 3), zorder=4)
                            fit_short = (f"fitted trend "
                                         f"({1/ft['fold_per_round']:.0f}x"
                                         f"/round)")
                            if not a.quiet:
                                print()
                                print("  FITTED TREND (censored/Tobit, all "
                                      f"{ft['n']} points)")
                                trendfit.report(ft)
        ax2.set_yscale("log"); ax2.set_ylabel("Amp:Kan")
        # A little headroom so the topmost marker does not sit on the frame.
        if has_bounds:
            _lo, _hi = ax2.get_ylim()
            ax2.set_ylim(_lo / 1.5, _hi * 1.5)
        # Say which panel holds what, in the title. Without data overlaid both
        # panels are model output and are labelled so; with data, only the
        # ratio panel mixes the two, and its legend blocks split them.
        ax2.set_title("Ratio across rounds" + (" — model vs measured"
                                               if overlay_drawn
                                               else "  (model)"),
                      loc="left", weight="bold")

        # ---- legend: model and measured kept visibly apart ----------------
        # With data on the panel, one flat list of keys does not say which
        # entries are the model's output and which are plate readings -- the
        # single most important distinction on this figure. Two separate
        # titled blocks say it without a caption. Each measured SHAPE gets its
        # own entry, explained, because a bare triangle tells a reader
        # nothing about which direction it bounds.
        if overlay_drawn:
            from matplotlib.lines import Line2D

            model_handles = [
                Line2D([], [], marker="o", color=c_D, lw=2.4),
                Line2D([], [], marker="o", color=MODEL_PARTICLES, lw=2.0,
                       linestyle="--"),
            ]
            model_labels = ["plate (COLONIES)", "tube (particles)"]

            # Colour names the lineage, so the key lists the constructs.
            cmap = labdata.lineage_colours()
            meas_handles, meas_lineage_labels = [], []
            seen_lin = []
            for m in mk:
                if m["lineage"] in seen_lin:
                    continue
                seen_lin.append(m["lineage"])
                meas_handles.append(
                    Line2D([], [], marker="o", linestyle="-", lw=1.2,
                           markersize=7.0, color=cmap[m["lineage"]],
                           markeredgewidth=1.6))
                meas_lineage_labels.append(m["lineage_label"])

            # No unicode glyph prefixes: the legend handle already draws the
            # shape, so a glyph beside it just doubles the symbol and widens
            # the box enough to cover the curves it sits over.
            meas_labels = list(meas_lineage_labels)
            if fit_short:
                meas_handles.append(
                    Line2D([], [], color=FITCOL, lw=2.0, dashes=(7, 3)))
                meas_labels.append(fit_short)

            # The model curves keep a SMALL in-panel key -- only two entries,
            # and the panel's top-right corner is empty. The measured key,
            # which now names nine lineages, goes below the figure.
            #
            # These were one flat legend with faked bold headers; with this
            # many entries matplotlib's column fill orphaned the "MEASURED"
            # header at the foot of the model column, away from everything it
            # labelled. Two real legends cannot drift like that.
            leg_model = ax2.legend(
                model_handles, model_labels, loc="upper right",
                title="MODEL (simulated)", frameon=True, framealpha=0.93,
                facecolor="white", edgecolor="none", fontsize=8,
                title_fontsize=8.5, borderpad=0.5)
            leg_model.get_title().set_fontweight("bold")
            leg_model.set_zorder(10)

            blank = Line2D([], [], linestyle="none")
            leg = fig.legend([blank] + meas_handles,
                             ["$\\bf{MEASURED}$ (plates)"] + meas_labels,
                             loc="lower center", ncol=4, frameon=False,
                             fontsize=8, bbox_to_anchor=(0.5, -0.01),
                             columnspacing=1.4, handletextpad=0.55)
            leg.set_zorder(10)

        # get_title() defaults to the CENTRE title, which is empty here -- the
        # panel titles are set with loc="left". Reading the wrong one silently
        # replaced this title with just "(model)".
        ax1.set_title(ax1.get_title(loc="left") + "  (model)",
                      loc="left", weight="bold")
        for ax in (ax1, ax2):
            ax.set_xticks(xs); ax.set_xticklabels([r["stage"] for r in res], rotation=45, ha="right")
            if not (ax is ax2 and overlay_drawn):
                ax.legend(frameon=False)
            if ps:
                ps.finish_axes(ax, legend=False)      # legends are drawn above
            else:
                ax.grid(alpha=0.2)
                for s_ in ("top", "right"):
                    ax.spines[s_].set_visible(False)
        fig.suptitle(f"{jumpstart_desc}, eps={a.eps:g} (hybrid a_h={a_h:.2g}), "
                     f"sigma={a.sigma:g}, size effect={a.size_effect:g}", fontsize=10, color="#555")
        # With the keys moved below, reserve space for them; without the
        # overlay there is no outside key and the old layout still applies.
        if overlay_drawn:
            fig.tight_layout(rect=(0, 0.10, 1, 1))
        else:
            fig.tight_layout()
        fig.savefig(a.plot, dpi=200, facecolor="white")
        print(f"\n  wrote {a.plot}")


if __name__ == "__main__":
    main()
