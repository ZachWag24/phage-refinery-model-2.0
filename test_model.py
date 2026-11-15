"""
test_model.py -- checks that guard the whole model against silent breakage.

    python3 test_model.py            everything (about a minute or two)
    python3 test_model.py --fast     skip the slow ones (long chains, command-line runs)
    python3 test_model.py -k phase2  only tests whose name contains "phase2"

Also works under pytest. Standard library plus numpy only.

WHAT KINDS OF TEST ARE HERE
------------------------------
LIMITS       Cases where the right answer is known without running anything:
             interchangeable lineages must not change the ratio; a species that
             is absent can never appear; total export of a dual cell is k2.
STRUCTURE    The equations do what the comments say (equation-level checks that
             do not depend on a whole simulation): exclusion is the same for
             every particle type, each infected cell exports at its own rate.
STATISTICS   passage() really is Poisson.
REGRESSION   Numbers the model gave at the last deliberate change. If one of
             these fails after an edit you did NOT intend to change the results
             with, that is the bug. If the edit was meant to change them, update
             the number here and say why.
DOCS         A few claims in the documentation that were wrong once and are cheap
             to keep right.

A FAILING TEST IS INFORMATION, NOT AN OBSTACLE: read the one-line description.
"""
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import jumpstart                     # noqa: E402
import jumpstart_mechanistic as jm   # noqa: E402
import phase2                        # noqa: E402
import recursion                     # noqa: E402
from passage import passage          # noqa: E402

_TESTS = []


def register(slow=False):
    def deco(fn):
        fn._slow = slow
        _TESTS.append(fn)
        return fn
    return deco


def ok(cond, msg=""):
    assert cond, msg


def close(a, b, rtol=1e-6, what=""):
    assert abs(a - b) <= rtol * max(abs(a), abs(b), 1e-300), f"{what}: {a!r} vs {b!r} (rtol {rtol})"


def p2_args(**o):
    """Positional arguments for phase2.derivatives, with overrides by name."""
    d = dict(k1=3e-11, h=1.0, k2p=3.0, k2h=3.0, k3=np.log(2) / 25, k4=0.6 * np.log(2) / 25,
             K=1e9, eps=1.0, a_h=1.0, sigma=0.05, rho=1.0)
    d.update(o)
    return tuple(d[k] for k in ("k1", "h", "k2p", "k2h", "k3", "k4", "K", "eps", "a_h", "sigma", "rho"))


def _ratio_pp_ph(Y):
    return (Y[-1, 3] + Y[-1, 4]) / (Y[-1, 5] + Y[-1, 6])


# =============================================================================
# JUMPSTART
# =============================================================================
@register()
def test_jumpstart_dual_cell_exports_k2_in_total():
    """A dual cell makes k2 particles/min IN TOTAL, split by phi -- not a full rate for each genome."""
    k2, phi = 3.0, 0.7
    d = jumpstart.derivatives(np.array([0.0, 1e9, 0.0, 0.0]), 3e-11, 1.0, k2, phi, 0.03, 0.02, 1e9)
    close(d[2] + d[3], k2 * 1e9, 1e-9, "total export")
    close(d[2], phi * k2 * 1e9, 1e-9, "phagemid share")


@register()
def test_jumpstart_ratio_follows_phi():
    """The Jumpstart's Amp:Kan is phi/(1-phi): the split is an explicit input, nothing hidden."""
    for phi in (0.5, 0.75, 0.9):
        t, Y = jumpstart.run(phi=phi)
        close(Y[-1, 2] / Y[-1, 3], phi / (1 - phi), 0.01, f"phi={phi}")


@register()
def test_jumpstart_default_phi_is_the_measured_split():
    """With no phi given the split is the MEASURED 4.88:1, not a derived number."""
    t, Y = jumpstart.run()
    close(Y[-1, 2] / Y[-1, 3], 4.88, 0.01)
    close(jumpstart.PHI_MEASURED, 4.88 / 5.88, 1e-12)


@register()
def test_jumpstart_negative_phi_selects_the_length_hypothesis():
    """phi<0 still gives the length-derived split (2.29:1, or 1:1 when Lp=Lh)."""
    close(jumpstart.run(phi=-1)[1][-1, 2] / jumpstart.run(phi=-1)[1][-1, 3], 8.7 / 3.8, 0.01)
    close(jumpstart.run(phi=-1, Lp=8.7)[1][-1, 2] / jumpstart.run(phi=-1, Lp=8.7)[1][-1, 3], 1.0, 0.01)


@register()
def test_recursion_length_phi_actually_uses_the_length_hypothesis():
    """recursion.run_recursion(phi="length") must give jumpstart.py's length-derived split, not the measured default."""
    R = recursion.run_recursion(rounds=0, phi="length", helper_edge=1.0)
    close(R[0]["ratio"], 8.7 / 3.8, 0.01, "recursion phi='length' vs jumpstart phi=-1")
    ok(abs(R[0]["ratio"] - 4.88) > 0.5, "must not have silently fallen back to the measured default")


@register()
def test_jumpstart_and_recursion_agree_on_the_measured_jumpstart():
    """Both files start from the same measurement: recursion's only difference is correcting for the helper edge."""
    R = recursion.run_recursion(rounds=0)
    close(R[0]["colony_ratio"], jumpstart.PHI_MEASURED / (1 - jumpstart.PHI_MEASURED), 2e-3,
          "recursion's Jumpstart plate ratio vs jumpstart.py's measured split")


@register()
def test_jumpstart_without_helper_makes_nothing():
    """No helper (moi=0) means no dual cell and no particle of either kind, ever."""
    t, Y = jumpstart.run(moi=0.0)
    ok(Y[:, 2:].max() == 0.0)


@register()
def test_jumpstart_doubling_export_doubles_output_not_ratio():
    """Output scales with the export rate; the ratio is set by phi, not by how fast cells export."""
    a = jumpstart.run(phi=0.8)[1][-1]
    b = jumpstart.run(phi=0.8, k2=6.0)[1][-1]
    ok(1.7 < (b[2] + b[3]) / (a[2] + a[3]) < 2.3, "total output should roughly double")
    close(a[2] / a[3], b[2] / b[3], 0.02, "ratio")


@register()
def test_jumpstart_kan_removes_only_uninfected_cells():
    """Kanamycin kills phagemid-only cells and leaves dual cells alone."""
    t, Y = jumpstart.run(kan_at=75.0)
    i = int(np.argmax(t >= 75.0)) + 2
    ok(Y[i, 0] == 0.0 and Y[-1, 0] == 0.0, "phagemid-only cells should be gone after Kan")
    ok(Y[i, 1] > 0.0 and Y[-1, 1] > 0.9e9, "dual cells should be untouched and fill the culture")
    ok(jumpstart.run(kan_at=None)[1][i, 0] > 0.0, "without Kan they should still be there")


@register()
def test_jumpstart_bounds_and_convergence():
    """Cells never exceed K, nothing goes negative, and the answer does not depend on the step size."""
    t, Y = jumpstart.run()
    ok(Y.min() >= 0 and (Y[:, 0] + Y[:, 1]).max() <= 1e9 * (1 + 1e-3))
    r = [jumpstart.run(steps=s)[1][-1] for s in (1000, 4000, 16000)]
    close(r[1][2] / r[1][3], r[2][2] / r[2][3], 1e-3, "ratio vs step size")


# =============================================================================
# JUMPSTART_MECHANISTIC (MODEL_REPORT.md section 9, built fresh)
# =============================================================================
@register()
def test_mechanistic_collapses_to_alpha_times_nres_when_rates_equal():
    """IDENTIFIABILITY: with rp==rh (replication_ratio=1), output ratio == alpha*Nres exactly, regardless of Ncap/T."""
    for Nres, alpha in ((5.0, 1.0), (20.0, 2.0), (50.0, 0.7)):
        t, Y, _ = jm.run(Nres=Nres, replication_ratio=1.0, alpha=alpha, Ncap=300.0)
        close(Y[-1, 2] / Y[-1, 3], alpha * Nres, 1e-2, f"Nres={Nres} alpha={alpha}")
    # must ALSO hold at a different Ncap -- the collapse does not depend on capacity
    t, Y, _ = jm.run(Nres=10.0, replication_ratio=1.0, alpha=1.0, Ncap=2000.0)
    close(Y[-1, 2] / Y[-1, 3], 10.0, 1e-2, "collapse independent of Ncap")


@register()
def test_mechanistic_ratio_rises_with_culture_length():
    """When the phagemid's origin genuinely replicates faster (replication_ratio>1), longer incubation should favour it more, not less."""
    kw = dict(Nres=10.0, replication_ratio=1.1, alpha=1.0, Ncap=300.0)
    r6 = jm.run(T=6 * 60.0, **kw)[1][-1]
    r18 = jm.run(T=18 * 60.0, **kw)[1][-1]
    r36 = jm.run(T=36 * 60.0, **kw)[1][-1]
    ratios = [p / h for p, h in ((r6[2], r6[3]), (r18[2], r18[3]), (r36[2], r36[3]))]
    ok(ratios[0] < ratios[1] < ratios[2], f"ratio should rise monotonically with T: {ratios}")


@register()
def test_mechanistic_ratio_rises_with_ncap():
    """More shared replication capacity gives the faster replicator more room to pull ahead."""
    kw = dict(Nres=10.0, replication_ratio=1.1, alpha=1.0, T=1080.0)
    ratios = [jm.run(Ncap=nc, **kw)[1][-1, 2] / jm.run(Ncap=nc, **kw)[1][-1, 3]
              for nc in (50.0, 300.0, 2000.0)]
    ok(ratios[0] < ratios[1] < ratios[2], f"ratio should rise monotonically with Ncap: {ratios}")


@register(slow=True)
def test_mechanistic_calibration_round_trips():
    """resolve_alpha must find an alpha whose OWN simulated output reproduces the target ratio, holding the real Nres fixed."""
    for target in (2.0, 4.88, 20.0):
        t, Y, alpha = jm.run(jump_ratio=target)
        close(Y[-1, 2] / Y[-1, 3], target, 1e-3, f"target={target}")
        ok(alpha > 0, "calibrated alpha must be positive")


@register()
def test_mechanistic_predictive_mode_does_not_recalibrate():
    """Given an explicit alpha, the model must NOT silently adjust it to match jump_ratio -- that would defeat the point of a prediction."""
    t, Y, alpha_used = jm.run(alpha=1.0, jump_ratio=4.88)
    ok(alpha_used == 1.0, "explicit alpha must be used as-is, not recalibrated")
    ratio = Y[-1, 2] / Y[-1, 3]
    ok(abs(ratio - 4.88) > 1.0, "alpha=1 at the real Nres=18 should NOT reproduce jump_ratio=4.88 by coincidence")


@register(slow=True)
def test_mechanistic_sensitivity_tables_hold_alpha_fixed():
    """sensitivity_tables must calibrate alpha ONCE at baseline, not re-fit it at every sweep point (which would trivially flatten every table to jump_ratio)."""
    t_table, ncap_table, alpha = jm.sensitivity_tables(jump_ratio=4.88, Ncap=300.0, T=1080.0)
    close(dict(t_table)[18], 4.88, 1e-2, "the calibration point itself should hit the target")
    ok(dict(t_table)[6] != dict(t_table)[36], "other points must NOT all collapse onto jump_ratio")
    close(dict(ncap_table)[300], 4.88, 1e-2, "the calibration point itself should hit the target")
    ok(alpha > 0, "calibrated alpha must be positive")


@register()
def test_mechanistic_no_helper_calibration_fails_clearly():
    """moi=0 means Ph stays exactly 0 forever regardless of alpha -- calibration has nothing to solve and must say so, not crash on a bisection-bracket error."""
    try:
        jm.run(moi=0.0)
        ok(False, "expected a ValueError")
    except ValueError as e:
        ok("moi=0" in str(e) or "no helper" in str(e).lower(), str(e))


@register(slow=True)
def test_mechanistic_default_regression():
    """REGRESSION: calibrated alpha (at the real Nres=18) and the two sensitivity tables at the module's documented defaults."""
    t, Y, alpha = jm.run()
    close(alpha, 0.2121, 2e-3, "calibrated alpha")
    t_table, ncap_table, _ = jm.sensitivity_tables()
    want_t = {6: 4.6442, 18: 4.8800, 36: 4.9119}
    want_ncap = {50: 4.4794, 300: 4.8800, 2000: 5.3024}
    for hours, r in t_table:
        close(r, want_t[hours], 2e-3, f"T={hours}h")
    for ncap, r in ncap_table:
        close(r, want_ncap[ncap], 2e-3, f"Ncap={ncap}")


# =============================================================================
# PASSAGE
# =============================================================================
@register()
def test_passage_high_counts_is_plain_division():
    """At huge counts passage() is just the 0.5 mL into 5.5 mL dilution."""
    close(passage(1e12, rng=np.random.default_rng(0))[0], 1e12 * 0.5 / 5.5, 1e-9)


@register()
def test_passage_is_poisson():
    """The number of particles pipetted has Poisson mean and spread."""
    rng = np.random.default_rng(1)
    n = np.array([passage(2e4, rng=rng)[0] for _ in range(20000)]) * 5.5     # 1e4 particles transferred
    close(n.mean(), 1e4, 0.01, "mean")
    close(n.std(), 100.0, 0.05, "spread")


@register()
def test_passage_extinction_probability_is_exp_minus_mean():
    """A pool of 4 particles/mL is missed entirely in the 0.5 mL pipetted about 13.5% of the time."""
    rng = np.random.default_rng(2)
    z = np.mean([passage(4.0, rng=rng)[0] == 0 for _ in range(40000)])
    close(z, np.exp(-2.0), 0.05, "P(zero)")


@register()
def test_passage_pools_are_sampled_independently():
    """Each pool is its own set of particles: draws from two pools are uncorrelated."""
    rng = np.random.default_rng(3)
    a = np.array([passage(4e3, 4e3, rng=rng) for _ in range(20000)])
    ok(abs(np.corrcoef(a[:, 0], a[:, 1])[0, 1]) < 0.03)


@register()
def test_passage_reproducible_and_takes_four_pools():
    """Same seed, same draw; four pools in, four out."""
    r1 = passage(50.0, 60.0, 70.0, 80.0, rng=np.random.default_rng(5))
    r2 = passage(50.0, 60.0, 70.0, 80.0, rng=np.random.default_rng(5))
    ok(r1 == r2 and len(r1) == 4)
    ok(passage(50.0, 60.0, 70.0, 80.0, rng=np.random.default_rng(6)) != r1)


# =============================================================================
# PHASE 2
# =============================================================================
@register()
def test_phase2_interchangeable_lineages_preserve_any_input_ratio():
    """With no difference between the lineages, whatever ratio goes in comes out -- including mixed coat states."""
    for rin in (0.1, 1, 10, 100):
        t, Y = phase2.run(rin * 1e11, 1e11)
        close(_ratio_pp_ph(Y), rin, 1e-6, f"{rin}:1")
    t, Y = phase2.run(3e10, 2e10, Pt_in=5e10, Hh_in=5e10)
    close(_ratio_pp_ph(Y), 8 / 7, 1e-6, "mixed-coat input 8:7")


@register()
def test_phase2_absent_species_never_appears():
    """No helper in means no helper ever; no phagemid in means no phagemid ever."""
    t, Y = phase2.run(2e11, 0.0)
    ok(Y[:, 5].max() == 0 and Y[:, 6].max() == 0 and Y[:, 2].max() == 0)
    t, Y = phase2.run(0.0, 4e10)
    ok(Y[:, 3].max() == 0 and Y[:, 4].max() == 0 and Y[:, 1].max() == 0)


@register()
def test_phase2_infected_cells_export_at_their_own_rate():
    """An Ip cell makes fully modified phagemid at k2p; an Ih cell makes hybrid helper at k2h; nothing else."""
    a = p2_args(k2p=4.0, k2h=3.0)
    d = phase2.derivatives(np.array([0, 1e9, 0, 0, 0, 0, 0.0]), *a)
    close(d[4], 4.0e9, 1e-9, "Ip -> Pt")
    ok(d[6] == 0.0 and d[3] == 0.0 and d[5] == 0.0)
    d = phase2.derivatives(np.array([0, 0, 1e9, 0, 0, 0, 0.0]), *a)
    close(d[6], 3.0e9, 1e-9, "Ih -> Hh")
    ok(d[4] == 0.0 and d[3] == 0.0 and d[5] == 0.0)


@register()
def test_phase2_exclusion_is_the_same_for_every_particle_type():
    """An infected cell adsorbs ANY particle at sigma times a naive cell's rate, and per-particle weights are as documented."""
    sigma, eps, a_h_val, edge = 0.05, 0.3, 0.6, 1.5
    a = p2_args(eps=eps, a_h=a_h_val, sigma=sigma, rho=1 / edge, k2p=0.0, k2h=0.0)
    rate = {}
    for cls, cn in ((0, "N"), (1, "Ip"), (2, "Ih")):
        for pool, pn in ((3, "Pw"), (4, "Pt"), (5, "Hw"), (6, "Hh")):
            y = np.zeros(7)
            y[cls], y[pool] = 1e9, 1e6
            rate[(cn, pn)] = -phase2.derivatives(y, *a)[pool] / 1e6
    for pn in ("Pw", "Pt", "Hw", "Hh"):
        close(rate[("Ip", pn)] / rate[("N", pn)], sigma, 1e-9, f"Ip/{pn}")
        close(rate[("Ih", pn)] / rate[("N", pn)], sigma, 1e-9, f"Ih/{pn}")
    close(rate[("N", "Pw")] / rate[("N", "Hw")], 1 / edge, 1e-3, "helper edge")
    close(rate[("N", "Pt")] / rate[("N", "Pw")], eps, 1e-3, "coat gap (phagemid)")
    close(rate[("N", "Hh")] / rate[("N", "Hw")], a_h_val, 1e-3, "hybrid coat (helper)")


@register()
def test_phase2_hybrid_advantage_over_fully_modified_is_as_documented():
    """One hybrid helper particle counts ~21.9x on a plate against one fully modified phagemid particle (eps .05, defaults m13ko7_copies=20/hp17_copies=8, edge 1.5)."""
    eps, edge = 0.05, 1.5
    a_h = phase2.hybrid_infectivity(eps)
    close(a_h, 0.72857, 1e-4, "a_h at defaults")
    amp, kan = phase2.plate_counts(0.0, 1.0, 0.0, 1.0, eps, a_h, edge)
    close(kan / amp, edge * a_h / eps, 1e-12)
    close(kan / amp, 21.857, 1e-3)


@register()
def test_phase2_bigger_infectivity_gap_means_lower_plate_ratio():
    """Lowering eps, or raising the helper edge, always lowers the plate Amp:Kan."""
    def plate(**kw):
        Y = phase2.run(2e11, 4e10, **kw)[1]
        eps = kw.get("eps", 1.0)
        a_h = phase2.hybrid_infectivity(eps, kw.get("m13ko7_copies", 20.0), kw.get("hp17_copies", 8.0))
        a, k = phase2.plate_counts(*Y[-1, 3:7], eps, a_h, kw.get("helper_edge", 1.0))
        return a / k
    ok(plate(eps=1.0) > plate(eps=0.2) > plate(eps=0.05))
    ok(plate(helper_edge=1.0) > plate(helper_edge=1.5) > plate(helper_edge=3.0))


@register()
def test_phase2_size_effect_switch():
    """size_effect=0 is exactly 'same length'; size_effect=1 makes the phagemid lineage out-produce the helper."""
    a = phase2.run(2e11, 4e10, size_effect=0.0)[1]
    b = phase2.run(2e11, 4e10, Lp=8.7)[1]
    ok(np.allclose(a, b, rtol=1e-12, atol=0), "size_effect=0 must equal Lp=Lh")
    t, Y = phase2.run(1e11, 1e11, size_effect=1.0)
    ok(_ratio_pp_ph(Y) > 1.5, "a 2.29x production advantage should show")


@register()
def test_phase2_origin_penalty_switch():
    """origin_penalty=1 is exactly 'no penalty'; a penalty < 1 slows Ih's own export rate, favouring the phagemid lineage."""
    a = phase2.run(2e11, 4e10, origin_penalty=1.0)[1]
    b = phase2.run(2e11, 4e10)[1]   # phase2.py's own default is neutral (1.0) too
    ok(np.allclose(a, b, rtol=1e-12, atol=0), "origin_penalty defaults to neutral in phase2.py itself")
    t, Y = phase2.run(2e11, 4e10, eps=1.0, helper_edge=1.0, origin_penalty=0.9)
    ok(_ratio_pp_ph(Y) > 1.0, "a penalty on Ih's export rate, with everything else neutral, should favour the phagemid lineage")


@register()
def test_phase2_bounds_and_convergence():
    """Cells never exceed K, nothing goes negative, and the answer does not depend on the step size."""
    t, Y = phase2.run(2e11, 4e10, eps=0.05, helper_edge=1.5)
    ok(Y.min() >= 0 and Y[:, :3].sum(axis=1).max() <= 1e9 * (1 + 1e-3))

    a_h = phase2.hybrid_infectivity(0.05)

    def plate(steps):
        Yc = phase2.run(2e11, 4e10, eps=0.05, helper_edge=1.5, steps=steps)[1]
        a, k = phase2.plate_counts(*Yc[-1, 3:7], 0.05, a_h, 1.5)
        return a / k
    close(plate(4000), plate(16000), 1e-4, "step size")


# =============================================================================
# RECURSION
# =============================================================================
@register()
def test_recursion_jumpstart_reports_the_measured_plate_ratio():
    """The Jumpstart starts at the measured PLATE ratio; its particle ratio is that times the helper edge (default 1.0, no longer an active effect -- see MODEL_REPORT.md section 7.10)."""
    R = recursion.run_recursion(rounds=0)
    close(R[0]["colony_ratio"], 4.88, 1e-3, "plate")
    close(R[0]["ratio"], 4.88 * 1.0, 1e-3, "particle")
    R = recursion.run_recursion(rounds=0, helper_edge=2.0, jump_ratio=8.0)
    close(R[0]["colony_ratio"], 8.0, 1e-3)
    close(R[0]["ratio"], 16.0, 1e-3)


@register()
def test_phase2_cli_default_inputs_match_what_recursion_hands_it():
    """phase2.py's built-in -P/-H defaults must not drift from the Jumpstart output / 11 that recursion.py delivers."""
    R = recursion.run_recursion(rounds=0)
    a = phase2.build_parser().parse_args([])
    close(a.Pw_in, R[0]["Pp"] / 11, 0.02, "Pw_in default")
    close(a.Hw_in, R[0]["Ph"] / 11, 0.02, "Hw_in default")


@register(slow=True)
def test_recursion_default_chain_regression():
    """REGRESSION: the default chain's plate Amp:Kan at the last deliberate change -- m13ko7_copies 200 -> 20 and origin_penalty 0.9 -> 1/1.05 on 2026-11-14 (aligned to the published wiki)."""
    want = [4.8796, 0.86987, 0.087116, 0.010451, 0.0015503, 0.00024109]
    for i, (x, w) in enumerate(zip(recursion.run_recursion(rounds=5), want)):
        close(x["colony_ratio"], w, 2e-3, f"stage {i}")


@register(slow=True)
def test_recursion_edge_one_regression():
    """REGRESSION: with no helper edge and phi=0.83 the chain is the earlier four-pool model's (values reset for m13ko7_copies 200 -> 20 and origin_penalty 0.9 -> 1/1.05, 2026-11-14)."""
    want = [4.8820, 0.87028, 0.087153, 0.010455, 0.0015510]
    R = recursion.run_recursion(rounds=4, phi=0.83, helper_edge=1.0)
    for i, (x, w) in enumerate(zip(R, want)):
        close(x["colony_ratio"], w, 2e-3, f"stage {i}")


@register(slow=True)
def test_recursion_no_advantage_anywhere_never_moves_the_ratio():
    """With no coat gap, no edge, and no origin penalty the ratio stays exactly where the Jumpstart left it."""
    q = [x["colony_ratio"] for x in recursion.run_recursion(rounds=5, eps=1.0, helper_edge=1.0, origin_penalty=1.0)]
    ok(np.allclose(q, q[0], rtol=1e-6), str(q))


@register(slow=True)
def test_recursion_default_trend():
    """The documented trend: plate ratio falls every stage, is below 1 by round 2, and the helper persists."""
    R = recursion.run_recursion(rounds=8)
    col = [x["colony_ratio"] for x in R]
    tube = [x["ratio"] for x in R]
    ok(all(np.diff(col) < 0), "plate ratio must fall at every stage")
    ok(col[2] < 1.0, "plate ratio below 1 by round 2")
    ok(min(x["Ph"] for x in R) > 3e11, "helper must persist")
    ok(all(c <= t_ + 1e-12 for c, t_ in zip(col[1:], tube[1:])), "plate ratio cannot exceed tube ratio here")


@register(slow=True)
def test_recursion_mechanistic_stays_in_the_same_ballpark_as_fixed_phi():
    """--mechanistic (calibrated) should land close to the fixed-phi model at the same jump_ratio -- that was the whole point of wiring it in."""
    fixed = recursion.run_recursion(rounds=3, jump_ratio=4.88)
    mech = recursion.run_recursion(rounds=3, jump_ratio=4.88, mechanistic=True)
    for f, m in zip(fixed, mech):
        rel_diff = abs(f["colony_ratio"] - m["colony_ratio"]) / f["colony_ratio"]
        ok(rel_diff < 0.02, f"{f['stage']}: fixed={f['colony_ratio']:.4f} mech={m['colony_ratio']:.4f} ({rel_diff*100:.1f}% apart)")


@register()
def test_recursion_mechanistic_predictive_mode_need_not_match_jump_ratio():
    """Giving mech_alpha explicitly must produce a genuine prediction, not silently recalibrate back to jump_ratio."""
    res = recursion.run_recursion(rounds=0, jump_ratio=4.88, mechanistic=True, mech_alpha=1.0)
    ok(abs(res[0]["colony_ratio"] - 4.88) > 1.0,
       "alpha=1 at the real Nres=18 should generally not reproduce jump_ratio by coincidence")


@register(slow=True)
def test_recursion_mechanistic_band_still_shows_real_spread():
    """--band --mechanistic: each of low/center/high must calibrate its OWN alpha (Nres is real and shared), not share one -- otherwise the band collapses to a single line."""
    band = recursion.run_band(low=3.0, high=16.0, center=4.88, rounds=0, mechanistic=True)
    row = band[0]
    ok(row["colony_low"] < row["colony_center"] < row["colony_high"],
       f"band should still show real spread under --mechanistic: {row}")


@register(slow=True)
def test_recursion_mechanistic_default_regression():
    """REGRESSION: the mechanistic chain at defaults, at the last deliberate change -- m13ko7_copies 200 -> 20 and origin_penalty 0.9 -> 1/1.05 on 2026-11-14 (aligned to the published wiki)."""
    res = recursion.run_recursion(rounds=3, mechanistic=True)
    want = [4.8800, 0.87625, 0.087870, 0.010531]
    for r, w in zip(res, want):
        close(r["colony_ratio"], w, 5e-3, r["stage"])


@register()
def test_recursion_band_uses_the_same_seed_across_scenarios():
    """run_band's low/center/high must share one seed -- any difference between them must come from jump_ratio alone."""
    band = recursion.run_band(low=3.0, high=16.0, center=4.88, rounds=3, seed=7)
    # center should exactly match a plain run_recursion call with the same seed and jump_ratio
    plain = recursion.run_recursion(rounds=3, seed=7, jump_ratio=4.88)
    for b, p in zip(band, plain):
        close(b["colony_center"], p["colony_ratio"], 1e-12, "band center must equal a plain run with the same seed")


@register()
def test_recursion_band_is_ordered_low_to_high():
    """The Jumpstart's own plate ratio (least perturbed by downstream dynamics) must satisfy low <= center <= high."""
    band = recursion.run_band(low=3.0, high=16.0, center=4.88, rounds=0)
    row = band[0]
    ok(row["colony_low"] <= row["colony_center"] <= row["colony_high"],
       f"low/center/high out of order: {row}")


@register(slow=True)
def test_recursion_band_pads_extinct_scenarios_with_nan_not_truncation():
    """A scenario that goes extinct one round earlier than another must contribute NaN at that one row, not cut the whole table short to the shortest scenario's length."""
    band = recursion.run_band(low=3.0, high=16.0, center=4.88, rounds=40, seed=1)
    lengths = {jr: len(recursion.run_recursion(rounds=40, seed=1, jump_ratio=jr))
               for jr in (3.0, 4.88, 16.0)}
    # here the three scenarios go extinct within a round of each other, not all at once --
    # the table's length must be the LONGEST scenario's, not the shortest's, with whichever
    # scenario(s) stopped earlier NaN-padded for the remaining row(s). Which scenario is
    # shortest is incidental to current parameters, not something to hard-code.
    ok(len(lengths.values()) > 1 and len(set(lengths.values())) > 1,
       f"expected staggered extinction across scenarios, got {lengths}")
    ok(len(band) == max(lengths.values()),
       f"table length {len(band)} should equal the longest scenario's length {max(lengths.values())}")
    ok(any(np.isnan(r["colony_low"]) or np.isnan(r["colony_center"]) or np.isnan(r["colony_high"])
           for r in band),
       "at least one scenario that stopped earlier must show NaN at the row(s) it has no data for")


@register()
def test_fmt_ratio_distinguishes_nan_from_true_infinity():
    """NaN (no data here) and a genuine infinite ratio (extinct helper) must render differently -- conflating them in the band table would be misleading."""
    ok(recursion._fmt_ratio(float("nan")) != recursion._fmt_ratio(float("inf")))
    ok(recursion._fmt_ratio(float("inf")) == "inf")


@register(slow=True)
def test_recursion_band_default_regression():
    """REGRESSION: the default --band chain (jump_ratio 3-16, center 4.88) at the last deliberate change -- m13ko7_copies 200 -> 20 and origin_penalty 0.9 -> 1/1.05 on 2026-11-14 (aligned to the published wiki)."""
    band = recursion.run_band(rounds=5)
    want = [
        (2.9999, 4.8796, 15.993),
        (0.54141, 0.86987, 2.7812),
        (0.057416, 0.087116, 0.25603),
        (0.0073674, 0.010451, 0.025556),
        (0.0011100, 0.0015503, 0.0035400),
        (0.00017302, 0.00024109, 0.00054438),
    ]
    for i, (row, (lo, ct, hi)) in enumerate(zip(band, want)):
        close(row["colony_low"], lo, 2e-3, f"stage {i} low")
        close(row["colony_center"], ct, 2e-3, f"stage {i} center")
        close(row["colony_high"], hi, 2e-3, f"stage {i} high")


@register(slow=True)
def test_recursion_only_the_phagemid_goes_extinct():
    """In a long chain the phagemid, and only the phagemid, is eventually lost (Poisson sampling in passage)."""
    R = recursion.run_recursion(rounds=40, seed=1)
    ext = [r for r in R if r.get("extinct")]
    ok(len(ext) == 1 and ext[0]["extinct"] == "pHK002.2", str([r.get("extinct") for r in R]))
    n = int(ext[0]["stage"].split()[1])
    ok(15 <= n <= 35, f"extinction at round {n}")


@register(slow=True)
def test_recursion_is_deterministic_at_current_titers():
    """At these titers passage() noise is ~1e-6, so different seeds give the same early rounds."""
    a = recursion.run_recursion(rounds=4, seed=1)
    b = recursion.run_recursion(rounds=4, seed=2)
    for x, y in zip(a, b):
        close(x["colony_ratio"], y["colony_ratio"], 1e-5)


# =============================================================================
# DOCS AND FILES
# =============================================================================
def _read(fn):
    with open(os.path.join(HERE, fn)) as f:
        return f.read()


@register()
def test_docs_do_not_make_claims_that_were_wrong():
    """Claims that were wrong once: saturation 'left out' of the Jumpstart, old pool names in passage.py."""
    j = _read("jumpstart.py")
    ok("no saturation" not in j, "jumpstart.py says saturation is absent; it is implemented")
    ok("NOT left out: Holling" in j)
    ok("(Pt, Pm, Ph)" not in _read("passage.py"), "passage.py still names the old pools")


@register()
def test_all_four_plotting_scripts_use_the_shared_style():
    """plotstyle.py exists so one edit changes every figure; it was once imported by nothing."""
    ok(os.path.exists(os.path.join(HERE, "plotstyle.py")))
    for f in ("jumpstart.py", "phase2.py", "recursion.py", "jumpstart_mechanistic.py"):
        ok("import plotstyle" in _read(f), f)


@register(slow=True)
def test_command_line_scripts_run_and_write_outputs():
    """Every script runs with default arguments and writes its csv and figure."""
    with tempfile.TemporaryDirectory() as tmp:
        runs = [["jumpstart.py", "-q", "--csv", "j.csv", "--plot", "j.png"],
                ["phase2.py", "-q", "-e", "0.05", "--csv", "p.csv", "--plot", "p.png"],
                ["recursion.py", "-q", "-r", "2", "--csv", "r.csv", "--plot", "r.png"],
                ["jumpstart_mechanistic.py", "-q", "--csv", "m.csv", "--plot", "m.png"],
                ["passage.py"]]
        for args in runs:
            r = subprocess.run([sys.executable, os.path.join(HERE, args[0])] + args[1:],
                               cwd=tmp, capture_output=True, text=True)
            ok(r.returncode == 0, f"{args[0]} failed: {r.stderr[-300:]}")
        for fn in ("j.csv", "p.csv", "r.csv", "m.csv"):         # a 2-round recursion csv is only ~250 bytes
            ok(os.path.getsize(os.path.join(tmp, fn)) > 100, f"{fn} missing or empty")
        for fn in ("j.png", "p.png", "r.png", "m.png"):
            ok(os.path.getsize(os.path.join(tmp, fn)) > 5000, f"{fn} missing or empty")


# =============================================================================
# RUNNER
# =============================================================================
def main(argv):
    fast = "--fast" in argv
    key = argv[argv.index("-k") + 1] if "-k" in argv else None
    chosen = [t for t in _TESTS if (not fast or not t._slow) and (key is None or key in t.__name__)]
    failed, t0 = [], time.time()
    for fn in chosen:
        t1 = time.time()
        try:
            fn()
            status = "PASS"
        except Exception as e:                       # noqa: BLE001 -- report any failure, keep going
            status = "FAIL"
            failed.append((fn.__name__, e))
        first = (fn.__doc__ or "").strip().splitlines()[0] if fn.__doc__ else ""
        print(f"  [{status}] {fn.__name__:62s} {time.time() - t1:5.1f}s  {first if status == 'FAIL' else ''}")
        if status == "FAIL":
            print(f"         -> {failed[-1][1]}")
    print(f"\n{len(chosen) - len(failed)} of {len(chosen)} passed in {time.time() - t0:.0f}s"
          + ("" if not failed else f"   FAILED: {[n for n, _ in failed]}"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
