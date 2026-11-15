"""
Phase 2: one overnight round in HP17KO7-bearing cells -- TWO LINEAGES.

THE PICTURE
-----------
Every Phase 2 cell carries HP17KO7, the helper plasmid that supplies the phage
assembly machinery AND the engineered (tagged) coat. Two kinds of particle
start two lineages, and each lineage makes its own kind:

    N --pHK002.2 particle--> Ip --> makes pHK002.2 phage, FULLY MODIFIED coat  (Pt)
    N --M13KO7   particle--> Ih --> makes M13KO7   phage, HYBRID coat          (Hh)

    N   naive cell (HP17KO7 only)
    Ip  cell infected by a pHK002.2 particle
    Ih  cell infected by an M13KO7 particle
    Pt  pHK002.2 particles made in Phase 2. The phagemid has no genes of its
        own, so its particles are assembled from HP17KO7's coat alone: fully
        modified.
    Hh  M13KO7 particles made in Phase 2 ("H" for helper, not "phagemid").
        M13KO7 brings its own wild-type coat proteins, so its particles wear
        an M13KO7/HP17 HYBRID coat.

Two more pools exist, but only briefly. The particles coming out of the
JUMPSTART were made without HP17KO7, so they wear plain WILD-TYPE coat:

    Pw  pHK002.2 particles from the Jumpstart   (wild-type coat)
    Hw  M13KO7   particles from the Jumpstart   (wild-type coat)

They infect at full strength (infectivity 1) and simply carry over, diluted 1:11
at every passage; nothing in Phase 2 makes more of them. Each wild-type
particle starts the same lineage as its Phase 2 counterpart: Pw makes Ip, Hw
makes Ih.

WHY THE AMP:KAN RATIO FALLS
------------------------------
The two lineages differ in the thing that matters here: how well their
particles infect. The hybrid coat infects better than the fully modified
coat. `eps` is the fully modified coat's infectivity relative to wild-type.
The hybrid coat's infectivity, a_h, is now derived from real coat
composition rather than fit as an exponent -- see THE HYBRID COAT'S
INFECTIVITY, DERIVED below. Each round the two kinds of particle compete
for the same pool of fresh cells. The better infector captures more of
them, so more cells make ITS kind, and it enters the next round in greater
proportion. The plate readout compounds this: Amp colonies count Pt
weighted by eps, Kan colonies count Ph weighted by a_h.

With NO difference in infectivity and NO difference in production (eps=1,
size_effect=0) the two lineages are interchangeable and the ratio never
changes -- checked in the tests, and it is what the defaults here give.

THE HYBRID COAT'S INFECTIVITY, DERIVED
-------------------------------------------
An Ih cell already carries HP17KO7 (resident, supplying MODIFIED pIII and
pVIII) and gets infected by an incoming M13KO7 particle (bringing its OWN
wild-type pIII and pVIII genes). Both gene sets get expressed, and Hh
particles assemble from that mixed pool -- genuinely stochastic at the
single-particle level, since a virion carries only ~5 copies of pIII (the
tip protein that actually does receptor binding, so functionally the one
that sets infectivity) versus ~2700 copies of pVIII (the bulk structural
coat, whose mixing barely matters for infectivity at that copy number).
This file tracks the POPULATION-AVERAGE wild-type fraction, the same kind
of simplification made everywhere else here (no per-particle distribution,
same as no eclipse delay, no cell-to-cell variation).

Given roughly equal promoter strength for the two pIII genes (relayed by
the user), pIII output per gene copy is about equal, so the wild-type
share of pIII protein in an Ih cell is set by GENE DOSAGE alone:

    wt_frac = m13ko7_copies / (m13ko7_copies + hp17_copies)
    a_h     = wt_frac + (1 - wt_frac) * eps          (linear mix)

m13ko7_copies defaults to 20 -- M13KO7's OWN copy number via its p15A
origin, NOT wild-type M13's f1 rolling-circle number. The two are easy to
conflate and the difference is an order of magnitude. A wild-type M13
infection drives its replicative form to ~130-300 copies off an intact f1
origin. M13KO7's f1 origin is the thing that is interrupted -- by the KanR
insertion, which is exactly what makes it a helper, since the damaged origin
causes it to preferentially replicate any co-resident wild-type f1 origin
(a phagemid's) over its own. That interruption is compensated by a separate
p15A origin, which replicates M13KO7 like an ordinary low-copy plasmid.
Published p15A estimates cluster around 10-30 copies/cell.

WHY 20 AND NOT 10. 20 is the high end of that published range, chosen because
the measured round-1 plate data favours it: at 10 the chain's round 1 lands
just above 1 (1.06), at 20 it lands below (0.870), and every measured round 1
is below 1. Both values are defensible from the literature alone -- the data
is what picks between them, which also makes this a falsifiable choice rather
than a free one.

HISTORY, worth knowing because the figure moved twice. It was 200 until
2026-11-14 (wild-type M13's number, simply the wrong quantity for this
construct), then 10 (the low end of the real p15A range), then 20. All three
rest on the same underlying fact -- M13KO7 replicates via p15A, not the
damaged f1 origin -- and only the point estimate within that range changed.
A direct qPCR of M13KO7's copy number in this host background is the one
measurement that would settle it.

hp17_copies defaults to 8 -- HP17's Benchling map shows a Rep101 feature,
the replication protein for the pSC101 origin, one of the lowest-copy
origins in routine use. Published pSC101 copy numbers range from ~5
copies/cell (the original characterization) to ~10-12 copies/cell (an
independent single-cell measurement); 8 sits between them.

At these defaults, wt_frac ~ 0.71: an Ih cell's pIII pool is still mostly
wild-type M13KO7 protein rather than HP17's modified version, so a_h (0.729)
sits clearly above the midpoint between eps and full wild-type strength --
though not close to 1, the way the old 200-copy default made it look
(wt_frac 0.96, a_h 0.963).

THIS INPUT IS NOT ONE THE MODEL IS INSENSITIVE TO, and an earlier version of
the docs wrongly implied it was. That was true at 200, where M13KO7's copy
number swamped HP17's whatever its exact value. Now that the two are
comparable (~20 vs ~8) BOTH genuinely move the answer: across the published
10-30 range the round-2 plate ratio moves 1.73x, and hp17_copies across 5-12
moves it 1.47x. That is an honest consequence of fixing the bigger error --
the corrected model is more parameter-sensitive, not less.

HELPER EDGE (helper_edge) -- NO LONGER AN ACTIVE DEFAULT, KEPT AS AN OPTION
------------------------------------------------------------------------------
A separate, coat-independent difference: an M13KO7 particle may simply infect
better than a pHK002.2 particle WEARING THE SAME COAT -- the natural virus versus
an engineered phagemid carrying a different genome. `helper_edge` is how many
times better (1 = no difference). It divides the infectivity of EVERY
pHK002.2 particle, wild-type or fully modified, so it would stack with the coat gap:

    pHK002.2, wild-type coat        1 / helper_edge
    pHK002.2, fully modified coat   eps / helper_edge
    M13KO7,   wild-type coat        1
    M13KO7,   hybrid coat           a_h  (see THE HYBRID COAT'S INFECTIVITY, DERIVED)

On review this never had a real mechanism behind it, and does not survive
scrutiny: infection (pIII binding the F-pilus tip, retraction, entry) is
mediated entirely by the coat, so if two particles genuinely carry an
identical coat, there is no obvious biophysical reason one should attach and
enter better than the other -- the packaged genome is not even read until
after entry. `helper_edge` therefore now DEFAULTS TO 1 (no effect) at every
level, including recursion.py, rather than being an active, unexplained edge
stacked on top of the coat gap. It remains available as a knob (e.g. to
explore a genuine post-entry establishment difference between a
self-sufficient phage genome and a complementation-dependent phagemid, which
is a real but DIFFERENT mechanism than "infects better" -- see
MODEL_REPORT.md section 7.10), but is no longer treated as an established
effect. Compare origin_penalty, below: a coat-independent, M13KO7-specific
effect that does have a stated mechanism (the KanR interruption), and is the
one now carrying that role in the model.

WHEN THE RATIO FLIPS. In the FIRST Phase 2 round the particles are still
wild-type, so the two lineages start at the Jumpstart's own split (adjusted by
helper_edge): the coat gap cannot act until Phase 2 has made particles with the
new coats. The made
particles then compete in round 2, and that is where the better-infecting
hybrid coat takes over. The wild-type carry-over also props up the phagemid's
plate count in round 1, because it is weighted at full infectivity. Every
advantage is therefore an explicit choice, not something hidden in the
structure.

WHAT IS LEFT OUT ON PURPOSE
------------------------------
DUAL-INFECTED CELLS. A cell already infected by one phage is very unlikely to
be infected by the other (superinfection exclusion), so the lineages do not mix
inside a cell. That removes everything dual infection needed in an earlier
version of this file: genome template pools, a packaging split between the two
genomes (alpha), a cap on copies per cell (Ncap), and a coat rule for phagemid
made in dual cells. That version is kept as phase2_dual_infection.py.

Also left out: an eclipse delay (an infected cell exports at full rate at
once), any dependence of export on growth state, and within-round randomness.

EXCLUSION (sigma)
--------------------
An infected cell (Ip or Ih) has lost most of its F-pili, so it adsorbs any
particle at only sigma times the rate of a naive cell (default 0.05,
unmeasured). A particle that lands on one is simply lost: its DNA goes in but
the cell's state does not change.

PRODUCTION
------------
An infected cell exports particles at a constant rate. Ih cells export at k2.
Ip cells export at k2 * (Lh/Lp)**size_effect: size_effect=1 is the hypothesis
that a shorter genome is packaged faster (rate ~ 1/length, 2.29x for the
phagemid); size_effect=0 (default) removes any production difference. The
phagemid MAY well have a replication or length advantage; this switch lets you
say so, but the effect of the coat is what the model is built to show.

INFECTION
----------
Mass action with a Holling type II ceiling, as in jumpstart.py: an infected
particle removes itself, and a cell can take up at most 1/h particles per
minute however much phage is around. Each particle type is weighted by its own
infectivity.

APPROXIMATIONS TO KNOW ABOUT
------------------------------
1. eps, m13ko7_copies, hp17_copies and sigma are not this construct's own
   measurements (m13ko7_copies and hp17_copies come from general M13/pSC101
   literature, not a pHK002.2/HP17KO7-specific assay -- see THE HYBRID
   COAT'S INFECTIVITY, DERIVED). The size of the fall is set by
   a_h/eps (about 19x at eps=0.05, the defaults' a_h~0.96); its existence
   needs only that the hybrid coat infects better.
2. Constant export from the moment of infection, no ramp-up.
3. Absolute titers matter for how long the wild-type carry-over lasts: it is
   consumed only as particles infect cells, so at the very high titers this
   model produces most of it survives the round.

STARTING CONDITIONS
----------------------
N(0)=B0 (fresh mid-log cells), Ip(0)=Ih(0)=0, and the four particle pools
Pw, Pt, Hw, Hh from whatever passage() delivered. From the Jumpstart only Pw
and Hw are nonzero.
"""
import numpy as np


def hybrid_infectivity(eps, m13ko7_copies=20.0, hp17_copies=8.0):
    """
    a_h: infectivity of the HYBRID coat (an Ih cell's M13KO7 particles),
    derived from coat-protein gene dosage rather than fit as an exponent.
    See module docstring, "THE HYBRID COAT'S INFECTIVITY, DERIVED".

        wt_frac = m13ko7_copies / (m13ko7_copies + hp17_copies)
        a_h     = wt_frac + (1 - wt_frac) * eps
    """
    wt_frac = m13ko7_copies / (m13ko7_copies + hp17_copies)
    return wt_frac + (1.0 - wt_frac) * eps


def derivatives(y, k1, h, k2p, k2h, k3, k4, K, eps, a_h, sigma, rho):
    """
    y = [N, Ip, Ih, Pw, Pt, Hw, Hh]; returns d(y)/dt.

    Infectivity weights (rho = 1/helper_edge divides every pHK002.2 particle):
        Pw (phagemid, wild-type coat)   rho
        Pt (phagemid, fully modified)   rho * eps
        Hw (helper, wild-type coat)     1
        Hh (helper, hybrid coat)        a_h  (see hybrid_infectivity())
    """
    N, Ip, Ih, Pw, Pt, Hw, Hh = np.maximum(y, 0.0)

    a_w = rho                      # phagemid, wild-type coat
    a_t = rho * eps                # phagemid, fully modified coat
    # a_h passed in directly -- computed once per run() call, not per RK4 substep

    # Shared saturation: every particle type competes for the same pili.
    sat = 1.0 / (1.0 + k1 * h * (a_w * Pw + a_t * Pt + Hw + a_h * Hh))

    infect_p = k1 * (a_w * Pw + a_t * Pt) * sat   # per-cell rate: a phagemid particle is caught
    infect_h = k1 * (Hw + a_h * Hh) * sat         # per-cell rate: an M13KO7 particle is caught
    N_to_Ip = infect_p * N
    N_to_Ih = infect_h * N

    room = 1.0 - (N + Ip + Ih) / K
    dN = -(N_to_Ip + N_to_Ih) + k3 * N * room
    dIp = N_to_Ip + k4 * Ip * room
    dIh = N_to_Ih + k4 * Ih * room

    # Each particle is lost when it lands on a naive cell (that is the infection)
    # or, at sigma times the rate, on an infected cell (its DNA goes in, nothing changes).
    sink = k1 * sat * (N + sigma * (Ip + Ih))
    dPw = -sink * a_w * Pw
    dPt = -sink * a_t * Pt + k2p * Ip
    dHw = -sink * Hw
    dHh = -sink * a_h * Hh + k2h * Ih
    return np.array([dN, dIp, dIh, dPw, dPt, dHw, dHh])


def plate_counts(Pw, Pt, Hw, Hh, eps, a_h, helper_edge=1.0):
    """
    What a plate counts: particles WEIGHTED BY INFECTIVITY, since a colony only
    forms if the particle infects an indicator cell.

        Amp colonies ~ (Pw + eps * Pt) / helper_edge   (pHK002.2: wild-type + fully modified)
        Kan colonies ~ Hw + a_h * Hh                    (M13KO7:   wild-type + hybrid)

    a_h is required, not defaulted -- see hybrid_infectivity() to derive it
    from eps and coat-protein gene dosage, rather than guessing a number here.
    Only the Amp:Kan RATIO is meaningful (volume plated, dilution and plating
    efficiency cancel). Returns (amp, kan).
    """
    return (Pw + eps * Pt) / helper_edge, Hw + a_h * Hh


def run(Pw_in, Hw_in, Pt_in=0.0, Hh_in=0.0, B0=3.0e8, T=900.0, steps=4000,
        k1=3.0e-11, h=1.0, Lp=3.8, Lh=8.7, k2=3.0, k3=np.log(2) / 25.0,
        K=1.0e9, growth_penalty=0.6, eps=1.0, m13ko7_copies=20.0, hp17_copies=8.0,
        sigma=0.05, size_effect=0.0, helper_edge=1.0, origin_penalty=1.0):
    """
    Integrate one overnight round (RK4, fixed step).

    Pw_in, Hw_in : wild-type-coat pHK002.2 and M13KO7 particles/mL at the start
                   (what the Jumpstart delivers, after passage()).
    Pt_in, Hh_in : fully modified pHK002.2 and hybrid M13KO7 particles/mL
                   carried in from a previous Phase 2 round (0 for round 1).
    eps          : infectivity of the fully modified coat relative to
                   wild-type. eps=1 = no difference.
    m13ko7_copies, hp17_copies : coat-protein gene dosage inside an Ih cell,
                   used to derive the hybrid coat's infectivity a_h -- see
                   hybrid_infectivity() and the module docstring, "THE
                   HYBRID COAT'S INFECTIVITY, DERIVED".
    sigma        : adsorption rate of an infected cell relative to a naive one.
    helper_edge  : how many times better an M13KO7 particle infects than a
                   pHK002.2 particle with the same coat (1 = no difference).
    size_effect  : exponent on (Lh/Lp) for the phagemid's export rate.
                   0 = same as helper (default); 1 = rate ~ 1/length.
    origin_penalty : multiplier on an Ih cell's own export rate (k2h),
                   reflecting M13KO7's KanR-interrupted origin slowing its
                   OWN plus-strand replication -- unlike jumpstart.py, Phase 2
                   cells are never dual-infected (see module docstring), so
                   this is not a competition effect the way it is in
                   jumpstart_mechanistic.py's replication_ratio; it is a
                   standalone assumption that the interrupted origin still
                   costs M13KO7 some export speed even racing alone.
                   1 = no penalty (Ih exports exactly as fast as Ip).
                   UNMEASURED; default 0.9 is a deliberately modest guess --
                   see MODEL_REPORT.md section 7.13.
    Other parameters match jumpstart.py's run().

    Returns (t, Y): Y has one row per timepoint, 7 columns
    [N, Ip, Ih, Pw, Pt, Hw, Hh].
    """
    k2h = k2 * origin_penalty
    k2p = k2 * (Lh / Lp) ** size_effect
    k4 = growth_penalty * k3
    a_h = hybrid_infectivity(eps, m13ko7_copies, hp17_copies)
    dt = T / steps
    Y = np.zeros((steps + 1, 7))
    Y[0] = [B0, 0.0, 0.0, Pw_in, Pt_in, Hw_in, Hh_in]
    args = (k1, h, k2p, k2h, k3, k4, K, eps, a_h, sigma, 1.0 / helper_edge)
    for i in range(steps):
        y = Y[i]
        s1 = derivatives(y, *args)
        s2 = derivatives(y + 0.5 * dt * s1, *args)
        s3 = derivatives(y + 0.5 * dt * s2, *args)
        s4 = derivatives(y + dt * s3, *args)
        Y[i + 1] = np.maximum(y + (dt / 6.0) * (s1 + 2 * s2 + 2 * s3 + s4), 0.0)
    return np.arange(steps + 1) * dt, Y


# ==============================================================================
# COMMAND-LINE INTERFACE
#   python3 phase2.py                 # neutral defaults: no difference between lineages
#   python3 phase2.py -e 0.05         # hybrid coat infects better than fully modified
#   python3 phase2.py -e 0.05 --plot out.png --csv trace.csv
# ==============================================================================
import argparse
import csv as csvlib


def build_parser():
    p = argparse.ArgumentParser(
        description="Phase 2: one overnight round, two lineages (pHK002.2 vs M13KO7) in HP17KO7 cells.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-P", "--Pw-in", type=float, default=1.946e11, dest="Pw_in",
                   help="wild-type-coat pHK002.2 particles/mL entering from the Jumpstart "
                        "(default: recursion.py's Jumpstart output / 11; test_model.py "
                        "checks it has not drifted)")
    p.add_argument("-H", "--Hw-in", type=float, default=3.987e10, dest="Hw_in",
                   help="wild-type-coat M13KO7 particles/mL entering from the Jumpstart "
                        "(default: recursion.py's Jumpstart output / 11; test_model.py "
                        "checks it has not drifted)")
    p.add_argument("--Pt-in", type=float, default=0.0, dest="Pt_in",
                   help="fully modified pHK002.2 particles/mL carried in from an earlier Phase 2 round")
    p.add_argument("--Hh-in", type=float, default=0.0, dest="Hh_in",
                   help="hybrid M13KO7 particles/mL carried in from an earlier Phase 2 round")
    p.add_argument("-e", "--eps", type=float, default=1.0,
                   help="infectivity of the fully modified coat relative to wild-type "
                        "(1 = no difference between the lineages)")
    p.add_argument("--m13ko7-copies", type=float, default=20.0, dest="m13ko7_copies",
                   help="M13KO7 genome copies in an infected cell. Default 20: its own "
                        "p15A copy number, NOT wild-type M13's ~130-300 f1 rolling-circle "
                        "number (M13KO7's f1 origin is the interrupted one). Published "
                        "p15A range ~10-30; 20 is the end of it the measured round-1 data "
                        "favours. Unmeasured for this construct and the answer is "
                        "sensitive to it (1.73x across 10-30) -- see the docstring")
    p.add_argument("--hp17-copies", type=float, default=8.0, dest="hp17_copies",
                   help="HP17's own copy number (Rep101/pSC101 origin, typical published "
                        "range ~5-12; not this construct's own measurement)")
    p.add_argument("--helper-edge", type=float, default=1.0, dest="helper_edge",
                   help="how many times better an M13KO7 particle infects than a pHK002.2 "
                        "particle with the SAME coat (1 = no difference)")
    p.add_argument("--sigma", type=float, default=0.05,
                   help="adsorption rate of an infected cell relative to a naive cell")
    p.add_argument("--size-effect", type=float, default=0.0, dest="size_effect",
                   help="exponent on (Lh/Lp) for the phagemid's export rate (0 = none, 1 = ~1/length)")
    p.add_argument("--origin-penalty", type=float, default=1.0, dest="origin_penalty",
                   help="multiplier on an Ih cell's export rate, for M13KO7's interrupted "
                        "origin slowing its own replication even without competition "
                        "(1 = no penalty). UNMEASURED")
    p.add_argument("-L", "--Lp", type=float, default=3.8, help="phagemid genome length, kb")
    p.add_argument("--Lh", type=float, default=8.7, help="M13KO7 genome length, kb")
    p.add_argument("-k", "--k2", type=float, default=3.0, help="export rate of an infected cell, particles/min")
    p.add_argument("--k1", type=float, default=3.0e-11, help="adsorption rate constant, mL/min")
    p.add_argument("--h", type=float, default=1.0, help="pilus handling time, min")
    p.add_argument("--k3", type=float, default=np.log(2) / 25.0, help="growth rate of naive cells, /min")
    p.add_argument("--growth-penalty", type=float, default=0.6, dest="growth_penalty",
                   help="growth rate of infected cells as a fraction of naive")
    p.add_argument("--K", type=float, default=1.0e9, help="carrying capacity, cells/mL")
    p.add_argument("--B0", type=float, default=3.0e8, help="starting naive cells/mL")
    p.add_argument("-T", "--time", type=float, default=900.0, help="round length, min")
    p.add_argument("--steps", type=int, default=4000, help="integration steps")
    p.add_argument("--csv", metavar="FILE", default=None, help="save the full time-course")
    p.add_argument("--plot", metavar="FILE", default=None, help="save a figure")
    p.add_argument("-q", "--quiet", action="store_true", help="skip the printed summary")
    return p


def main():
    a = build_parser().parse_args()
    kw = dict(B0=a.B0, T=a.time, k1=a.k1, h=a.h, Lp=a.Lp, Lh=a.Lh, k2=a.k2, k3=a.k3,
              K=a.K, growth_penalty=a.growth_penalty, eps=a.eps,
              m13ko7_copies=a.m13ko7_copies, hp17_copies=a.hp17_copies,
              sigma=a.sigma, size_effect=a.size_effect, helper_edge=a.helper_edge,
              origin_penalty=a.origin_penalty)
    t, Y = run(a.Pw_in, a.Hw_in, Pt_in=a.Pt_in, Hh_in=a.Hh_in, steps=a.steps, **kw)
    N, Ip, Ih, Pw, Pt, Hw, Hh = Y[-1]
    Pp, Ph = Pw + Pt, Hw + Hh
    a_h = hybrid_infectivity(a.eps, a.m13ko7_copies, a.hp17_copies)
    amp, kan = plate_counts(Pw, Pt, Hw, Hh, a.eps, a_h, a.helper_edge)
    if not a.quiet:
        ratio_in = (a.Pw_in + a.Pt_in) / (a.Hw_in + a.Hh_in)
        wt_frac = a.m13ko7_copies / (a.m13ko7_copies + a.hp17_copies)
        print("=" * 60 + "\nPHASE 2, TWO LINEAGES\n" + "=" * 60)
        print(f"  in : wild-type Pw={a.Pw_in:.3e}  Hw={a.Hw_in:.3e}   (Amp:Kan particles {ratio_in:.3g})")
        print(f"  fully modified coat eps = {a.eps:g};  hybrid a_h = {a_h:.3g}  "
              f"(wt_frac={wt_frac:.3g} from m13ko7_copies={a.m13ko7_copies:g}, hp17_copies={a.hp17_copies:g})"
              f"  -> helper out-infects phagemid by {a_h / a.eps:.3g}x")
        print(f"  helper edge = {a.helper_edge:g}x (M13KO7 vs pHK002.2 with the same coat)")
        print(f"  origin penalty = {a.origin_penalty:g}x on Ih export (UNMEASURED, interrupted-origin cost)")
        print(f"  sigma={a.sigma:g}, size_effect={a.size_effect:g}")
        print("-" * 60)
        print(f"  final N  = {N:.3e}   Ip = {Ip:.3e}   Ih = {Ih:.3e}")
        print(f"  final Pw = {Pw:.3e}   Pt = {Pt:.3e}   (pHK002.2: carried-over wild-type, fully modified)")
        print(f"  final Hw = {Hw:.3e}   Hh = {Hh:.3e}   (M13KO7:   carried-over wild-type, hybrid)")
        print(f"\n  Amp:Kan (particles) = {Pp / Ph if Ph > 0 else float('inf'):.4g}   <- what is physically in the tube")
        print(f"  Amp:Kan (COLONIES)  = {amp / kan if kan > 0 else float('inf'):.4g}   <- what a PLATE would show")
        print(f"  cell count check: total = {N + Ip + Ih:.4e} (K={a.K:.2e})")
        print("\n  convergence (colony ratio vs integration steps):")
        for s in (1000, 4000, 16000):
            _, Yc = run(a.Pw_in, a.Hw_in, Pt_in=a.Pt_in, Hh_in=a.Hh_in, steps=s, **kw)
            am, ka = plate_counts(Yc[-1, 3], Yc[-1, 4], Yc[-1, 5], Yc[-1, 6], a.eps, a_h, a.helper_edge)
            print(f"    steps={s:6d}  {am / ka if ka > 0 else float('inf'):.6g}")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csvlib.writer(f)
            w.writerow(["t_min", "N", "Ip", "Ih", "Pw", "Pt", "Hw", "Hh"])
            for tt, row in zip(t, Y):
                w.writerow([tt] + list(row))
        print(f"\n  wrote {a.csv}")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            import plotstyle as ps          # shared palette; figures still work without it
            ps.apply_style()
            c_N, c_P, c_H = ps.NAIVE, ps.PHAGEMID, ps.HELPER
        except ImportError:
            ps = None
            c_N, c_P, c_H = "#8a8580", "#1f6fb4", "#d1495b"
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 8.6))
        for col, lab, c in ((0, "N (naive)", c_N), (1, "Ip (pHK002.2)", c_P), (2, "Ih (M13KO7)", c_H)):
            m = Y[:, col] > 0
            ax1.plot(t[m], Y[m, col], color=c, lw=2, label=lab)
        ax1.set_yscale("log"); ax1.set_ylim(1, 3 * a.K); ax1.set_xlabel("time (min)"); ax1.set_ylabel("cells / mL")
        ax1.set_title("Cell populations", loc="left", weight="bold"); ax1.legend(frameon=False)
        ax2.plot(t, Y[:, 4], color=c_P, lw=2.2, label="Pt (pHK002.2, fully modified)")
        ax2.plot(t, Y[:, 6], color=c_H, lw=2.2, label="Hh (M13KO7, hybrid)")
        ax2.plot(t, Y[:, 3], "--", color=c_P, lw=1.6, alpha=0.7, label="Pw (carried-over wild-type)")
        ax2.plot(t, Y[:, 5], "--", color=c_H, lw=1.6, alpha=0.7, label="Hw (carried-over wild-type)")
        ax2.set_yscale("log"); ax2.set_xlabel("time (min)"); ax2.set_ylabel("particles / mL")
        ax2.set_title("Free particles (raw count)", loc="left", weight="bold"); ax2.legend(frameon=False)

        # Pw/Hw look flat in ax2 above -- they are NOT: they decline by roughly
        # 10-15% over the round (see phase2.py's WHAT IS LEFT OUT / APPROXIMATIONS
        # section, and MODEL_REPORT.md section 7.3), it is just invisible against
        # a log axis spanning 5+ orders of magnitude. This panel is the SAME two
        # curves, normalized to their own starting value on a LINEAR axis, so the
        # real decline is visible instead of hidden by the other panel's scale.
        Pw0, Hw0 = Y[0, 3], Y[0, 5]
        if Pw0 > 0:
            ax3.plot(t, 100 * Y[:, 3] / Pw0, color=c_P, lw=2.2, label="Pw (% of round-start Pw)")
        if Hw0 > 0:
            ax3.plot(t, 100 * Y[:, 5] / Hw0, color=c_H, lw=2.2, label="Hw (% of round-start Hw)")
        if Pw0 > 0:
            ax3.annotate(f"{100 * Y[-1, 3] / Pw0:.1f}%", (t[-1], 100 * Y[-1, 3] / Pw0),
                        xytext=(6, 4), textcoords="offset points", color=c_P, fontsize=9.5, weight="bold")
        if Hw0 > 0:
            ax3.annotate(f"{100 * Y[-1, 5] / Hw0:.1f}%", (t[-1], 100 * Y[-1, 5] / Hw0),
                        xytext=(6, -12), textcoords="offset points", color=c_H, fontsize=9.5, weight="bold")
        ax3.set_ylim(0, 105); ax3.set_xlabel("time (min)"); ax3.set_ylabel("% of round-start value")
        ax3.set_title("Carry-over survival (linear)", loc="left", weight="bold"); ax3.legend(frameon=False, loc="lower left")

        # RAW particle count (ax2) can be misleading on its own: Pt piles up despite
        # being a weak infector, because production (k2p*Ip) does not depend on eps
        # at all, and a LOW a_t also means Pt is consumed slowly once made. This
        # panel is what a PLATE would actually show -- the same trajectory, weighted
        # by infectivity -- which is why it can invert relative to ax2.
        amp_t, kan_t = plate_counts(Y[:, 3], Y[:, 4], Y[:, 5], Y[:, 6], a.eps, a_h, a.helper_edge)
        ax4.plot(t, amp_t, color=c_P, lw=2.2, label="Amp (pHK002.2, infectivity-weighted)")
        ax4.plot(t, kan_t, color=c_H, lw=2.2, label="Kan (M13KO7, infectivity-weighted)")
        ax4.set_yscale("log"); ax4.set_xlabel("time (min)"); ax4.set_ylabel("plate-equivalent count / mL")
        ax4.set_title("What a plate would show", loc="left", weight="bold"); ax4.legend(frameon=False)

        for ax in (ax1, ax2, ax3, ax4):
            if ps:
                ps.finish_axes(ax, legend=False)      # legends are drawn above
            else:
                ax.grid(alpha=0.2)
                for s_ in ("top", "right"):
                    ax.spines[s_].set_visible(False)
        fig.suptitle(f"eps={a.eps:g} (hybrid a_h={a_h:.2g}), sigma={a.sigma:g}, "
                     f"size_effect={a.size_effect:g}", fontsize=10, color="#555")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200, facecolor="white")
        print(f"\n  wrote {a.plot}")


if __name__ == "__main__":
    main()
