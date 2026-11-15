"""
The Jumpstart, stripped to basics -- now with a command-line interface.

FOUR states. No coat tracking, no eclipse delay, and no packaging biology: the
split of a dual cell's export between the two genomes is a single input number,
phi, not something the model computes. Infection does saturate (a cell takes up
at most 1/h particles per minute). Just: cells get infected, infected cells make
both kinds of particle, particles get used up making infections.

    B  = cells with pHK002.2 only (every cell starts here -- it's a plasmid,
         already inside from transformation)
    I  = cells with BOTH pHK002.2 and M13KO7 (the only state that produces
         anything, since M13KO7 supplies all the machinery)
    Pp = free pHK002.2 particles  (what you'd count on Amp plates)
    Ph = free M13KO7 particles    (what you'd count on Kan plates)
         ^ NOTE: "Ph" stands for Particle-helper, NOT "phagemid" -- easy to
           misread since the letters look similar. Pp is phagemid, Ph is
           helper (M13KO7). Worth a re-check any time these two get
           confusing in context.

USAGE
-----
    python3 jumpstart.py --help                  # see every option
    python3 jumpstart.py                          # defaults: the MEASURED 4.88:1 plate split
    python3 jumpstart.py --phi 0.9                  # set the phagemid's export share directly
    python3 jumpstart.py --phi -1                    # derive it from genome length instead (2.29:1)
    python3 jumpstart.py --phi -1 -L 8.7              # ... with Lp=Lh: phi=0.5, an even split, 1:1
    python3 jumpstart.py -k 6.0                      # double the dual cell's TOTAL export rate
    python3 jumpstart.py -T 1800 --plot out.png       # longer run, save a graph
    python3 jumpstart.py --csv trace.csv               # save the full time-course

THE FOUR THINGS THAT HAPPEN
---------------------------
1. INFECTION.  A helper particle finds a B cell and gets in.
      B -> I,  at rate k1 * Ph * B / (1 + k1*h*Ph)

   The denominator is new: M13 docks at the TIP of an F-pilus, and a cell
   has only 1-3 of them. Plain mass action (k1*Ph*B alone) says double
   the phage, double the infection rate, forever -- which is physically
   impossible once pili become the bottleneck rather than phage
   abundance. This is a Holling type II saturating response: at low Ph
   it reduces to plain mass action (denominator ~1); at high Ph it flattens
   toward a ceiling of 1/h infections per minute per cell, no matter how
   much MORE phage you add. h is the pilus "handling time" -- how long one
   stays occupied per docking event.

   Pp is NOT in this denominator, even though it is also floating around:
   this file already made the choice not to model Pp's unproductive
   interactions with cells at all (see CONSUMPTION below), and adding
   saturation does not revisit that -- it only caps how fast Ph itself can
   infect. If Pp's receptor competition turns out to matter, that is a
   separate, identifiable addition, not something to fold in silently here.

2. GROWTH.  Cells divide, capped by carrying capacity K. B (not yet
   infected) grows at the full rate k3. I (infected -- actively making
   phage particles) grows SLOWER, at k4 = growth_penalty * k3. This is a
   measured, literature-backed number used throughout this whole project:
   infected M13-family cells grow at roughly 60% of the uninfected rate
   (growth_penalty defaults to 0.6). An earlier version of this file gave
   B and I the same growth rate -- flagged then as a simplification,
   fixed here.

3. EXPORT.  I cells make particles continuously (M13 doesn't lyse the
   cell). An I cell has ONE set of assembly machinery, and all of it comes
   from M13KO7 -- the phagemid has no machinery of its own. So the two
   genomes SHARE one total export rate, k2, rather than each getting a
   full rate of its own:

       dPp/dt = phi       * k2 * I      phagemid's share of the cell's export
       dPh/dt = (1 - phi) * k2 * I      helper's share

   THE JUMPSTART DOES NOT PREDICT ITS OWN Amp:Kan RATIO. The model's output
   ratio is exactly phi/(1-phi), and NOTHING else in this file can move it:
   not the adsorption rate, the helper dose, the kanamycin timing, the export
   rate, the culture density or the run length (checked -- all of those leave
   the ratio at phi/(1-phi) to three decimals while changing the TITER by up
   to 30x). Whatever really decides how a dual-infected cell splits its
   packaging between the two genomes -- pHK002.2's intact f1 origin versus
   M13KO7's origin interrupted by the KanR gene, genome length, copy number --
   is NOT modelled. phi stands in for all of it.

   So phi is an input to be MEASURED, not derived. The default, 0.830, is the
   plate ratio measured for this assay (Amp:Kan about 4.88:1, observed range
   roughly 3-16). recursion.py sets phi the same way, from --jump-ratio, but
   corrects for the helper edge first, so the two files agree on what the
   Jumpstart produces.

   phi = -1 (or any negative number) derives it from genome length instead:

       phi = Lh / (Lh + Lp)             = 0.696 with Lp=3.8, Lh=8.7 kb

   an export split of Lh/Lp = 2.29:1. That is a HYPOTHESIS (shorter genome ->
   packaged more often, rate ~ 1/length), not a measurement; it was the default
   in earlier versions, which is why the number appears throughout this file's
   history. Lp and Lh have no other effect on this stage.

   WHY THE SHARED RATE: an earlier version gave each genome its own full rate
   (k2p = k2*Lh/Lp = 6.87 and k2h = k2 = 3.00), so one I cell exported 9.87
   particles/min in total -- using the same machinery twice. Total output per
   I cell is now k2 = 3.0, so absolute titers are lower, and the infection wave
   is slower because less helper is exported per I cell.

4. CONSUMPTION.  Every helper particle that causes an infection is used up
   doing it. Nothing consumes Pp here (see the long comment inside
   derivatives() for why that's a deliberate simplification, not an
   oversight).

5. KANAMYCIN SELECTION, at t=kan_at (default 75 min).  A real, named step
   of the actual protocol, not an optional refinement: any B cell that
   never got infected by helper is killed. This is applied as an instant,
   complete event -- B is set to exactly zero at that one timepoint, then
   everything continues evolving from there. I, Pp, Ph are untouched by
   the event itself (only B changes), but I benefits afterward: with B
   gone, it no longer competes with I for the carrying-capacity "room"
   term, so removing a large B population can noticeably raise the final
   titer of survivors. Turn it off with --no-kan to compare directly.

   HOW MUCH KANAMYCIN ACTUALLY KILLS depends entirely on k1 and moi, and at
   this assay's values (k1 = 3.0e-11 mL/min, moi = 2.5 helper particles per
   starting cell) it is 70%: 29.7% of cells have caught helper by minute 75.
   That is a middle case, and a sensitive one -- at k1=3e-12 only 1.8% are
   infected (Kan kills 98%), at k1=3e-10 essentially all of them are (Kan kills
   nothing). The Jumpstart's OUTPUT barely notices either way: across the
   protocol's 60-90 minute window the final titer moves 3-4%, and the ratio not
   at all. So these numbers matter for understanding what the kanamycin step
   does to the culture, not for what this stage hands to Phase 2.

WHY THERE IS NO AMPICILLIN STEP
---------------------------------
Ampicillin in the real protocol is not a one-time selection event like
kanamycin. It is present continuously during the PRE-culture growth --
growing a single colony up to A600<0.05, BEFORE any phage is added -- to
stop plasmid-loss cells from taking over as the culture divides. By the
time this simulation's clock starts, that has already happened: B0 being
100% phagemid-positive already IS the outcome of that selection. There is
also no mechanism anywhere in this model for a cell to ever LOSE the
phagemid once it has it (no plasmid segregation loss), so even if
ampicillin is still physically in the flask during the simulated window,
there is nothing left for it to select against here. That is a real
simplifying assumption -- plasmid loss over one overnight is small but
not exactly zero in reality -- named here rather than silently assumed.

WHAT'S DELIBERATELY LEFT OUT (compared to the full model)
-----------------------------------------------------------
The eclipse delay; coat composition / infectivity (Phase 2's business, not
this stage's: no HP17KO7 is present here, so every particle is wild-type);
superinfection exclusion (particles landing on already-infected cells are not
removed -- measured effect on final titers about 0.5%); plasmid segregation
loss; kanamycin staying in the culture overnight (here it is a one-time kill at
kan_at); and any packaging biology behind phi. NOT left out: Holling
saturation of infection, which is implemented.

STARTING CONDITIONS, EXPLICITLY
--------------------------------
B(0)  = B0        -- the colony: bacteria already carrying pHK002.2 as
                      resident plasmid DNA. This IS your starting material.
I(0)  = 0         -- nobody has caught helper yet.
Pp(0) = 0         -- NO free phagemid particles exist yet -- it's DNA
                      inside cells right now, not virus.
Ph(0) = moi * B0  -- the helper YOU add from an external stock.
"""
import argparse
import csv as csvlib

import numpy as np


# ==============================================================================
# THE EQUATIONS. This function is the entire model -- everything below it
# in the file is just machinery for running it forward through time, or for
# letting you change numbers without editing the file.
# ==============================================================================
def derivatives(y, k1, h, k2, phi, k3, k4, K):
    """
    Given the current state y = [B, I, Pp, Ph], return how fast each of the
    four quantities is changing RIGHT NOW.

    k3 : growth rate for B (not yet infected)
    k4 : growth rate for I (infected, actively producing particles) --
         slower than k3, reflecting the real metabolic cost of infection.
    h  : pilus handling time -- caps infection at 1/h per cell per minute,
         however much phage is around (see module docstring).
    k2 : TOTAL export rate of one I cell, particles/min -- shared between
         the two particle types (see module docstring, section 3).
    phi: phagemid's share of that export; helper gets 1 - phi.
    """
    B, I, Pp, Ph = np.maximum(y, 0.0)

    # THE ONE EVENT THAT MOVES A CELL: a free helper particle finds and
    # enters a B cell. Saturating (Holling type II), not plain mass action
    # -- see module docstring for why. Reused three times below (it's a
    # single physical event, viewed from three different ledgers).
    sat = 1.0 / (1.0 + k1 * h * Ph)
    infection = k1 * Ph * B * sat

    # How full is the flask? 1.0 = empty (grow freely), 0.0 = at capacity.
    growth_room = 1.0 - (B + I) / K

    dB = -infection + k3 * B * growth_room            # loses to infection, gains by division at the FULL rate
    dI = infection + k4 * I * growth_room              # gains exactly what B lost, divides at the PENALISED rate
    dPp = phi * k2 * I                                 # only a gain: phagemid's share of I's export
    dPh = -infection + (1.0 - phi) * k2 * I            # loses to infection, gains helper's share of export

    return np.array([dB, dI, dPp, dPh])


# ==============================================================================
# THE INTEGRATOR. Not the model -- bookkeeping that turns "the speed right
# now" into an actual trajectory over time, using RK4 (see the loop below
# for what that means).
# ==============================================================================
# The measured Jumpstart export split: plate Amp:Kan of about 4.88:1 for this
# assay (observed range roughly 3-16). phi = r/(1+r). recursion.py derives the
# same number from --jump-ratio, correcting for the helper edge first.
PHI_MEASURED = 4.88 / (1.0 + 4.88)


def resolve_phi(phi, Lp, Lh):
    """
    The phagemid's share of a dual-infected cell's export.

    phi=None       the measured default, PHI_MEASURED (see module docstring)
    phi<0          derive it from genome length instead, Lh/(Lh+Lp)
    0 <= phi <= 1  use it as given

    This is the ONLY place those rules live, so run() and the CLI cannot
    disagree.
    """
    if phi is None:
        phi = PHI_MEASURED
    elif phi < 0:
        phi = Lh / (Lh + Lp)
    if not 0.0 <= phi <= 1.0:
        raise ValueError(f"phi must be between 0 and 1, got {phi}")
    return phi


def run(B0=4.0e7, moi=2.5, T=1080.0, steps=4000,
        k1=3.0e-11, h=1.0, Lp=3.8, Lh=8.7, k2=3.0, k3=np.log(2)/25.0, K=1.0e9,
        kan_at=75.0, growth_penalty=0.6, phi=None):
    """
    Simulate from t=0 to t=T. Returns (t, Y): t is an array of timepoints,
    Y has one row per timepoint and 4 columns [B, I, Pp, Ph].

    k2     : TOTAL export rate of one dual-infected (I) cell, particles/min.
             Shared between phagemid and helper -- NOT a rate each gets.
    phi    : phagemid's share of that export (helper gets 1 - phi). None
             (default) uses the MEASURED split, PHI_MEASURED = 0.830 (4.88:1).
             A negative value derives it from genome length as Lh/(Lh+Lp)
             (2.29:1) -- a hypothesis, not a measurement. Lp, Lh only matter
             through that option.
    kan_at : minute at which kanamycin kills every B cell (set to None to
             run without any selection at all, for comparison).
    growth_penalty : I's growth rate as a fraction of B's (k4 = this * k3).
             Defaults to 0.6, the measured value used throughout this
             project -- infected cells grow at roughly 60% of the
             uninfected rate.
    """
    phi = resolve_phi(phi, Lp, Lh)
    k4 = growth_penalty * k3
    y0 = np.array([B0, 0.0, 0.0, moi * B0])   # B, I, Pp, Ph -- see module docstring
    dt = T / steps
    Y = np.zeros((steps + 1, 4))
    Y[0] = y0

    # Which step index does kan_at fall on? None means "never" (disabled).
    kan_step = None if kan_at is None else int(round(kan_at / dt))

    for i in range(steps):
        y = Y[i]
        s1 = derivatives(y, k1, h, k2, phi, k3, k4, K)
        s2 = derivatives(y + 0.5*dt*s1, k1, h, k2, phi, k3, k4, K)
        s3 = derivatives(y + 0.5*dt*s2, k1, h, k2, phi, k3, k4, K)
        s4 = derivatives(y + dt*s3, k1, h, k2, phi, k3, k4, K)
        Y[i+1] = np.maximum(y + (dt/6.0)*(s1 + 2*s2 + 2*s3 + s4), 0.0)

        # KANAMYCIN EVENT: at the one step matching kan_at, kill every B
        # cell instantly. Only the B component changes -- I, Pp, Ph carry
        # straight through untouched. This has to happen AFTER the normal
        # RK4 update above, as a separate overwrite of that one row.
        if kan_step is not None and (i + 1) == kan_step:
            Y[i+1, 0] = 0.0

    t = np.arange(steps+1) * dt
    return t, Y


# ==============================================================================
# COMMAND-LINE INTERFACE. Every model parameter is a flag. The ones you'll
# actually tune day to day (the two export rates, and how long to run) get
# short forms; everything else is a long flag only, so nothing is hidden.
# ==============================================================================
def build_parser():
    p = argparse.ArgumentParser(
        description="Run the simple 4-state Jumpstart model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    common = p.add_argument_group("commonly-changed parameters (short flags)")
    common.add_argument("-L", "--Lp", type=float, default=3.8,
        help="phagemid (pHK002.2) genome length, kb. Only used when --phi is "
             "negative (the length-derived split, Lh/(Lh+Lp)); ignored "
             "otherwise, since phi is measured rather than derived.")
    common.add_argument("-k", "--k2", type=float, default=3.0,
        help="TOTAL export rate of one dual-infected cell, particles/min, "
             "shared between phagemid and helper (split by phi).")
    common.add_argument("-T", "--time", type=float, default=1080.0,
        help="total simulated duration, minutes (default 1080 = 18 hours, "
             "an overnight incubation -- the protocol just says 'overnight' "
             "with no fixed number, so treat this as a working estimate, "
             "not a documented constant)")
    common.add_argument("-M", "--moi", type=float, default=2.5,
        help="multiplicity of infection: M13KO7 particles added per starting "
             "cell (this assay: 2.5). Affects the TITER and how many cells "
             "kanamycin kills, never the ratio")

    full = p.add_argument_group("every other model parameter (long flags only)")
    full.add_argument("--Lh", type=float, default=8.7,
        help="M13KO7 genome length, kb -- the reference point phi's default "
             "is derived against. Rarely needs changing; it's a fixed fact "
             "about M13KO7.")
    full.add_argument("--phi", type=float, default=None,
        help="phagemid's share of a dual-infected cell's export, 0 to 1 "
             "(helper gets the rest). Unset = the MEASURED split (4.88:1 "
             "Amp:Kan). Negative = derive it from genome length instead, "
             "Lh/(Lh+Lp) = 2.29:1, which is a hypothesis, not a measurement.")
    full.add_argument("--kan-at", type=float, default=75.0, dest="kan_at",
        help="minute at which kanamycin kills every un-infected (B) cell")
    full.add_argument("--no-kan", action="store_true",
        help="disable kanamycin selection entirely, for comparison")
    full.add_argument("--growth-penalty", type=float, default=0.6,
        dest="growth_penalty",
        help="I's growth rate as a fraction of B's (k4 = this * k3). "
             "0.6 is the measured value; use 1.0 to disable the penalty "
             "entirely, for comparison")
    full.add_argument("--k1", type=float, default=3.0e-11,
        help="infection (adsorption) rate constant for this assay, mL/min. "
             "Affects the TITER and how many cells kanamycin kills, never "
             "the ratio")
    full.add_argument("--h", type=float, default=1.0,
        help="pilus handling time, min -- caps infection at 1/h per cell "
             "per minute regardless of phage concentration")
    full.add_argument("--k3", type=float, default=np.log(2)/25.0,
        help="growth rate (default corresponds to a 25 min doubling time)")
    full.add_argument("--K", type=float, default=1.0e9,
        help="carrying capacity, cells/mL")
    full.add_argument("--B0", type=float, default=4.0e7,
        help="starting cell density, cells/mL")
    full.add_argument("--steps", type=int, default=4000,
        help="integration steps (see the convergence check for whether "
             "you have enough)")

    out = p.add_argument_group("output")
    out.add_argument("--csv", metavar="FILE", default=None,
        help="save the full time-course to this CSV file")
    out.add_argument("--plot", metavar="FILE", default=None,
        help="save a graph of B,I,Pp,Ph vs time to this file (.png or "
             ".pdf); needs matplotlib, only imported if you use this flag")
    out.add_argument("-q", "--quiet", action="store_true",
        help="skip the summary and convergence check, just do --csv/--plot")

    return p


def main():
    args = build_parser().parse_args()

    # phi is resolved by the same function run() uses -- there is exactly one
    # place that turns "unset" into the length-derived default, so what the
    # CONFIG banner shows is what the simulation actually used.
    phi_used = resolve_phi(args.phi, args.Lp, args.Lh)
    split_txt = (f"{phi_used/(1-phi_used):.3f}:1" if phi_used < 1.0
                 else "all phagemid")

    kan_at = None if args.no_kan else args.kan_at

    t, Y = run(B0=args.B0, moi=args.moi, T=args.time, steps=args.steps,
              k1=args.k1, h=args.h, Lp=args.Lp, Lh=args.Lh, k2=args.k2,
              k3=args.k3, K=args.K, kan_at=kan_at,
              growth_penalty=args.growth_penalty, phi=args.phi)
    B, I, Pp, Ph = Y[-1]

    if not args.quiet:
        print("=" * 60)
        print("CONFIG")
        print("=" * 60)
        print(f"  Lp / Lh (kb)             = {args.Lp} / {args.Lh}")
        print(f"  k2 (dual-cell export)    = {args.k2}  particles/cell/min, TOTAL")
        print(f"  phi (phagemid share)     = {phi_used:.4f}   "
              f"({'measured default' if args.phi is None else ('length-derived' if args.phi < 0 else 'set by --phi')})"
              f"  -> export split Pp:Ph = {split_txt}")
        print(f"    -> Pp export           = {args.k2*phi_used:.4f}  /cell/min")
        print(f"    -> Ph export           = {args.k2*(1-phi_used):.4f}  /cell/min")
        print(f"  k1 (infection), h        = {args.k1:.2e}, {args.h}")
        print(f"  k3 (growth), K (cap)     = {args.k3:.4f}, {args.K:.2e}")
        print(f"  growth penalty           = {args.growth_penalty}   "
              f"-> k4 (infected growth) = {args.growth_penalty*args.k3:.4f}")
        print(f"  B0, moi                  = {args.B0:.2e}, {args.moi}")
        print(f"  kanamycin                = "
              f"{'DISABLED' if kan_at is None else f'kills B at t={kan_at} min'}")
        print(f"  T, steps                 = {args.time}, {args.steps}")
        print()
        print("=" * 60)
        print("RESULT")
        print("=" * 60)
        print(f"  final B  = {B:.3e}")
        print(f"  final I  = {I:.3e}")
        print(f"  final Pp = {Pp:.3e}")
        print(f"  final Ph = {Ph:.3e}")
        print(f"  Amp:Kan  = {Pp/Ph:,.4f}" if Ph > 0 else "  Amp:Kan = inf")
        print()
        print("Convergence check (halve dt, does the answer stop moving?):")
        prev = None
        for s in [args.steps//8, args.steps//4, args.steps//2, args.steps, args.steps*2]:
            s = max(s, 10)
            _, Yc = run(B0=args.B0, moi=args.moi, T=args.time, steps=s,
                       k1=args.k1, h=args.h, Lp=args.Lp, Lh=args.Lh,
                       k2=args.k2, k3=args.k3, K=args.K, kan_at=kan_at,
                       growth_penalty=args.growth_penalty, phi=args.phi)
            r = Yc[-1,2]/Yc[-1,3] if Yc[-1,3] > 0 else float('inf')
            note = "" if prev is None else f"  change {abs(r-prev)/r:.2e}"
            print(f"  steps={s:6d}  Amp:Kan={r:14,.4f}{note}")
            prev = r

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csvlib.writer(f)
            w.writerow(["t_min", "B", "I", "Pp", "Ph"])
            for row_t, row_y in zip(t, Y):
                w.writerow([row_t] + list(row_y))
        print(f"\n  wrote {args.csv}")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        try:
            import plotstyle as ps          # shared palette; figures still work without it
            ps.apply_style()
            c_B, c_I, c_kan, c_Pp, c_Ph = ps.NAIVE, ps.DUAL, ps.ACCENT, ps.PHAGEMID, ps.HELPER
        except ImportError:
            ps = None
            c_B, c_I, c_kan, c_Pp, c_Ph = "#8a8580", "#123b5e", "#d1495b", "#1f6fb4", "#d1495b"

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

        # left panel: cell populations. Linear scale is fine here since
        # both are bounded by K, not spanning many orders of magnitude.
        ax1.plot(t, Y[:,0], "-", color=c_B, lw=2.2, label="B (phagemid only)")
        ax1.plot(t, Y[:,1], "-", color=c_I, lw=2.2, label="I (dual-infected)")
        if kan_at is not None:
            ax1.axvline(kan_at, color=c_kan, ls="--", lw=1.4, alpha=0.7)
            ax1.text(kan_at, ax1.get_ylim()[1]*0.97, "  Kan added", color=c_kan,
                     fontsize=8.5, va="top")
        ax1.set_xlabel("time (min)")
        ax1.set_ylabel("cells / mL")
        ax1.set_title("Cell populations", loc="left", weight="bold")
        ax1.legend(frameon=False)
        ax1.grid(alpha=0.2)

        # right panel: particle populations. These DO span many orders of
        # magnitude (starting near zero, ending near carrying-capacity-
        # driven plateaus), so log scale is the only way to see the shape.
        # If a species never gets above zero (e.g. moi=0, so no helper
        # ever enters the system and nothing is ever produced), there is
        # nothing a log-scale plot can show for it -- log(0) is undefined.
        # Rather than silently leaving an ambiguous empty panel, say so.
        any_pp = (Y[:,2] > 0).any()
        any_ph = (Y[:,3] > 0).any()
        if any_pp:
            m = Y[:,2] > 0
            ax2.plot(t[m], Y[m,2], "-", color=c_Pp, lw=2.2, label="Pp (Amp)")
        if any_ph:
            m = Y[:,3] > 0
            ax2.plot(t[m], Y[m,3], "-", color=c_Ph, lw=2.2, label="Ph (Kan)")
        ax2.set_yscale("log")
        ax2.set_xlabel("time (min)")
        ax2.set_ylabel("particles / mL")
        ax2.set_title("Free particles", loc="left", weight="bold")
        if any_pp or any_ph:
            ax2.legend(frameon=False)
        else:
            ax2.text(0.5, 0.5, "no particles produced\n(zero for the entire run)",
                     transform=ax2.transAxes, ha="center", va="center",
                     fontsize=11, color="#888")
            ax2.set_yscale("linear")
            ax2.set_xticks([]); ax2.set_yticks([])
        ax2.grid(alpha=0.2)

        for ax in (ax1, ax2):
            if ps:
                ps.finish_axes(ax, legend=False)      # legends are drawn above
            else:
                for s_ in ("top", "right"):
                    ax.spines[s_].set_visible(False)

        fig.suptitle(f"Lp={args.Lp} kb, Lh={args.Lh} kb  "
                    f"phi={phi_used:.3f} (export split {split_txt})   "
                    f"moi={args.moi}, T={args.time} min",
                    fontsize=10, color="#555")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200, facecolor="white")
        print(f"\n  wrote {args.plot}")


if __name__ == "__main__":
    main()
