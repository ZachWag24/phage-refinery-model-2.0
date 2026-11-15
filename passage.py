"""
Passage: the transfer between rounds.

500 uL of supernatant goes into 5 mL of fresh culture -- an 11x dilution.
Standalone piece, used between EVERY pair of rounds (Jumpstart -> Phase 2,
and later Phase 2 -> Phase 2), so it lives in its own file rather than
being bundled into whichever stage happens to produce its input.

THE OBVIOUS WAY TO DO THIS IS WRONG AT LOW COUNTS
--------------------------------------------------
The naive approach: multiply every concentration by (vol_transfer/vol_new).
That's exact when there are huge numbers of particles, but breaks down once
counts get low -- you cannot transfer "0.3 of a particle", and dividing by
11 will never produce a TRUE zero no matter how few particles remain, which
means a species can never actually go extinct in a naive model. It just
decays toward an arbitrarily tiny concentration that keeps re-growing next
round, which is not physical.

THE FIX: POISSON SAMPLING
--------------------------
Pipetting a fixed volume out of a well-mixed suspension is a random draw --
you get whatever happens to be in that 500 uL, not a precise fraction of
the whole. The number of particles that end up in the transferred volume
follows a Poisson distribution with mean = concentration x volume:

    n_transferred ~ Poisson(concentration * vol_transfer)

This can genuinely land on zero. Once a species is down to a handful of
particles in the transfer volume, that possibility is not a rounding
error -- it is the actual, correct answer, and it is exactly how a real
species can go extinct between rounds.

WHEN DOES THIS ACTUALLY MATTER?
---------------------------------
Only when the expected count in the transferred volume is small. At
1e12 particles/mL, 500 uL contains ~5e11 particles -- the Poisson spread
around that mean is a rounding error, so Poisson sampling and simple
division give indistinguishable answers. The two approaches only diverge
once a concentration has fallen far enough that the transferred volume
holds single-digit particles.
"""
import numpy as np


def passage(*concentrations, vol_transfer=0.5, vol_new=5.5, rng=None):
    """
    Transfer vol_transfer mL of supernatant into vol_new mL of fresh
    culture, sampling EACH pool independently.

    *concentrations : one or more concentrations (particles/mL) BEFORE
                    transfer. Each is sampled separately, because each is a
                    physically distinct pool of particles that dilutes on
                    its own. Originally this took exactly two (Pp, Ph);
                    generalised so coat state can be carried across rounds
                    as four pools (Pw, Pt, Hw, Hh) -- phagemid and helper,
                    each wild-type or made in Phase 2 -- without needing a
                    second near-identical function.
    vol_transfer  : volume taken from the old culture, mL
    vol_new       : TOTAL volume after adding fresh culture, mL (so the
                    dilution factor is vol_transfer/vol_new, e.g. 0.5/5.5)
    rng           : a numpy random generator. Pass one you created yourself
                    (e.g. np.random.default_rng(seed)) to get a reproducible
                    run; omit it for a fresh random draw each call.

    Returns a tuple of the same length as the input, with concentrations
    AFTER transfer and dilution.

    NOTE ON SAMPLING EACH POOL SEPARATELY: this is the physically correct
    thing to do -- the particles really are distinct objects and a pipette
    really does draw each independently. It is NOT the same as sampling a
    combined total and then splitting it proportionally, which would
    understate how much the small pools fluctuate.
    """
    if rng is None:
        rng = np.random.default_rng()

    out = []
    for conc in concentrations:
        mean_n = conc * vol_transfer   # expected number of particles in the pipetted volume
        # Poisson sampling only matters once counts are small enough for the
        # randomness to show (see module docstring). At a mean of 1e7 the
        # relative spread is already 3e-4, so above that the two approaches
        # are indistinguishable and this skips straight to the mean. That is a
        # convenience, not a necessity: numpy can sample large means fine.
        n = rng.poisson(mean_n) if mean_n < 1e7 else mean_n
        out.append(n / vol_new)        # back to a concentration in the new total volume

    return tuple(out)


if __name__ == "__main__":
    print("PASSAGE, standalone demonstration\n")

    print("1. At high particle counts, Poisson sampling matches simple")
    print("   division almost exactly -- no meaningful difference.")
    print("   (1.7e12 is above passage()'s 1e7 shortcut, so it just divides; 1.9e7 is")
    print("   below it, so passage() really samples. A REAL Poisson draw is shown for both.)")
    rng = np.random.default_rng(0)
    for conc in (1.7e12, 1.9e7):
        naive = conc * 0.5 / 5.5
        real = rng.poisson(conc * 0.5) / 5.5          # a genuine Poisson draw
        via = passage(conc, rng=rng)[0]
        print(f"   conc={conc:.1e}: real Poisson draw differs from division by "
              f"{abs(real-naive)/naive:.1e}; passage() by {abs(via-naive)/naive:.1e}")
    print()

    print("2. At low counts, extinction becomes a real, frequent outcome --")
    print("   not something a naive model could ever produce.")
    print(f"   {'conc/mL':>10} {'mean in 0.5 mL':>16} {'P(zero)':>10}")
    for conc in [1e12, 1e6, 100, 10, 4, 1]:
        mean_n = conc * 0.5
        p_zero = np.exp(-mean_n)
        print(f"   {conc:10.0e} {mean_n:16.3g} {p_zero*100:9.1f}%")

    print("\n3. Chaining Jumpstart output through passage(), previewing what")
    print("   Phase 2 will receive as its starting concentrations:")
    from jumpstart import run
    t, Y = run()
    B, I, Pp_final, Ph_final = Y[-1]
    print(f"   Jumpstart output : Pp={Pp_final:.4e}  Ph={Ph_final:.4e}  "
          f"(ratio {Pp_final/Ph_final:.3f})")
    rng = np.random.default_rng(1)
    Pp_in, Ph_in = passage(Pp_final, Ph_final, rng=rng)
    print(f"   after passage    : Pp={Pp_in:.4e}  Ph={Ph_in:.4e}  "
          f"(ratio {Pp_in/Ph_in:.3f})")
    print("   -> the ratio is preserved (both species diluted by the same")
    print("      factor); only the ABSOLUTE concentrations drop.")
