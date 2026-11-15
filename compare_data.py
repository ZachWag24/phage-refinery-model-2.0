r"""
compare_data.py -- the raw plate data, and what it does and doesn't tell us.
=============================================================================
Holds every colony count from the R1-R3 platings and turns them into Amp:Kan
ratios with confidence intervals.

    python3 compare_data.py

Needs numpy and scipy (for the confidence interval) only -- no model import,
no plotting. This file is a DATA ARCHIVE, not a fitting exercise: it does not
try to find a model parameter that reproduces these numbers. See
modelling_report.md Section 7 for why we deliberately did not fit ("use the
model to explore, not to match data").

WHAT N010 AND N011 ACTUALLY ARE
--------------------------------
This matters and was previously ambiguous in team discussion:

    N010 = HP17_SORT3  -- the SORT3-TAGGED helper plasmid. The arm carrying
                           the coat modification whose infectivity effect
                           (H2) we are trying to test.
    N011 = HP17K07     -- the UNTAGGED control. Same protocol, same rounds,
                           no sortase tag.

So the direct test of "does the tag reduce infectivity" is N010 vs N011 at
the same round, not N010 alone across rounds. See compare_tagged_vs_control()
below -- the result is mixed, not a clean confirmation (round 2 goes the
predicted direction, round 3 does not: the untagged control drops to 0.00,
harder than the tagged arm).

WHY THE RAW COUNTS LIVE HERE
-----------------------------
So every number quoted anywhere (report, wiki, Slack) traces back to actual
colonies on a plate. Change the inclusion rules below and see what it does
to the answer.

PLATE LAYOUTS
-------------
R1 (07/25): 6x6 grid. Phage diluted 10^0..10^-5 left to right, bacteria
            diluted top to bottom. 5 uL spotted per square.
R2 (07/27), R3 (07/29): 5x5 grid, different scheme. Rows are phage dilutions
            at fixed bacteria; columns are dilutions of the whole mixture.
            A "phage only, no bacteria" control row was run and came back
            clean (all zeros) on every plate.

Amp and Kan plates were spotted from the SAME mixture at the SAME dilution,
so matched (row, column) cells are directly comparable and their ratio is
the ratio of transducing units.

INCLUSION RULES
----------------
  - LAWN squares are excluded: confluent growth cannot be counted.
  - A cell is used only if Amp + Kan >= MIN_COUNTS (default 3). Below that,
    Poisson noise swamps the signal.
  - Ratios are pooled across included cells rather than averaged, which
    weights by information content rather than treating a 3-colony square
    the same as a 50-colony one.

KNOWN DATA PROBLEM -- R1
--------------------------
The R1 kanamycin plates failed for BOTH lineages -- N010 came back entirely
zero across the full grid, and N011 shows the same left-side dropout (zero
at the most concentrated corner, a peak in the middle, which is backwards).
Confirmed by the PI. R1 is reported here but excluded from every summary
statistic.

KNOWN DATA PROBLEM -- R2 to R3 titre collapse
-----------------------------------------------
Total colonies (Amp+Kan) fell 7.3x (N010) and 11.2x (N011) between rounds 2
and 3, in BOTH arms. Team discussion (as of the round-3 Slack thread) ruled
out incubation time as the cause -- plaque infection time is fixed at 10 min
with no incubation step, and R1/R2/R3 were run on different dates rather
than in the planned same-day format. The cause is not resolved. Because both
arms dropped together, and the untagged control dropped harder than the
tagged arm, this is NOT simply consistent with random per-plate error --
see compare_tagged_vs_control() below.
"""

import numpy as np

# scipy is used ONLY for the Jeffreys confidence intervals. The rest of this
# file -- the plate counts themselves, and the pooled ratios -- is plain
# arithmetic. The model is deliberately numpy-only and scipy is not installed
# alongside it, so importing scipy at module level would make the raw DATA
# dict unreachable from code that just wants the counts (labdata.py does).
# Import it lazily instead: the counts always load, and only the CI degrades.
try:
    from scipy import stats
    HAVE_SCIPY = True
except ImportError:
    stats = None
    HAVE_SCIPY = False

LAWN = None          # sentinel for an uncountable square
MIN_COUNTS = 3       # minimum Amp+Kan in a square for it to be used

TAGGED = 'N010'       # HP17_SORT3 -- carries the coat modification
UNTAGGED = 'N011'     # HP17K07    -- untagged control


# ==============================================================================
# THE RAW PLATE COUNTS
# Keys are (row, column) as they appear on the plate.
# ==============================================================================
DATA = {
    # ---------------- R1, 07/25/2026 -- KAN PLATES COMPROMISED --------------
    (TAGGED, 1): {
        'amp': {(0,0):26, (0,1):24, (0,2):18, (0,3):16, (0,4):14, (0,5):0,
                (1,0):20, (1,1):25, (1,2):20, (1,3):8,  (1,4):0,  (1,5):0,
                (2,0):5,  (2,1):5,  (2,2):6,  (2,3):0,
                (3,0):0,  (3,1):2,  (3,2):0},
        'kan': {(0,0):0,  (0,1):0,  (0,2):0,  (0,3):0,  (0,4):0,  (0,5):0,
                (1,0):0,  (1,1):0,  (1,2):0,  (1,3):0,  (1,4):0,  (1,5):0,
                (2,0):0,  (2,1):0,  (2,2):0,  (2,3):0,
                (3,0):0,  (3,1):0,  (3,2):0},
        'valid': False,
        'note': 'Kan plate entirely zero -- plating failure confirmed by PI',
    },
    (UNTAGGED, 1): {
        'amp': {(0,0):LAWN, (0,1):24, (0,2):23, (0,3):8, (0,4):7, (0,5):0,
                (1,0):30, (1,1):21, (1,2):15, (1,3):8,  (1,4):0,
                (2,0):9,  (2,1):18, (2,2):10, (2,3):3,
                (3,0):2,  (3,1):3,
                (4,2):2},
        'kan': {(0,0):0,  (0,1):5,  (0,2):2,  (0,3):0, (0,4):0, (0,5):0,
                (1,0):6,  (1,1):6,  (1,2):2,  (1,3):0, (1,4):0,
                (2,0):3,  (2,1):12, (2,2):0,  (2,3):0,
                (3,0):0,  (3,1):0,
                (4,2):0},
        'valid': False,
        'note': 'Kan plate shows the same left-side dropout: 0 at the most '
                'concentrated corner, peak in the middle',
    },

    # ---------------- R2, 07/27/2026 ---------------------------------------
    (TAGGED, 2): {
        'amp': {(0,0):22, (0,1):52, (0,2):1,  (0,3):0, (0,4):0,
                (1,0):45, (1,1):14, (1,2):4,  (1,3):0, (1,4):0,
                (2,0):29, (2,1):1,  (2,2):0,
                (3,0):4,  (3,1):1,
                (4,0):0},
        'kan': {(0,0):LAWN, (0,1):53, (0,2):12, (0,3):2, (0,4):0,
                (1,0):LAWN, (1,1):50, (1,2):12, (1,3):0, (1,4):0,
                (2,0):47,   (2,1):21, (2,2):3,
                (3,0):37,   (3,1):1,
                (4,0):3},
        'valid': True,
        'note': '',
    },
    (UNTAGGED, 2): {
        'amp': {(0,0):51, (0,1):37, (0,2):4, (0,3):0, (0,4):0,
                (1,0):41, (1,1):12, (1,2):2, (1,3):0, (1,4):0,
                (2,0):27, (2,1):5,  (2,2):2,
                (3,0):21, (3,1):4,
                (4,0):3},
        'kan': {(0,0):LAWN, (0,1):58, (0,2):4, (0,3):2, (0,4):0,
                (1,0):LAWN, (1,1):31, (1,2):4, (1,3):2, (1,4):1,
                (2,0):42,   (2,1):4,  (2,2):1,
                (3,0):15,   (3,1):0,
                (4,0):3},
        'valid': True,
        'note': '',
    },

    # ---------------- R3, 07/29/2026 ---------------------------------------
    (TAGGED, 3): {
        'amp': {(0,0):5, (0,1):0, (0,2):0,
                (1,0):4, (1,1):0, (1,2):1,
                (2,0):3, (2,1):0,
                (3,0):0,
                (4,0):0},
        'kan': {(0,0):12, (0,1):1, (0,2):0,
                (1,0):14, (1,1):0, (1,2):0,
                (2,0):5,  (2,1):0,
                (3,0):4,
                (4,0):1},
        'valid': True,
        'note': '',
    },
    (UNTAGGED, 3): {
        'amp': {(0,0):10, (0,1):0, (0,2):0, (0,3):0,
                (1,0):0,  (1,1):0,
                (2,0):6,  (2,1):0, (2,2):0,
                (3,0):0,  (3,1):0,
                (4,0):0},
        'kan': {(0,0):LAWN, (0,1):4, (0,2):0, (0,3):2,
                (1,0):LAWN, (1,1):3,
                (2,0):LAWN, (2,1):3, (2,2):3,
                (3,0):LAWN, (3,1):3,
                (4,0):9},
        'valid': True,
        'note': 'Amp counts only in column 0, where Kan is uncountable',
    },
}


# ==============================================================================
# TURNING COUNTS INTO A RATIO
# ==============================================================================
def pooled_ratio(entry, min_counts=MIN_COUNTS):
    """
    Pool matched squares and return (amp_total, kan_total, ratio, lo, hi, cells).

    Only squares where BOTH plates were countable and where Amp+Kan >=
    min_counts contribute. The confidence interval treats the Amp fraction
    as a binomial proportion with a Jeffreys prior, which behaves sensibly
    when one side is zero -- as it is for the untagged arm at R3.
    """
    amp, kan = entry['amp'], entry['kan']
    a_tot = k_tot = 0
    used = []
    for key in sorted(set(amp) & set(kan)):
        a, k = amp[key], kan[key]
        if a is LAWN or k is LAWN:
            continue
        if a + k < min_counts:
            continue
        a_tot += a
        k_tot += k
        used.append((key, a, k))
    if a_tot + k_tot == 0:
        return a_tot, k_tot, np.nan, np.nan, np.nan, used
    ratio = a_tot / k_tot if k_tot > 0 else np.inf
    if not HAVE_SCIPY:
        # The ratio is still exact; only the interval is unavailable. Returning
        # nan rather than raising keeps every caller that wants the counts
        # working on a machine without scipy.
        return a_tot, k_tot, ratio, np.nan, np.nan, used
    lo_f, hi_f = stats.beta.ppf([0.025, 0.975], a_tot + 0.5, k_tot + 0.5)
    lo = lo_f / (1 - lo_f) if lo_f < 1 else np.inf
    hi = hi_f / (1 - hi_f) if hi_f < 1 else np.inf
    return a_tot, k_tot, ratio, lo, hi, used


def summarise():
    """Print every lineage/round with its pooled ratio. Returns the valid rows."""
    print("=" * 74)
    print("PLATE DATA -> Amp:Kan RATIOS")
    print("=" * 74)
    print(f"Squares used only where both plates countable and Amp+Kan >= "
          f"{MIN_COUNTS}\n")
    rows = []
    for (lineage, rnd), entry in sorted(DATA.items(), key=lambda x: (x[0][1], x[0][0])):
        a, k, r, lo, hi, used = pooled_ratio(entry)
        tag = "" if entry['valid'] else "   [EXCLUDED]"
        rstr = f"{r:8.3f}" if np.isfinite(r) else "     inf"
        cistr = f"{lo:.3f}-{hi:.3f}" if np.isfinite(lo) and np.isfinite(hi) else "--"
        label = "tagged (SORT3)" if lineage == TAGGED else "untagged control"
        print(f"  {lineage} round {rnd}  [{label}]{tag}")
        print(f"     squares used : {len(used)}")
        print(f"     Amp total    : {a}")
        print(f"     Kan total    : {k}")
        print(f"     ratio        : {rstr} : 1     95% CI {cistr}")
        if entry['note']:
            print(f"     note         : {entry['note']}")
        print()
        if entry['valid']:
            rows.append((lineage, rnd, r, lo, hi, a, k))
    return rows


def titre_trend(rows):
    """Total colony counts per round -- is the whole system washing out?"""
    print("=" * 74)
    print("TOTAL COLONY COUNTS BY ROUND")
    print("=" * 74)
    print("  Both particle types together. A falling total means the system is")
    print("  losing phage overall, not just shifting composition.\n")
    print(f"  {'lineage':>8} {'round':>6} {'Amp':>6} {'Kan':>6} {'total':>7}")
    by_lineage = {}
    for lineage, rnd, r, lo, hi, a, k in rows:
        print(f"  {lineage:>8} {rnd:6d} {a:6d} {k:6d} {a+k:7d}")
        by_lineage.setdefault(lineage, []).append((rnd, a + k))
    print()
    for lineage, pts in by_lineage.items():
        pts.sort()
        if len(pts) >= 2:
            (r0, t0), (r1, t1) = pts[0], pts[-1]
            if t1 > 0:
                print(f"  {lineage}: total fell {t0/t1:.1f}x between rounds "
                      f"{r0} and {r1}")
            else:
                print(f"  {lineage}: total fell to zero by round {r1}")
    print()


def compare_tagged_vs_control(rows):
    """
    THE ACTUAL TEST OF H2 (does the tag reduce infectivity).

    N010 (tagged) vs N011 (untagged) at the SAME round is the direct
    comparison. If the tag hurts infectivity, tagged should read lower than
    untagged at every round. It does at round 2 and does not at round 3.
    """
    print("=" * 74)
    print("TAGGED (N010) vs UNTAGGED CONTROL (N011) -- the direct H2 test")
    print("=" * 74)
    print("  If the sortase tag reduces infectivity, tagged should read LOWER")
    print("  than the untagged control at every round.\n")
    by_round = {}
    for lineage, rnd, r, lo, hi, a, k in rows:
        by_round.setdefault(rnd, {})[lineage] = (r, a, k)
    print(f"  {'round':>6} {'tagged (N010)':>15} {'untagged (N011)':>17} "
          f"{'tag reads lower?':>18}")
    for rnd in sorted(by_round):
        d = by_round[rnd]
        if TAGGED in d and UNTAGGED in d:
            rt, at, kt = d[TAGGED]
            ru, au, ku = d[UNTAGGED]
            consistent = "yes" if rt < ru else "NO"
            print(f"  {rnd:6d} {rt:15.3f} {ru:17.3f} {consistent:>18}")
    print()
    print("  Round 2 is consistent with H2 (tagged reads lower).")
    print("  Round 3 is NOT: the UNTAGGED control drops to 0.00, harder than")
    print("  the tagged arm at 0.34. A control with no tag crashing worse than")
    print("  the tagged sample argues against infectivity being the driver of")
    print("  the round-3 result specifically. Two rounds is not enough to")
    print("  call this confirmed or refuted either way.")
    print()


if __name__ == "__main__":
    rows = summarise()
    titre_trend(rows)
    compare_tagged_vs_control(rows)
