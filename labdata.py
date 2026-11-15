"""
labdata.py -- the measured plate ratios, as reference markers for the plots.

PURPOSE: show WHERE THE MODEL WAS ANCHORED, not score how well it fits.
recursion.py's defaults are a TREND PRESET -- chosen to reproduce the shape of
the plate data, not fitted to it -- so these markers are context for that
choice. Do not read closeness to a marker as validation.

SOURCE
------
Transcribed from `Assay_Data_.xlsx` (Sheet1), the curated assay summary. That
sheet, not the raw squares, is the authority here: it already resolves the
dilution weighting, carries 95% intervals, and gives explicit BOUNDS for the
plates that cannot yield a point estimate.

The table is embedded as a literal rather than read from the .xlsx on purpose.
Reading it would need openpyxl, and this model is deliberately dependency-light
(numpy only -- scipy is not even installed beside it). Thirteen rows transcribe
cleanly and stay readable and auditable in place. Check them against the sheet
if it is ever revised.

compare_data.py still holds the RAW per-square counts for run A and remains the
place to re-derive anything from scratch. Its pooled ratios differ slightly from
this sheet (e.g. 4.97 vs 5.08 at the Jumpstart, 0.722 vs 0.701 at round 1)
because the sheet weights dilutions differently. Where they disagree, the sheet
wins -- it is the reconciled version.

ROUND ALIGNMENT -- the thing most easily got wrong
--------------------------------------------------
The model's Jumpstart IS the first plated stage. The sheet names its stages
"Jumpstart", "Phase 2 rd 1", "Phase 2 rd 2", which map straight onto a
recursion results list:

    sheet stage     lab plate   results index
    Jumpstart       R1          0
    Phase 2 rd 1    R2          1
    Phase 2 rd 2    R3          2

So lab round N sits at index N-1, NOT index N. Getting this off by one slides
every marker a full round sideways, which on a curve falling ~10x per round
changes the apparent agreement by an order of magnitude.

recursion.py's own preset confirms it: jump_ratio defaults to 4.88, documented
as "the Jumpstart's PLATE Amp:Kan (measured: about 5)", and the sheet's
Jumpstart HP17KO7 value is 5.08. The model starts FROM that number rather than
predicting it -- which is what the hollow marker in the plot records.

TWO RUNS, SEVEN CONSTRUCTS -- READ THIS BEFORE QUOTING THE SPREAD
------------------------------------------------------------------
Run A is the N010/N011 pair (Sort3 vs the untagged HP17KO7 control). Run B is a
separate experiment carrying five further coat variants: Strep3, Sort8,
Sort8Sort3, Sort8Strep3, plus its own Sort3 and HP17KO7 arms.

The plot draws every one of these as a SINGLE undifferentiated series. That is
deliberate -- it answers "what did we measure" -- but it means several points at
one round are DIFFERENT CONSTRUCTS IN DIFFERENT EXPERIMENTS, not replicates of
one measurement. The vertical spread at Phase 2 rd 1 spans 0.70 down to <0.016,
and almost all of that is construct-to-construct and run-to-run difference, not
measurement error. Quote it as a range of observed values; never as precision,
and never as an error bar on a single quantity. `run` and `construct` are kept
on every point so the split can be restored without re-deriving it.

CENSORED VALUES
---------------
Four rows are one-sided. "> 63" at the Jumpstart is a plate whose Kan count was
zero in every square, so the ratio has a floor but no finite value; the "<"
rows are plates whose Amp count was zero. These are real results and are drawn
as arrows at the bound, not dropped and not silently turned into points. An
omitted bound reads as data never collected, which would be wrong.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Transcribed from Assay_Data_.xlsx, Sheet1, in sheet order.
#   stage       : the sheet's stage name
#   construct   : coat variant, as the sheet labels it
#   run         : 'A' or 'B'
#   value       : the number, or the bound for a censored row
#   bound       : None for a point estimate
#                 'lower' for a "> x" row  (true value is ABOVE value)
#                 'upper' for a "< x" row  (true value is BELOW value)
#   lo, hi      : 95% interval, or None where the sheet gives "--"
# ---------------------------------------------------------------------------
ASSAY_ROWS = [
    ("Jumpstart",    "HP17KO7",        "A",  5.08,  None,    3.21, 10.5),
    ("Jumpstart",    "Sort3",          "A", 63.0,   "lower", None, None),
    ("Jumpstart",    "all constructs", "B",  1.53,  None,    1.00, 2.80),
    ("Phase 2 rd 1", "HP17KO7",        "A",  0.701, None,    0.543, 1.06),
    ("Phase 2 rd 1", "Sort3",          "A",  0.44,  None,    0.139, 0.715),
    ("Phase 2 rd 1", "Sort3",          "B",  0.295, None,    0.062, 0.810),
    ("Phase 2 rd 1", "Strep3",         "B",  0.068, None,    0.012, 0.159),
    ("Phase 2 rd 1", "Sort8Strep3",    "B",  0.034, "upper", None, None),
    ("Phase 2 rd 1", "Sort8Sort3",     "B",  0.031, "upper", None, None),
    ("Phase 2 rd 1", "HP17KO7",        "B",  0.021, "upper", None, None),
    ("Phase 2 rd 1", "Sort8",          "B",  0.016, "upper", None, None),
    ("Phase 2 rd 2", "Sort3",          "A",  0.351, None,    0.182, 0.550),
    ("Phase 2 rd 2", "HP17KO7",        "A",  0.146, "upper", None, None),
]

# Sheet stage name -> index into a recursion results list. See ROUND ALIGNMENT.
STAGE_INDEX = {
    "Jumpstart":    0,
    "Phase 2 rd 1": 1,
    "Phase 2 rd 2": 2,
}

# The single value jump_ratio was calibrated from. Still flagged on every
# point as `is_anchor`, though the figures no longer distinguish it visually
# -- see the note in markers().
ANCHOR = ("Jumpstart", "HP17KO7", "A")

# Rows present in the sheet but kept OFF the figures, with the reason. They
# stay in ASSAY_ROWS so the transcription still matches the spreadsheet
# row-for-row; filtering happens at plot time instead of by deletion, so
# nothing silently disappears from the record.
#
#   ("Jumpstart", "all constructs", "B") -- a pooled figure over every
#   construct in run B, not a lineage of its own. It cannot be followed
#   across rounds like the others and sits oddly beside them, since it is an
#   average of things plotted separately at Phase 2 rd 1.
EXCLUDE_FROM_PLOT = {
    ("Jumpstart", "all constructs", "B"),
}

# ---------------------------------------------------------------------------
# Colour per LINEAGE (construct + run), not per point.
#
# Thirteen arbitrary colours would be thirteen things to look up. Colouring by
# lineage instead gives each curve-of-one-sample its own colour, so the two
# lineages that span all three rounds -- HP17KO7 run A and Sort3 run A -- can
# be followed across the plot, and the big spread at Phase 2 rd 1 is revealed
# as what it actually is: different constructs, not scatter on one.
#
# Okabe-Ito based, which stays distinguishable for the common forms of colour
# blindness. Black and grey are deliberately absent -- the model's two curves
# own those -- and the pale yellow of that palette is dropped as too faint on
# white. Assignment is by first appearance in ASSAY_ROWS, so it is stable.
# ---------------------------------------------------------------------------
LINEAGE_PALETTE = [
    "#D55E00",   # vermillion
    "#0072B2",   # blue
    "#009E73",   # bluish green
    "#CC79A7",   # reddish purple
    "#E69F00",   # orange
    "#56B4E9",   # sky blue
    "#7C3AED",   # violet
    "#A52A2A",   # brown
    "#B8860B",   # dark goldenrod
]


def lineage_key(m):
    """The grouping a colour is assigned to: one construct in one run."""
    return (m["construct"], m["run"])


def lineage_label(m):
    """Readable name for that grouping, e.g. 'Sort3 (run A)'."""
    return f"{m['construct']} (run {m['run']})"


def lineage_colours():
    """
    Map lineage key -> colour, in order of first appearance.

    If there were ever more lineages than palette entries the colours would
    repeat, so this raises instead -- silently reusing a colour would make two
    different constructs look like one.
    """
    keys = []
    for stage, construct, run, *_ in ASSAY_ROWS:
        if (construct, run) not in keys:
            keys.append((construct, run))
    if len(keys) > len(LINEAGE_PALETTE):
        raise ValueError(
            f"labdata: {len(keys)} lineages but only "
            f"{len(LINEAGE_PALETTE)} palette colours. Add more to "
            "LINEAGE_PALETTE rather than letting colours repeat.")
    return {k: LINEAGE_PALETTE[i] for i, k in enumerate(keys)}


def markers():
    """
    Every measured ratio, from both runs and all constructs, in sheet order.

    Colour distinguishes lineage; shape does not. Every point is drawn
    identically, INCLUDING the anchor and the six one-sided bounds, so
    nothing on a marker says whether it is a measurement, a bound, or the
    value the model was calibrated from. `is_anchor` and `bound` are still
    carried here, and a caption that does not mention them will mislead.

    Rows in EXCLUDE_FROM_PLOT are left out; `excluded()` reports them.

    Each point is a dict:
        stage      sheet stage name
        construct  coat variant
        run        'A' or 'B'
        index      position in a recursion results list
        value      the ratio, or the bound for a censored row
        bound      None | 'lower' | 'upper'
        lo, hi     95% interval, or None
        is_anchor  True for the one point jump_ratio was set from
    """
    out = []
    for stage, construct, run, value, bound, lo, hi in ASSAY_ROWS:
        if (stage, construct, run) in EXCLUDE_FROM_PLOT:
            continue
        if stage not in STAGE_INDEX:
            # A stage the model has no counterpart for should be noticed, not
            # quietly skipped -- a renamed stage in a revised sheet lands here.
            raise KeyError(
                f"labdata: stage {stage!r} has no entry in STAGE_INDEX. "
                "If Assay_Data_.xlsx added a stage, map it there.")
        out.append({
            "stage": stage,
            "construct": construct,
            "run": run,
            "index": STAGE_INDEX[stage],
            "value": value,
            "bound": bound,
            "lo": lo,
            "hi": hi,
            "is_anchor": (stage, construct, run) == ANCHOR,
            "lineage": (construct, run),
            "lineage_label": f"{construct} (run {run})",
        })
    return out


def excluded():
    """Rows deliberately kept off the figures, as (stage, construct, run)."""
    return sorted(EXCLUDE_FROM_PLOT)


def unplottable():
    """
    Kept for callers that report dropped points. Nothing is dropped now: the
    censored rows are drawn as bounds rather than discarded, so this is empty.
    """
    return []


if __name__ == "__main__":
    print(__doc__)
    print("POINTS (all drawn as one series; construct and run not "
          "distinguished)\n")
    hdr = f"  {'stage':<14} {'construct':<16} {'run':<4} {'Amp:Kan':>10}  interval"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for m in markers():
        if m["bound"] == "lower":
            val = f"> {m['value']:g}"
        elif m["bound"] == "upper":
            val = f"< {m['value']:g}"
        else:
            val = f"{m['value']:g}"
        iv = "--" if m["lo"] is None else f"{m['lo']:g} - {m['hi']:g}"
        tag = "   [ANCHOR for jump_ratio]" if m["is_anchor"] else ""
        print(f"  {m['stage']:<14} {m['construct']:<16} {m['run']:<4} "
              f"{val:>10}  {iv}{tag}")
    for stage, construct, run in excluded():
        print(f"\n  excluded from figures: {construct} (run {run}) at "
              f"{stage} -- pooled over constructs, not a lineage")
    n_b = sum(1 for m in markers() if m["bound"])
    print(f"\n  {len(markers())} points, {n_b} of them one-sided bounds.")
    print("  Bounds are drawn exactly like measured values -- nothing on the "
          "marker\n  distinguishes them, so say so in the caption.")
