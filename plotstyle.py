"""
Shared plot styling for the SymPleco model figures.

WHY THIS IS ITS OWN FILE
---------------------------
jumpstart.py, phase2.py and recursion.py all draw figures, and
before this they each carried their own colours and rcParams. That meant
three places to edit for any style change, and they had already drifted
apart. Putting it here means one edit changes every figure, and the
species colours stay consistent across files -- phagemid is always the
same blue whether you are looking at a single round or the whole chain.

USAGE
-----
    from plotstyle import apply_style, finish_axes, PHAGEMID, HELPER, ...

    apply_style()                       # once, before creating the figure
    fig, ax = plt.subplots(...)
    ax.plot(..., color=PHAGEMID, ...)
    finish_axes(ax)                     # once per axes, after plotting

Every function here is cosmetic only -- nothing in this module can change
a number the model produces.
"""
import matplotlib as mpl


# ==============================================================================
# PALETTE. One colour per biological entity, used identically everywhere.
# Chosen for reasonable contrast in both colour and greyscale, since these
# figures may end up printed or in a poster.
# ==============================================================================
PHAGEMID = "#2563eb"   # pHK002.2 / Amp -- blue
HELPER   = "#e11d48"   # M13KO7 / Kan   -- rose
NAIVE    = "#94a3b8"   # uninfected     -- slate, deliberately recessive
DUAL     = "#0f172a"   # dual-infected  -- near-black, the "end state"
ACCENT   = "#ea580c"   # events (e.g. Kan addition) -- amber
INK      = "#1e293b"   # titles, strong text
MUTED    = "#64748b"   # secondary text, captions
GRID     = "#e2e8f0"   # gridlines


def apply_style():
    """
    Set global matplotlib defaults. Call once before building a figure.

    Deliberately conservative: a clean sans-serif, generous whitespace,
    and gridlines light enough to read past. pdf.fonttype 42 embeds real
    fonts rather than outlines, so text stays selectable and searchable if
    a figure goes into a document.
    """
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 10.5,
        "axes.labelsize": 11,
        "axes.labelcolor": MUTED,
        "axes.titlesize": 13,
        "axes.titlecolor": INK,
        "axes.titleweight": "bold",
        "axes.titlepad": 12,
        "axes.edgecolor": "#cbd5e1",
        "axes.linewidth": 1.0,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.frameon": False,
        "legend.fontsize": 9.5,
        "lines.linewidth": 2.2,
        "lines.markersize": 6.5,
        "lines.markeredgewidth": 0.9,
        "lines.solid_capstyle": "round",
        "grid.color": GRID,
        "grid.linewidth": 0.9,
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })


def finish_axes(ax, legend=True, legend_loc="best"):
    """
    Tidy one axes after plotting: drop the top/right frame, put a light
    grid behind the data, and soften what remains.

    legend=False for axes that have nothing to label (e.g. a single line
    whose meaning is already in the title), which keeps the figure from
    carrying redundant keys.
    """
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cbd5e1")

    ax.grid(True, which="major", linestyle="-", alpha=0.7)
    ax.set_axisbelow(True)          # data draws OVER the grid, not under
    ax.tick_params(length=0)        # ticks are redundant with gridlines

    if legend and ax.get_legend_handles_labels()[0]:
        ax.legend(loc=legend_loc, borderpad=0.4, handlelength=1.8,
                  labelcolor=INK)


def marker_style(color):
    """
    Consistent marker look: filled with the line colour, ringed in white so
    overlapping points stay distinguishable on a dense log axis.
    """
    # linestyle is set EXPLICITLY here: without it, any caller that also
    # passes a style string risks the line rendering oddly. And the edge
    # width must stay BELOW the line width -- a wider ring punches a white
    # gap through the line at every marker, which reads as a dashed line
    # even though it is solid.
    return dict(marker="o", linestyle="-", markerfacecolor=color,
                markeredgecolor="white", markeredgewidth=0.9)


def caption(fig, text):
    """
    Put a small muted note under the figure -- used to record the exact
    parameter set a figure came from, so a saved PNG can always be traced
    back to the run that produced it.
    """
    fig.text(0.5, -0.02, text, ha="center", va="top",
             fontsize=9, color=MUTED)
