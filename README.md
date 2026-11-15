# SymPleco Phage Refinery Model

Princeton iGEM 2026. Models the phage refinery — a Jumpstart infection followed
by repeated Phase 2 rounds of infection and dilution — and predicts how the
Amp:Kan colony ratio (pHK002.2 vs M13K07) moves round over round.

Full write-up: `MODEL_REPORT.md`, and the wiki at
https://2026.igem.wiki/princeton/model

## Requirements

Python 3 and numpy. matplotlib only for the `--plot` flags.

**Nothing else — scipy is deliberately not used.** The ODEs are integrated with
a hand-written RK4 step, the Jeffreys intervals and the censored regression in
the data-analysis files compute their own normal CDF via `math.erf`. Everything
here runs on a plain numpy install.

```
python3 recursion.py          # the whole assay, defaults
python3 test_model.py         # 52 checks, ~80 s (--fast for ~8 s)
```

## Files

| File | What it is |
| --- | --- |
| `recursion.py` | Orchestrates the whole chain. **Start here.** No model logic of its own |
| `jumpstart.py` | The Jumpstart stage, 4 states, fixed packaging split `phi` |
| `jumpstart_mechanistic.py` | Alternative Jumpstart: genomes competing to replicate, `phi(t)` derived rather than assumed |
| `phase2.py` | One Phase 2 round, 7 states, the coat-competition core |
| `passage.py` | The 500 µL → 5 mL dilution, Poisson-sampled per pool |
| `plotstyle.py` | Shared palette and axis styling |
| `test_model.py` | 52 automated checks: invariants, equations, regressions |

Data analysis, separate from the model itself:

| File | What it is |
| --- | --- |
| `compare_data.py` | Raw per-square plate counts, R1–R3, both arms |
| `labdata.py` | The curated assay summary as plottable reference points |
| `trendfit.py` | Censored (Tobit) fit of the measured decline |
| `Assay_Data_.xlsx` | The source spreadsheet `labdata.py` transcribes |

`figures/` holds the seven report figures, regenerated from this exact code.

## Common commands

```
python3 recursion.py -r 10 --plot chain.png --csv chain.csv
python3 recursion.py --band               # a range across the measured Jumpstart split (3-16)
python3 recursion.py --mechanistic        # derived phi instead of fixed
python3 recursion.py --overlay-data       # measured plate ratios on the ratio panel
python3 recursion.py --fit                # + a censored fit through those points
python3 recursion.py --compare fig.png    # presentation figure: measured beside model
python3 jumpstart_mechanistic.py --sensitivity   # the two falsifiable predictions
```

## Three things to know before quoting a number

**The defaults are a trend preset, not a fit.** They were chosen to reproduce
the *shape* of the plate data — Amp:Kan falling every round — not to minimise
error against the points. Only `jump_ratio` comes from our own measurement.

**The Jumpstart is a calibration input, not a prediction.** `jump_ratio = 4.88`
was set from the measured R1 plate. The model starts there; the genuine test is
round 1 onward.

**`eps` is unmeasured and it is the largest lever** (14.3× across a plausible
0.02–0.2). The model is mechanistically grounded and quantitatively uncertain;
`origin_penalty` and the mechanistic Jumpstart's `alpha` are the other two real
unknowns. `MODEL_REPORT.md` §5 ranks every parameter by how much it moves the
answer, and §9 lists what the model does not reproduce — above all the ~10×
titre collapse the plates show by round 3.
