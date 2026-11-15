# Phage Refinery Model — Report

## 1. Summary

This model simulates the phage refinery process — a starter Jumpstart infection followed by repeated Phase 2 rounds of infection and dilution — and predicts how the Amp:Kan colony ratio (pHK002.2 vs. M13K07) changes round over round. It reproduces the trend seen on our plates: the ratio starts above 1 and falls every round as the engineered coat's infectivity advantage compounds, flipping below 1 already within the first Phase 2 round.

**Files:** `jumpstart.py`, `jumpstart_mechanistic.py` (an alternative account of the Jumpstart, §8), `passage.py`, `phase2.py`, `recursion.py` (orchestrates the whole chain), `plotstyle.py`, `test_model.py` (52 automated checks).

**Entry point:** `python3 recursion.py` runs the whole assay. See §13 for the full command reference.

---

## 2. The experiment being modeled

**Part 1 — Jumpstart (one time).** A fresh pHK002.2 colony is grown in LB with ampicillin to A₆₀₀ ≈ 0.05, superinfected with M13K07 helper phage, shaken 60–90 min, then kanamycin is added and the culture grown overnight. Output: an overnight culture containing secreted particles carrying either the pHK002.2 or the M13K07 genome — both still wearing wild-type coat, since no engineered coat exists yet at this stage.

**Part 2 — Purification, titering, next round (repeats).** Centrifuge the culture, then from the same supernatant: (A) titer 500 µL on Amp and Kan plates against wild-type ER2738, (B) a negative-control plate (no host, confirms no bacterial carry-over), and (C) feed 500 µL into fresh ER2738 carrying HP17K07 for the next round.

**The readout** is the Amp:Kan colony ratio each round.

### The three genetic elements

| Element | What it is | Key facts |
|---|---|---|
| **pHK002.2** | The phagemid. Carries the M13/f1 origin and AmpR on a pBR322-type backbone. | Carries no M13 genes — cannot make particles alone. ~3.8 kb. Resident copy number ~15–20/cell (confirmed via its Benchling map: `ori` + `rop`, the pBR322 copy-number-repression signature). |
| **M13K07** | Helper phage. | Complete phage genome; its f1 origin is interrupted by the KanR insertion — a documented, designed feature of helper phages, causing it to preferentially replicate any co-resident wild-type f1 origin (like pHK002.2's) over its own. ~8.7 kb. Its replicative-form DNA reaches ~130–300 copies/cell once established. |
| **HP17K07** | Helper plasmid resident in every Phase 2 cell. | Supplies all assembly machinery and the engineered (sortase-tagged) coat — modified pIII and pVIII. Cannot itself be packaged. Backbone confirmed via Benchling: `Rep101`/pSC101, ~5–12 copies/cell. |

### Coat states

| Particle | Where made | Coat |
|---|---|---|
| Jumpstart pHK002.2 / M13K07 | Jumpstart (no HP17K07 present) | wild-type |
| Phase 2 pHK002.2 | `Ip` cell (HP17K07 + phagemid) | fully modified |
| Phase 2 M13K07 | `Ih` cell (HP17K07 + M13K07) | hybrid |

A cell already infected by one phage essentially never gets infected by the other (per the user's assessment), so Phase 2 has two parallel lineages, never a dual-infected state.

**Central mechanism:** the hybrid coat infects better than the fully modified coat. That's why the Amp:Kan ratio inverts.

---

## 3. Model structure

```
jumpstart.run() ──passage()──► phase2.run() ──passage()──► phase2.run() ── ...
   (once)                          (round 1)                   (round 2)
```

`recursion.py` orchestrates; it contains no model logic of its own.

**`jumpstart.py` — 4 states** (`B`: phagemid-only cells, `I`: dual-infected cells, `Pp`/`Ph`: free particles). A cell starts with just the phagemid and becomes dual-infected once it catches a helper particle; a shared export capacity, split by `phi`, produces both particle types. Kanamycin kills every remaining uninfected cell at `kan_at=75` min.

**`passage.py`** — 500 µL into 5 mL (5.5 mL total, a 1:11 dilution), each of the four particle pools sampled independently by Poisson (not plain division), so a species can genuinely reach exact extinction rather than asymptoting toward zero forever.

**`phase2.py` — 7 states** (`N`: naive, `Ip`/`Ih`: infected, `Pw`/`Pt`: phagemid wild-type/modified, `Hw`/`Hh`: helper wild-type/hybrid). Each round, naive cells get captured by whichever particle type is more concentrated × more infective; the plate readout weights every pool by its own infectivity.

**`recursion.py`** — runs the Jumpstart once, then N Phase 2 rounds, passaging all four pools between stages. `--band` reports a range across the measured Jumpstart split (3–16) instead of a single line; `--mechanistic` swaps in `jumpstart_mechanistic.py` (§8) for the Jumpstart stage.

---

## 4. The core equations

**Jumpstart (4 equations).** `sat` caps how fast infection can happen — a cell only has a few pili, no matter how many particles are around. `room` caps population growth at carrying capacity `K`.

```
sat  = 1 / (1 + k1·h·Ph)
inf  = k1·Ph·B·sat
room = 1 − (B+I)/K

dB/dt  = −inf + k3·B·room
dI/dt  =  inf + k4·I·room
dPp/dt =  φ·k2·I
dPh/dt = −inf + (1−φ)·k2·I
```

`k4 = growth_penalty × k3` (an infected cell divides slower — a real, measured cost of devoting resources to making phage). `k2` is the shared export rate, split by `φ`.

**Phase 2 (7 equations).** Each particle type infects at a rate set by its own weight (`a_w`, `a_t`, `1`, `a_h` — see §5). `sink` is the shared removal rate every pool draws from, scaled by its own weight.

```
sat = 1 / (1 + k1·h·(a_w·Pw + a_t·Pt + Hw + a_h·Hh))
infect_p = k1·(a_w·Pw + a_t·Pt)·sat
infect_h = k1·(Hw + a_h·Hh)·sat
room = 1 − (N+Ip+Ih)/K

dN/dt  = −(infect_p+infect_h)·N + k3·N·room
dIp/dt =  infect_p·N + k4·Ip·room
dIh/dt =  infect_h·N + k4·Ih·room

sink = k1·sat·(N + σ·(Ip+Ih))
dPw/dt = −sink·a_w·Pw
dPt/dt = −sink·a_t·Pt + k2p·Ip
dHw/dt = −sink·Hw
dHh/dt = −sink·a_h·Hh + k2h·Ih
```

`σ` (sigma) is superinfection exclusion — an already-infected cell absorbs further particles at only ~5% the rate of a naive cell, but still absorbs them (this term is identical in form for every particle type; a cell doesn't know what genome an incoming particle carries). `k2p = k2 × (Lh/Lp)^size_effect`; `k2h = k2 × origin_penalty` (§5).

**The hybrid coat's infectivity, `a_h`, is derived from real coat-protein gene dosage:**

```
wt_frac = m13ko7_copies / (m13ko7_copies + hp17_copies)
a_h     = wt_frac + (1 − wt_frac) × eps
```

An `Ih` cell carries HP17K07's resident (modified) pIII/pVIII genes and an infecting M13K07's wild-type ones. pIII (only ~5 copies/virion, versus pVIII's ~2700) is what sets infectivity, so the wild-type share of pIII protein — and hence `a_h` — is set by gene dosage. Given roughly equal promoter strength for the two pIII genes, at defaults (`m13ko7_copies=20`, `hp17_copies=8`) this gives `a_h ≈ 0.729` — nearer wild-type strength than `eps`, but well short of it.

**Plate readout:**
```
Amp colonies ~ (Pw + eps·Pt)
Kan colonies ~  Hw + a_h·Hh
```

---

## 5. Parameters: every active one, with status

### Measured / protocol-derived

| Parameter | Default | Source |
|---|---|---|
| `vol_transfer`, `vol_new` | 0.5, 5.5 mL | SOP (500 µL into 5 mL), confirmed against the wiki process diagram |
| `T1` (Jumpstart length) | 1080 min (18 h) | "overnight" |
| `T2` (Phase 2 round) | 900 min (15 h) | "overnight" |
| `kan_at` | 75 min | SOP midpoint of the 60–90 min window |
| `B0_jump` | 4.0e7 cells/mL | A₆₀₀ ≈ 0.05 |
| `Lp`, `Lh` | 3.8, 8.7 kb | genome lengths |
| `jump_ratio` | 4.88 | measured Jumpstart plate ratio (range observed ≈3–16) |
| `k1` | 3.0e-11 mL/min | stated for this assay |
| `moi` | 2.5 | derived: M13K07 stock titer (1e11 pfu/mL) × the protocol's 1:1000 dilution → moi=2.498, matching to 3 sig figs |

### Assumed — and they matter (ranked by round-2 plate-ratio sensitivity)

| Parameter | Default | What it is | Spread over a plausible range |
|---|---|---|---|
| `eps` | 0.05 | fully modified coat's infectivity relative to wild-type | **15.1×** (0.02→0.2) |
| `jump_ratio` | 4.88 | the Jumpstart's measured plate ratio, bounded (3–16) | **4.3×** (3→16) |
| `vol_transfer` | 0.5 mL | dilution factor at passage | **3.9×** (0.25→1.0 mL) |
| `origin_penalty` | 1/1.05 ≈ 0.952 | multiplier on an `Ih` cell's export rate, for M13K07's interrupted origin costing it export speed even without direct competition. Deliberately tied to the mechanistic Jumpstart's `replication_ratio` — one guess about the origin defect's size, not two independent ones | **1.9×** (0.7→1.0) |
| `size_effect` | 0.0 | exponent on (Lh/Lp) for the phagemid's *production* rate (0 = no advantage) | 1.8× (0.5→1.0) |

### Assumed, from general literature — and no longer robust to where within range

| Parameter | Default | What it is | Spread over the *published* range |
|---|---|---|---|
| `m13ko7_copies` | 20 | M13K07's copy number via its p15A origin — **not** wild-type M13's f1 rolling-circle number. Published p15A range ~10–30; 20 is the end of it the measured round-1 data favours | **1.73×** (10→30) |
| `hp17_copies` | 8 | HP17's own copy number, Rep101/pSC101 (published range ~5–12) | **1.47×** (5→12) |

These were robust at the old (incorrect) `m13ko7_copies=200`, where M13K07's copy number swamped HP17's whatever its exact value. Now that the two are comparable magnitude (~20 vs ~8), **both genuinely move the result**. That is a real consequence of fixing the bigger error: the corrected model is more parameter-sensitive, not less.

### Assumed — and they barely matter

`sigma` (0.05, exclusion), `K` (1e9 cells/mL), `k2` (3.0 particles/cell/min), `h` (1.0 min), `k3` (ln2/25 min), `growth_penalty` (0.6 — an infected cell's division rate as a fraction of an uninfected one's), `B0_jump`, `B0_phase2` (3e8, still just an assumed default with no protocol-derived anchor), `kan_at` — each moves the round-2 plate ratio by ≤1.2× across a wide range.

---

## 6. Current results

Default chain (`python3 recursion.py -r 8`):

| Stage | Pp (/mL) | Ph (/mL) | Tube Amp:Kan | **Plate Amp:Kan** |
|---|---|---|---|---|
| Jumpstart | 2.14e12 | 4.39e11 | 4.88 | **4.88** |
| round 1 | 2.23e12 | 4.27e11 | 5.23 | **0.870** |
| round 2 | 1.35e12 | 1.28e12 | 1.06 | **0.0871** |
| round 3 | 3.20e11 | 2.25e12 | 0.142 | **0.0105** |
| round 4 | 5.46e10 | 2.51e12 | 0.0218 | **0.00155** |
| round 5 | 8.79e9 | 2.55e12 | 0.00345 | **0.000241** |
| round 6 | 1.39e9 | 2.56e12 | 0.000545 | **3.78e-05** |
| round 7 | 2.20e8 | 2.56e12 | 8.60e-05 | **5.94e-06** |
| round 8 | 3.47e7 | 2.56e12 | 1.36e-05 | **9.33e-07** |

The phagemid goes extinct at round 18 (seed 1). The helper never does.

**Round 1 already flips below 1, within the round itself.** Amp's round-1 plate count is 63.2% carried-over wild-type phagemid; Kan's is 88.7% freshly-made hybrid particles, because the hybrid coat is weighted at 73% strength. The newly-made hybrid particles dominate Kan's own count already within round 1.

**Why the ratio keeps falling:** each round, the better-infecting hybrid helper captures more of the fresh cells, so more cells make helper, so it enters the next round in greater proportion — compounding. One hybrid helper particle counts **14.6×** a fully modified phagemid particle on a plate (`a_h/eps`).

**Absolute titers are not a prediction.** Total particles per round (Pp+Ph) settle near ~2.56e12/mL by round 4 onward — nowhere near the ~10× collapse the real plate data show by round 3. See §9, limitation 1.

---

## 7. Sensitivity (round-2 plate ratio, baseline 0.0600)

| Parameter | Range tested | Low / high | Spread |
|---|---|---|---|
| eps | 0.02 / 0.2 | 0.0261 / 0.3932 | 15.1× |
| jump_ratio | 3 / 16 | 0.0401 / 0.1728 | 4.3× |
| vol_transfer | 0.25 / 1.0 mL | 0.0347 / 0.1340 | 3.9× |
| origin_penalty | 0.7 / 1.0 | 0.0931 / 0.0500 | 1.9× |
| T1 | 12 h / 24 h | 0.0437 / 0.0765 | 1.8× |
| size_effect | 0.5 / 1.0 | 0.0983 / 0.1721 | 1.8× |
| T2 | 10 h / 20 h | 0.0821 / 0.0498 | 1.7× |
| k2 | 0.3 / 9.0 | 0.0541 / 0.0632 | 1.2× |
| B0_jump | 1e7 / 1e8 | 0.0515 / 0.0642 | 1.2× |
| B0_phase2 | 1e8 / 1e9 | 0.0638 / 0.0559 | 1.1× |
| k1 | 3e-12 / 3e-10 | 0.0511 / 0.0617 | 1.2× |
| growth_penalty | 0.3 / 0.9 | 0.0549 / 0.0616 | 1.1× |
| sigma | 0.01 / 0.2 | 0.0612 / 0.0569 | 1.1× |
| m13ko7_copies | 10 / 30 | 0.1289 / 0.0744 | **1.73×** |
| hp17_copies | 5 / 12 | 0.0729 / 0.1074 | **1.47×** |
| moi | 1.0 / 5.0 | 0.0581 / 0.0611 | 1.05× |
| kan_at | 60 / 90 min | 0.0584 / 0.0614 | 1.05× |

**The qualitative conclusion is robust.** Across the whole measured Jumpstart range (3–16), and even at a stronger origin penalty (0.7), the plate ratio is below 1 by round 1–2 and falls every round thereafter.

---

## 8. The mechanistic Jumpstart (`jumpstart_mechanistic.py`) — an alternative to a fixed input

`jumpstart.py`'s `phi` is a measured input, not a derived quantity. `jumpstart_mechanistic.py` instead models phagemid and helper genomes competing to replicate inside a dual-infected cell, toward a shared per-cell capacity (`Ncap`). Phagemid enters already established at `Nres` resident copies (18, the real, confirmed backbone value); helper enters as a single injected genome, replicating slightly slower (`replication_ratio=1.05`) due to its interrupted origin. `phi(t) = alpha·gp(t) / (alpha·gp(t) + gh(t))` — no longer constant, but drifting over the incubation.

Swap it into the full chain with `recursion.py --mechanistic`; it lands within ~2% of the fixed-`phi` results at every stage.

**`alpha` (≈0.21, calibrated) is the model's one honestly unresolved number.** Raw copy number (`Nres=18`) alone predicts a much higher ratio than measured; `alpha<1` is needed to bring it back down — a packaging-affinity discount with no confirmed mechanism. Five candidate explanations, none resolved:

1. **Copy number isn't packageable template.** Most of the ~18 resident copies may be ordinary double-stranded plasmid, never engaged in the rolling-circle synthesis that actually makes packageable single strands. Testable: an ssDNA-vs-total-DNA assay.
2. **`alpha` conflates packaging rate with particle quality.** It's inferred from colony counts, so it also captures how often a packaged particle is defective. Testable: physical (qPCR/nanoparticle tracking) vs. infectious particle counts.
3. **A genuinely weaker packaging signal** on pHK002.2's own f1 intergenic region. Testable: sequence it and compare to the canonical element.
4. **Competition for pV** favoring the helper's own newly-made strands in cis. Weakest of the four; no direct evidence either way.
5. **M13K07's origin defect may be larger than modeled** — the interruption is documented to cause the helper to be *outcompeted* by any co-resident wild-type origin, which could mean `replication_ratio` (currently a "weak" 1.05) is itself an underestimate.

This is genuinely open. It needs wet-lab discrimination this model cannot generate on its own.

---

## 9. Known limitations

1. **Titer collapse not reproduced.** Total particle output stays roughly constant round to round; the data show a ~10× crash in colonies by round 3. Nothing in the model currently removes productivity or capacity over successive rounds. Per the user's assessment, this may partly reflect plate-methodology issues (saturation, an inconsistently-applied lawn threshold in earlier runs) rather than being purely biological — but it is not reproduced either way.
2. **`eps`, `origin_penalty`, and the mechanistic Jumpstart's `alpha` are the model's real unknowns.** `eps` and `jump_ratio` remain the biggest active levers on the exact numbers (§7).
3. **`m13ko7_copies`/`hp17_copies` are general-literature values**, not measurements of these specific constructs, **and the model is not robust to exactly where within their ranges they sit** (§5, §7). The published p15A range for M13K07's copy number is ~10–30/cell — `a_h` moves from ~0.58 at the low end to ~0.73 at the high end, and the round-2 plate ratio by 1.73×. This parameter's history is worth recording precisely: an initial default of **200** was the wrong figure entirely (wild-type M13's rolling-circle number, not M13K07's p15A-driven one); corrected to **10** (the low end of the real range); corrected again to the current **20**, once it was clear the measured round-1 direction favours the high end of that same range. All three rest on one underlying fact — M13K07 replicates via p15A, not the damaged f1 origin — and only the point estimate moved, each time on better evidence. A direct qPCR of M13K07's copy number in this host background would settle it outright.
4. **The plate readout is idealized** — linear, infectivity-weighted, no saturation — while real plates saturate and are overdispersed.
5. **Not modeled at all:** eclipse delay, kanamycin persisting through the overnight (a one-time kill instead), plasmid segregation loss, antibiotic carry-over at passage, within-round stochasticity, HP17_SORT3 (a second, differently-behaving construct).
6. **`Lp`/`Lh` only act through `size_effect`** (0 by default); genome length does nothing at current settings.

---

## 10. Open questions

1. **Which of the five `alpha` hypotheses (§8) is right?** The wet-lab tests that would distinguish them are named there.
2. **What sets `B0_phase2`?** Unlike `B0_jump` (tied to a real A₆₀₀ target), this is still an assumed default. Low priority — §7 shows it barely matters.
3. **Were plating volumes and dilutions identical across rounds?** Needed to properly convert colony counts into titers and test the round-3 collapse (§9, limitation 1).
4. **Is HP17_SORT3 worth modelling** (it behaves differently — levels off near 0.4 rather than continuing to fall), or should it stay a comparison only?

---

## 11. Experimental data

Two constructs: **N011/HP17K07** (matches this model) and **N010/HP17_SORT3** (not modeled). Two independent runs gave Jumpstart ratios of ~5.08 and ~1.53 (by-row range 0.61–16.5); Phase 2 round 1 measured 0.004–0.70 across runs, all below 1, matching the model's direction. Round 3 shows a real colony crash in the data (e.g. N011: Amp 142→66→0, Kan 27→107→18) that the model does not reproduce (§9, limitation 1). Known data-quality caveats: plate saturation at high titers, a lawn threshold now fixed at >150 colonies, and substantial run-to-run variability in round 1 that may be lab-handling noise rather than a modeling target.

---

## 12. Testing

`test_model.py` — 52 automated checks (~80 s full, ~8 s `--fast`), covering: mass-balance and neutrality invariants (interchangeable lineages must not move the ratio; an absent species can never appear), the equations themselves (exclusion applies identically to every particle type; the hybrid-coat gene-dosage derivation), the passage step's Poisson statistics, the mechanistic Jumpstart's identifiability collapse, band-mode extinction handling, and exact regression values pinned to the current code.

---

## 13. Quick reference

```bash
python3 recursion.py                          # the whole assay, defaults
python3 recursion.py -r 10 --plot chain.png --csv chain.csv
python3 recursion.py --jump-ratio 3           # low end of the measured Jumpstart range
python3 recursion.py --band                   # report a range (3-16) instead of a single line
python3 recursion.py --mechanistic            # use the mechanistic Jumpstart instead of fixed phi
python3 recursion.py --band --mechanistic     # both together
python3 recursion.py --origin-penalty 0.7     # a stronger pullback than the default 0.9
python3 recursion.py --size-effect 1          # include the phagemid production advantage
python3 jumpstart.py                          # Jumpstart alone, measured 4.88:1 split
python3 jumpstart_mechanistic.py              # mechanistic Jumpstart alone, real Nres, calibrated alpha
python3 jumpstart_mechanistic.py --sensitivity --plot sens.png   # the two falsifiable predictions
python3 phase2.py -e 0.05 --origin-penalty 0.9   # one Phase 2 round
python3 passage.py                            # passage demonstration
python3 test_model.py                         # 52 checks
```
