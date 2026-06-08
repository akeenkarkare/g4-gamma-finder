# Results

Recreation of Okabe et al., *Tetris-inspired detector with neural network for
radiation mapping*, Nat. Commun. 15:3061 (2024), using our own hardware model:
**2-inch cylindrical LaBr3 crystals in 2 mm aluminium casings**, 4-crystal array,
1 mm lead inter-pixel padding, 0.5 MeV gammas, source coplanar with the detector.

Pipeline: Geant4 (per-pixel scoring, configurable shape/padding, random
(distance, angle) source sampling) -> aggregated per-config CSV -> a directional
net trained with the paper's **cyclic Wasserstein (earth-mover) loss**, ranked by
1st-Wasserstein angular error. Robust evaluation: best-epoch tracking, cosine LR
decay, mean over 5 seeds.

## Shape sweep (1000 configs/shape, 5000 events/config)

Lower angular error = better directional resolution.

| Rank | Shape  | Angular error (deg) |
|------|--------|---------------------|
| 1    | **S**  | **7.38 ± 1.41**     |
| 2    | square | 14.34 ± 2.53        |
| 3    | J      | 15.48 ± 1.91        |
| 4    | T      | 18.17 ± 2.80        |
| 5    | L      | 19.75 ± 1.96        |

### Interpretation
- **S-shape wins decisively** (~2x better than any other; its error band does not
  overlap any other shape's). This reproduces the paper's headline finding that
  the S-shape Tetromino gives the best directional resolution.
- The ranking is essentially a **symmetry ranking**. The square's 4-fold symmetry
  creates a mirror ambiguity *within* each quadrant (it resolves which quadrant a
  source is in, but struggles to distinguish e.g. 30 deg from 60 deg). The S-shape
  has no such symmetry, so every pixel sees a distinct signal at every angle.
- T ranks near the bottom, consistent with the paper: the T-shape's left/right
  pixels receive identical signals from a front source, reducing effective pixels.
- Caveat: square vs J is within ~1 deg with overlapping bands -- not separable.
  S vs everything else is robust.

**Design recommendation: build the S-shape arrangement.**

## Padding sweep (S-shape, lead thickness)

Holding shape = S, sweeping the lead inter-pixel padding thickness
(1000 configs each, same evaluation protocol).

| Padding (mm) | Angular error (deg) |
|--------------|---------------------|
| **1.0**      | **7.38 ± 1.41**     |
| 0.5          | 8.03 ± 1.66         |
| 2.0          | 8.05 ± 0.60         |
| 3.0          | 10.33 ± 1.86        |
| 5.0          | 10.95 ± 1.29        |

### Interpretation
- **1 mm padding (the paper's value) is optimal.** 1.0 and 2.0 mm are
  statistically tied (overlapping bands); both clearly beat thicker lead.
- Beyond ~2 mm, resolution degrades: thicker lead absorbs too many gammas
  (less signal) without adding useful inter-pixel contrast.
- No tuning gain available here -- the original configuration was already
  well-chosen.

**Best 4-crystal design: S-shape, 1 mm lead padding -> ~7-9 deg accuracy.**
(Superseded overall by the 5-crystal P-pentomino -- see the crystal-count
section below.)

## Number of crystals -- fair shape search (4 vs 5 vs 6)

A first attempt compared the optimized 4-shape (S) against ARBITRARY 5- and
6-shapes and wrongly concluded "more is not better". The fair test enumerates
candidate shapes for EVERY count and lets the data pick the best of each:
- 4 crystals: square, S, J, T, L (S won earlier at 9.0 deg)
- 5 crystals: ALL 12 free pentominoes (F,I,L,N,P,T,U,V,W,X,Y,Z)
- 6 crystals: 8 hexominoes (asymmetric-weighted + a 2x3 symmetric baseline)
Idealized conditions, 1000 configs each, 5 seeds, same protocol.

### Best shape per crystal count
| Crystals | Best shape        | Angular error (deg) |
|----------|-------------------|---------------------|
| **5**    | **P-pentomino**   | **4.69 ± 0.59**     |
| 6        | extended-N hexomino | 6.83 ± 0.53       |
| 4        | S-tetromino       | 9.00 ± 1.45         |

### Interpretation
- **The earlier conclusion was wrong.** With a fair shape search, the best
  5-crystal (P-pentomino, 4.7 deg) is ~2x better than the best 4-crystal (S,
  9.0 deg), and the best 6-crystal also beats the best 4. The original "more is
  not better" came from comparing optimized-4 vs arbitrary-5/6.
- **Shape dominates count.** Within 5 crystals the error ranges 4.7 deg (P) to
  47 deg (I, a straight line that is nearly 1-D and cannot resolve front/back) --
  a 10x spread. You cannot judge a crystal count without first finding its best
  shape.
- **Knee is at 5.** best-4 (9.0) -> best-5 (4.7) is a large real gain, but
  best-5 -> best-6 REGRESSES (4.7 -> 6.8): the 6th crystal enlarges the footprint
  and dilutes per-pixel contrast more than it adds information. The P-pentomino
  is the sweet spot.
- The P-pentomino is a 2x2 block plus one extra crystal on the side -- compact
  and asymmetric.
- Caveat: all 12 pentominoes were tested (exhaustive for 5), but only 8 of the
  35 hexominoes. A better hexomino may exist, but every tested hexomino lost to
  the P-pentomino.

**Revised optimal design: 5-crystal P-pentomino, 1 mm lead padding -> ~4.7 deg.**

## Realistic field performance (S-shape)

The numbers above use an idealized single-energy (0.5 MeV) source with no
background. To estimate real field performance we retrained on a dataset that
adds the actual measurement conditions:
- external **Eu-152 + Ba-133** source (discrete lines, weighted by intensity),
- **20% of events as 138La internal background** (the intrinsic LaBr3 activity),
- per-pixel readout taken **only within a [40, 450] keV photopeak ROI** (the
  energy cut you would apply to real data).

| Conditions | Angular error (deg) |
|------------|---------------------|
| Idealized (0.5 MeV mono, no background, full energy) | 7.38 ± 1.41 |
| **Realistic (Eu/Ba + 20% 138La background + ROI)**   | **8.71 ± 1.09** |

### Interpretation
- Realistic conditions cost only **~1.3 deg** (~18% relative); the error bands
  overlap, so the degradation is mild.
- It holds up because (1) the ROI keeps the low-energy directional source and
  rejects most 138La, and (2) the 138La that leaks through is symmetric across
  pixels, adding noise but not directional bias.
- **Expect ~8-9 deg directional accuracy in the field** with the S-shape, the
  Eu/Ba source, and a standard photopeak ROI cut. LaBr3 self-activity does not
  meaningfully degrade this provided the energy window is applied.

### Background discrimination (energy + symmetry)
`analyze_discrimination.py` quantifies the two separation levers on a mixed
measurement: the source puts ~79% of its counts in the low-energy ROI and is
strongly directional (one pixel ~0.94 of the in-ROI signal), while 138La is
~70% outside the ROI and uniform across pixels (~0.25 each) -- so it is
separable from the source by BOTH energy and symmetry.

## Reproduce

```bash
# 1. Generate datasets (Geant4), ~110 s per shape
cd build
for s in square S J T L; do ./exampleB1 gen_$s.mac; done

# 2. Train + rank (PyTorch venv)
cd ..
source ml_env/bin/activate
python3 train_compare.py \
  build/ds_square_nt_pixels_t0.csv build/ds_S_nt_pixels_t0.csv \
  build/ds_J_nt_pixels_t0.csv build/ds_T_nt_pixels_t0.csv \
  build/ds_L_nt_pixels_t0.csv

# 3. Realistic S-shape run (Eu/Ba source + 138La background + ROI)
cd build && ./exampleB1 gen_realistic_S.mac && cd ..
python3 train_compare.py build/dsRealS_nt_pixels.csv

# 4. Background-discrimination demo (run the two macros first to make spectra)
python3 analyze_discrimination.py build/mix_src build/mix_int
```
