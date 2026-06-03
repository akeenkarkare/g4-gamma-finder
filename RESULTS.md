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
```
