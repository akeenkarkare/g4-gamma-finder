# GPU sweeps on the RTX 4070 Windows laptop

Two heavy sweeps belong on the GPU machine. Both use the **MC-initialized U-Net**
(`train_unet.py`, auto-CUDA) -- the validated model -- not the old MLP. The driver
`sweep_unet.py` does everything per shape: generate dataset -> generate MC filters
-> train U-Net -> rank.

## Why these two
1. **Full 35-hexomino re-rank under the U-Net.** The MLP mis-ranks compact shapes
   (it rated h13 mid-pack; the U-Net found h13 = best of everything). So the
   MLP-based "best-6" is unreliable and must be redone with the U-Net.
2. **Footprint-controlled N=7.** A sprawling heptomino re-introduces the footprint
   confound (bigger array blurs angle). We restrict to heptominoes fitting a 3x3
   box (`--maxspan 2`, 7 shapes) so N=7 is comparable to pP(5) and h13(6).

## Prereqs
1. Pull the repo; build the Geant4 app -> `build\exampleB1.exe`.
2. `pip install torch numpy` (CUDA wheel from pytorch.org). Verify GPU:
   `python -c "import torch; print(torch.cuda.is_available())"` -> True.
3. `copy gen_config_step.mac build\`  (the sweep needs it in the build dir).
4. Quick smoke test (tiny, ~2 min): confirms the chain runs before the big jobs:
   `set SEEDS=1 & set EPOCHS=10 & python sweep_unet.py --exe build\exampleB1.exe --mode heptominoes --configs 50 --events 1000`

## Sweep A -- all 35 hexominoes (true best-6)
```
python sweep_unet.py --exe build\exampleB1.exe --mode hexominoes --configs 1000 --events 5000
```
Prints a ranking; top line = best 6-crystal shape under the U-Net.
(Compare to the lucky data point we already have: h13 = 2.58 deg.)

## Sweep B -- footprint-controlled N=7 (the 7 compact heptominoes)
```
python sweep_unet.py --exe build\exampleB1.exe --mode heptominoes --maxspan 2 --configs 1000 --events 5000
```
7 shapes, each fitting a 3x3 footprint. Top line = best compact 7-crystal shape.

## Tuning cost (env vars)
- `SEEDS` (default 3), `EPOCHS` (default 200) for train_unet.
- `CPU=1` forces CPU if needed.
- 5 seeds / 300 epochs for publication-grade precision.

## Send back
For each sweep: the best shape's tag + its cells (printed in the run) and the
mean +/- std (deg), so we can complete the curve:
```
  best-4  S            5.05 deg  (U-Net)
  best-5  pP           4.36 deg  (U-Net, exhaustive over 12 pentominoes)
  best-6  ???          (Sweep A)   -- h13 = 2.58 is the current lower bound
  best-7  ??? compact  (Sweep B)
```
Key question: does best-7 (compact) drop below h13's 2.58 deg? If yes, the
"density within a fixed footprint" trend extends to N=7 -- the headline result.
