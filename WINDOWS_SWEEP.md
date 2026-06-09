# Running the 7-crystal (heptomino) sweep on the RTX 4070 Windows laptop

The 108 heptominoes are too many to run on the Mac, so generate them here on the
GPU machine. Geometry uses the runtime `/det/cells` command (no recompile needed
per shape), and `train_compare.py` auto-uses CUDA if available.

## Prereqs
1. Pull the repo.
2. Build the Geant4 app (Windows): produces `build\exampleB1.exe`.
   (Standard Geant4 + CMake build; same CMakeLists as macOS.)
3. Python env with PyTorch (CUDA build) + numpy:
   `pip install torch numpy`  (use the CUDA wheel from pytorch.org)
4. Make sure `gen_config_step.mac` is in the `build` dir
   (`copy gen_config_step.mac build\`).

## Step 1 - generate all 108 heptomino datasets (Geant4, CPU-bound)
```
python heptomino_sweep.py --exe build\exampleB1.exe --configs 1000 --events 5000
```
This writes `build\dsHep_hep000_nt_pixels.csv` ... `hep107`. It prints each shape's
cell spec as it goes. ~108 x ~30 s ≈ under an hour on a fast CPU.
(Tip: test first with `--limit 3`.)

## Step 2 - train + rank all 108 (PyTorch, GPU)
```
python train_compare.py build\dsHep_*_nt_pixels.csv
```
With CUDA this is much faster than the Mac. The ranking's top line is the best
7-crystal shape; note its `hepNNN` tag and look up its cells in the step-1 output.

## Step 3 - send back the result
Report the best-7 mean +/- std (deg) and which hepNNN it was, so it can be added
to the crystal-count knee curve alongside best-4 (S, ~7-9), best-5 (P-pentomino,
~4.7), best-6 (from the Mac hexomino sweep).

## Force CPU (if no GPU): set env CPU=1 before train_compare.py.
## Knee context so far (idealized, mono 0.5 MeV, 1000 configs):
##   best-4 S         ~9.0 deg
##   best-5 P-pentomino ~4.7 deg   <- current optimum
##   best-6           (Mac sweep, pending)
##   best-7           (this sweep)
