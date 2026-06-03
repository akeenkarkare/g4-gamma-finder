# Reference notes from the paper's repo (RyotaroOKabe/radiation_mapping)

Cloned to /tmp/radiation_mapping_ref (depth 1). Key findings for our recreation.

## Data format (what the NN consumes)
- **Input x**: flattened detector readout (normalized).
  - Square 2x2 -> 4 values.
  - Tetromino -> 2x3 = 6 values, with the 2 vacant cells filled with 0.
- **Target y**: angular distribution over `seg_angles = 64` sectors covering [0, 360).
  Each sector = fraction of incident intensity from that direction (pie-chart rep).
- **Sampling** (gen_data_square.py): per data point, random
  - distance in [20, 500] cm
  - angle in [0, 360) deg (float)
  - energy 0.5 MeV (5e5 eV)
  - num_particles = 50000 per MC point; num_data = 3000.
- Normalization: square readout normalized to mean 0 / std 1 (per the paper).

## Model (utils/model.py)
- `MyNet2` = `Filterlayer2` (2 channels: far-field + near-field global filters,
  initialized from high-quality sims) -> `UNet(c1=32)` -> squeeze -> softmax.
- Two-filter-layer version is the main model; `MyNet1` is the single-filter baseline.
- Filters are (2, 64, h, w); fine-tuned with a LOWER LR than the rest.

## Loss (utils/emd_ring_torch.py)
- `emd_loss_ring(p, p_hat, r=2)`: CYCLIC Wasserstein / earth-mover distance.
  - Normalizes both to PMFs.
  - For each of n cyclic shifts, compute ECDF diff, take r-norm, then MIN over shifts.
  - r=2 for training (converges faster); r=1 for evaluation (= angle error in deg).
- Cost matrix is cyclic: `M[i,j] = min(|i-j|, j+n-i, i+n-j)`.

## Training (utils/model.py Model.train)
- Adam, lr 1e-3 for most params, 3e-5 for filter layer.
- 90/10 train/test split, 5-fold CV.
- Metric for ranking configs: 1st Wasserstein distance (avg angular error).

## How this maps to OUR Geant4 setup
- Our CSV rows are PER-EVENT; paper's x is PER-CONFIG aggregated readout.
  -> Post-process: group events by (angle bucket / config), sum per-pixel energy,
     normalize -> one training sample per source configuration.
- Our target y must be built as the 64-sector one-hot-ish distribution from the
  true source angle (eq. 1 in the paper: triangular interpolation into sectors).
- Tetromino shapes -> output a 2x3 readout (6 cols, 2 always zero).
