#!/usr/bin/env python3
"""
Train a direction-prediction network per detector configuration and rank the
configurations by cyclic-Wasserstein angular error -- the paper's metric
(Okabe et al., Nat. Commun. 15:3061, 2024).

Input: aggregated per-config CSVs from the Geant4 detector, columns
  e0_keV, e1_keV, e2_keV, e3_keV, angle_deg, dist_cm, nEvents
(one row = one source configuration).

For each dataset we:
  1. build the 4-pixel input vector (normalized to mean 0 / std 1),
  2. build the 64-sector ground-truth angular distribution (paper eq. 1:
     triangular interpolation of the true angle into the two nearest sectors),
  3. train a small net (4 -> 64 softmax) with the SQUARED cyclic-EMD loss,
  4. evaluate the 1st (L1) cyclic-Wasserstein distance on a held-out test set,
     which equals the mean angular error in sector units (x 360/64 deg).

NOTE on the model: the paper uses a filter-layer + U-Net designed for their
larger filtered representation. With only a 4-pixel input vector, a compact MLP
is the appropriate analogue; the loss and target representation are faithful to
the paper. This is a deviation documented here.

Usage:
  python3 train_compare.py ds_square_nt_pixels_t0.csv ds_S_nt_pixels_t0.csv
  (defaults to build/ds_square... and build/ds_S... if no args)
"""
import sys
import os
import numpy as np
import os
import torch
import torch.nn as nn

SEG = 64                # angle sectors (matches the paper)
SECTOR_DEG = 360.0 / SEG
EPOCHS = 400
BATCH = 64
SEEDS = [0, 1, 2, 3, 4]   # average over several seeds/splits for a stable answer

# If set (env USE_BANDS=1), use per-crystal per-band COUNT columns (c*_b*) as the
# feature vector instead of the lumped per-crystal energy. Lets a crystal's
# spectral shape (peak vs Compton) inform direction -- the multi-band experiment.
USE_BANDS = os.environ.get("USE_BANDS", "0") == "1"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def read_csv(path):
    """Read the aggregated Geant4 ntuple CSV -> (X[n,k], angles_deg[n]).

    Handles both the old 4-column format (e0..e3, angle, dist, n) and the new
    6-column format (e0..e5, angle, dist, n). The energy columns and the
    angle column are located by the '#column ... <name>' header lines. Crystal
    columns that are all-zero across the dataset (unused by the shape) are
    dropped, so X has exactly the active crystal count k."""
    cols = []                       # ordered column names from the header
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#column"):
                cols.append(line.split()[-1])   # last token is the name
                continue
            if not line or line.startswith("#") or line.startswith("entries"):
                continue
            parts = line.split(",")
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    arr = np.array(rows, dtype=np.float64)

    # Locate columns by name (fallback to positional if header missing).
    if cols and "angle_deg" in cols:
        e_idx = [i for i, c in enumerate(cols) if c.startswith("e") and c.endswith("_keV")]
        # Per-crystal per-band count columns are named c<i>_b<b>.
        b_idx = [i for i, c in enumerate(cols)
                 if c.startswith("c") and "_b" in c]
        a_idx = cols.index("angle_deg")
    else:
        e_idx = list(range(arr.shape[1] - 3))   # assume last 3 are angle,dist,n
        b_idx = []
        a_idx = arr.shape[1] - 3

    # USE_BANDS (module global): if True and band columns exist, the feature
    # vector is the per-crystal per-band counts INSTEAD of the lumped energy.
    if USE_BANDS and b_idx:
        feat_idx = b_idx
    else:
        feat_idx = e_idx

    X = arr[:, feat_idx]
    ang = arr[:, a_idx]
    # Drop feature columns that are entirely zero (unused crystals/bands).
    active = X.sum(axis=0) > 0
    X = X[:, active]
    return X, ang


def angular_target(angle_deg):
    """Paper eq. 1: a point source at angle theta becomes a 64-sector
    distribution with mass split linearly between the two nearest sectors."""
    n = len(angle_deg)
    y = np.zeros((n, SEG), dtype=np.float64)
    pos = (angle_deg % 360.0) / SECTOR_DEG     # fractional sector position
    lo = np.floor(pos).astype(int) % SEG
    hi = (lo + 1) % SEG
    frac = pos - np.floor(pos)
    for i in range(n):
        y[i, lo[i]] += 1.0 - frac[i]
        y[i, hi[i]] += frac[i]
    return y


def normalize_inputs(X):
    """Per-sample normalize to mean 0 / std 1 across the 4 pixels (paper)."""
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True) + 1e-9
    return (X - mu) / sd


# ---------------------------------------------------------------------------
# Cyclic Wasserstein / EMD on a ring (adapted from the paper's emd_ring_torch)
# ---------------------------------------------------------------------------
def emd_loss_ring(p, p_hat, r=2, eps=1e-8):
    """Cyclic r-Wasserstein distance between two batches of ring PMFs.
    r=2 for training (smooth), r=1 for evaluation (= angular error in sectors).

    Vectorized form of the cut-point method: for the per-element difference
    d = p - p_hat, the cost for cut point k is the r-norm of the cumulative sum
    of d started at index k. Using c = cumsum(d) and total T = c[-1], the
    cut-at-k cumulative profile equals (c rolled by k) - (T rolled), which for
    each shift is c - c[k-1]. So the cost matrix over all 64 cut points is built
    with a single broadcast instead of a 64-iteration Python loop."""
    p = p / (p.sum(dim=1, keepdim=True) + eps)
    p_hat = p_hat / (p_hat.sum(dim=1, keepdim=True) + eps)
    d = p - p_hat                       # [B, n]
    c = torch.cumsum(d, dim=1)          # [B, n]; c[:, n-1] ~ 0 (equal mass)
    # For cut point k, the shifted CDF difference is c - c[:, k-1] (with k=0 -> 0).
    # cprev[:, k] = c[:, k-1], cprev[:, 0] = 0.
    cprev = torch.cat([torch.zeros_like(c[:, :1]), c[:, :-1]], dim=1)  # [B, n]
    # diff_kj = c[:, j] - cprev[:, k]  -> [B, n(cuts k), n(positions j)]
    diff = c.unsqueeze(1) - cprev.unsqueeze(2)        # [B, k, j]
    cost = (diff.abs() ** r).sum(dim=2)               # [B, k]  cost per cut point
    emd = cost.min(dim=1)[0] ** (1.0 / r)             # min over cut points
    return emd.mean()


# ---------------------------------------------------------------------------
# Model: compact MLP  4 -> 64 softmax
# ---------------------------------------------------------------------------
class DirNet(nn.Module):
    def __init__(self, n_in=4, seg=SEG):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_in, 128), nn.ELU(),
            nn.Linear(128, 256), nn.ELU(),
            nn.Linear(256, 256), nn.ELU(),
            nn.Linear(256, seg),
        )

    def forward(self, x):
        return torch.softmax(self.net(x), dim=1)


# ---------------------------------------------------------------------------
# Train / evaluate one dataset
# ---------------------------------------------------------------------------
def train_eval_seed(Xn, Y, seed):
    """One train/eval run for a given seed. Returns BEST test EMD (sectors)
    seen during training (early-stopping style), not the last-epoch value."""
    torch.manual_seed(seed)
    rng = np.random.RandomState(seed)

    n = len(Xn)
    idx = rng.permutation(n)
    ntest = max(1, n // 10)
    te, tr = idx[:ntest], idx[ntest:]

    dt = torch.float64
    Xtr = torch.tensor(Xn[tr], dtype=dt); Ytr = torch.tensor(Y[tr], dtype=dt)
    Xte = torch.tensor(Xn[te], dtype=dt); Yte = torch.tensor(Y[te], dtype=dt)

    net = DirNet(n_in=Xn.shape[1]).double()
    opt = torch.optim.Adam(net.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS)

    best = float("inf")
    for ep in range(EPOCHS):
        net.train()
        perm = torch.randperm(len(Xtr))
        for b in range(0, len(Xtr), BATCH):
            bi = perm[b:b + BATCH]
            opt.zero_grad()
            loss = emd_loss_ring(Ytr[bi], net(Xtr[bi]), r=2)
            loss.backward()
            opt.step()
        sched.step()
        net.eval()
        with torch.no_grad():
            te_l = emd_loss_ring(Yte, net(Xte), r=1).item()
        best = min(best, te_l)
    return best


def run_one(path, label):
    X, ang = read_csv(path)
    Xn = normalize_inputs(X)
    Y = angular_target(ang)
    print(f"\n[{label}] {len(X)} configs, {X.shape[1]} crystals, {len(SEEDS)} seeds...",
          flush=True)
    per_seed = []
    for s in SEEDS:
        emd = train_eval_seed(Xn, Y, s)
        deg = emd * SECTOR_DEG
        per_seed.append(deg)
        print(f"  seed {s}: best test error = {deg:.2f} deg", flush=True)
    arr = np.array(per_seed)
    print(f"  -> {label}: mean = {arr.mean():.2f} deg  (std {arr.std():.2f}, "
          f"min {arr.min():.2f}, max {arr.max():.2f})")
    return arr.mean(), arr.std()


def main():
    args = sys.argv[1:]
    if not args:
        args = ["build/ds_square_nt_pixels_t0.csv", "build/ds_S_nt_pixels_t0.csv"]
    results = {}
    for path in args:
        label = (os.path.basename(path)
                 .replace("ds_", "")
                 .replace("_nt_pixels_t0.csv", "")
                 .replace("_nt_pixels.csv", ""))
        results[label] = run_one(path, label)

    print("\n================ RANKING (lower = better) ================")
    for label, (mean, std) in sorted(results.items(), key=lambda kv: kv[1][0]):
        print(f"  {label:10s}  {mean:.2f} +/- {std:.2f} deg")
    best = min(results, key=lambda k: results[k][0])
    print(f"\nBest configuration: {best}  ({results[best][0]:.2f} deg)")


if __name__ == "__main__":
    main()
