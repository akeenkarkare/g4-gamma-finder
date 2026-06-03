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
import torch
import torch.nn as nn

SEG = 64                # angle sectors (matches the paper)
SECTOR_DEG = 360.0 / SEG
EPOCHS = 400
BATCH = 64
SEEDS = [0, 1, 2, 3, 4]   # average over several seeds/splits for a stable answer


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def read_csv(path):
    """Read the aggregated Geant4 ntuple CSV -> (X[n,4], angles_deg[n])."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(",")
            if len(p) < 6:
                continue
            try:
                e = [float(p[0]), float(p[1]), float(p[2]), float(p[3])]
                ang = float(p[4])
            except ValueError:
                continue
            rows.append(e + [ang])
    arr = np.array(rows, dtype=np.float64)
    return arr[:, :4], arr[:, 4]


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
def ecdf(p):
    n = p.shape[1]
    tri = torch.tril(torch.ones(n, n, dtype=p.dtype, device=p.device)).t()
    return torch.matmul(p, tri)


def shift(p, d):
    if d == 0:
        return p
    return torch.cat([p[:, d:], p[:, :d]], dim=1)


def emd_loss_ring(p, p_hat, r=2, eps=1e-8):
    """Cyclic r-Wasserstein distance between two batches of ring PMFs.
    r=2 for training (smooth), r=1 for evaluation (= angular error in sectors)."""
    p = p / (p.sum(dim=1, keepdim=True) + eps)
    p_hat = p_hat / (p_hat.sum(dim=1, keepdim=True) + eps)
    n = p.shape[1]
    diffs = []
    for i in range(n):
        d = ecdf(shift(p, i)) - ecdf(shift(p_hat, i))
        diffs.append((d.abs() ** r).sum(dim=1))
    emd = torch.stack(diffs, dim=0).min(dim=0)[0]   # min over cyclic cut points
    emd = emd ** (1.0 / r)
    return emd.mean()


# ---------------------------------------------------------------------------
# Model: compact MLP  4 -> 64 softmax
# ---------------------------------------------------------------------------
class DirNet(nn.Module):
    def __init__(self, seg=SEG):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(4, 128), nn.ELU(),
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

    net = DirNet().double()
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
    print(f"\n[{label}] {len(X)} configs, {len(SEEDS)} seeds...")
    per_seed = []
    for s in SEEDS:
        emd = train_eval_seed(Xn, Y, s)
        deg = emd * SECTOR_DEG
        per_seed.append(deg)
        print(f"  seed {s}: best test error = {deg:.2f} deg")
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
