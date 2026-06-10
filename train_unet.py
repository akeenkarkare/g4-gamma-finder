#!/usr/bin/env python3
"""
Faithful re-implementation of the paper's directional model (Okabe et al.):
an MC-initialized filter layer + a 1D (cyclic) U-Net over the 64 angle sectors,
trained with the cyclic-Wasserstein loss. Used to check that the crystal-count
KNEE (best-5 P-pentomino < best-6 < best-4) is not an artifact of the simpler
MLP in train_compare.py.

Filter layer: for each field channel (far, near) and each of `seg` sectors, the
MC filter template is the expected normalized per-crystal readout for a source at
that sector (from gen_filters.py). The layer projects the detector readout onto
these templates -> a (2, seg) directional response -> U-Net -> softmax over seg.

Usage (per config):
  python3 train_unet.py --data build/dsShape_pP_nt_pixels.csv \
                        --filters build/filters_pP.npy --label pP
Compares to train_compare.py's MLP numbers on the same data.
"""
import argparse
import os
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# Reuse the validated data loader + cyclic-Wasserstein loss + target builder.
from train_compare import (read_csv, angular_target, normalize_inputs,
                           emd_loss_ring, SEG, SECTOR_DEG, DEVICE)

# Lean defaults: 3 seeds / 200 epochs is plenty to confirm a RANKING (the knee).
# Bump back to 5 seeds / 300+ epochs for publication-grade precision (e.g. on GPU).
EPOCHS = int(os.environ.get("EPOCHS", "200"))
BATCH = 64
SEEDS = list(range(int(os.environ.get("SEEDS", "3"))))


# ---------------------------------------------------------------------------
# Model: MC-initialized filter layer + 1D cyclic U-Net  (paper's MyNet2)
# ---------------------------------------------------------------------------
class FilterLayer(nn.Module):
    """Project the n-crystal readout onto MC angle templates for 2 channels.
    filt: numpy (2, seg, n_in). Output: (B, 2, seg)."""
    def __init__(self, filt):
        super().__init__()
        seg, n_in = filt.shape[1], filt.shape[2]
        # weight[c] is (n_in, seg): out = x @ weight  -> (B, seg).
        w0 = torch.tensor(filt[0].T.copy(), dtype=torch.double)
        w1 = torch.tensor(filt[1].T.copy(), dtype=torch.double)
        self.weight0 = nn.Parameter(w0)
        self.weight1 = nn.Parameter(w1)
        self.b0 = nn.Parameter(torch.zeros(1, seg, dtype=torch.double))
        self.b1 = nn.Parameter(torch.zeros(1, seg, dtype=torch.double))
        self.wn0 = nn.Parameter(torch.ones(1, dtype=torch.double))
        self.wn1 = nn.Parameter(torch.ones(1, dtype=torch.double))

    def forward(self, x):
        o0 = torch.matmul(x, self.weight0) / self.wn0 + self.b0
        o1 = torch.matmul(x, self.weight1) / self.wn1 + self.b1
        return torch.stack([o0, o1], dim=1)          # (B, 2, seg)


class DoubleConv(nn.Module):
    def __init__(self, ci, co):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(ci, co, 3, padding=1, padding_mode='circular'), nn.ELU(),
            nn.Conv1d(co, co, 3, padding=1, padding_mode='circular'), nn.ELU())

    def forward(self, x): return self.net(x)


class UNet1D(nn.Module):
    """Small cyclic 1D U-Net over the seg-axis (2 down/up levels keep seg=64
    divisible). Circular padding everywhere because angle is periodic."""
    def __init__(self, c1=32):
        super().__init__()
        self.d1 = DoubleConv(2, c1)
        self.d2 = DoubleConv(c1, c1 * 2)
        self.bott = DoubleConv(c1 * 2, c1 * 4)
        self.pool = nn.MaxPool1d(2)
        self.up2 = nn.ConvTranspose1d(c1 * 4, c1 * 2, 2, stride=2)
        self.c2 = DoubleConv(c1 * 4, c1 * 2)
        self.up1 = nn.ConvTranspose1d(c1 * 2, c1, 2, stride=2)
        self.c1c = DoubleConv(c1 * 2, c1)
        self.out = nn.Conv1d(c1, 1, 3, padding=1, padding_mode='circular')

    def forward(self, x):
        x1 = self.d1(x)
        x2 = self.d2(self.pool(x1))
        xb = self.bott(self.pool(x2))
        u2 = self.c2(torch.cat([self.up2(xb), x2], dim=1))
        u1 = self.c1c(torch.cat([self.up1(u2), x1], dim=1))
        return self.out(u1).squeeze(1)               # (B, seg)


class PaperNet(nn.Module):
    def __init__(self, filt):
        super().__init__()
        self.fl = FilterLayer(filt)
        self.unet = UNet1D(c1=32)

    def forward(self, x):
        return torch.softmax(self.unet(self.fl(x)), dim=1)


# ---------------------------------------------------------------------------
def train_eval_seed(Xn, Y, filt, seed):
    torch.manual_seed(seed)
    rng = np.random.RandomState(seed)
    n = len(Xn); idx = rng.permutation(n); ntest = max(1, n // 10)
    te, tr = idx[:ntest], idx[ntest:]
    dt = torch.float64
    Xtr = torch.tensor(Xn[tr], dtype=dt, device=DEVICE)
    Ytr = torch.tensor(Y[tr], dtype=dt, device=DEVICE)
    Xte = torch.tensor(Xn[te], dtype=dt, device=DEVICE)
    Yte = torch.tensor(Y[te], dtype=dt, device=DEVICE)

    net = PaperNet(filt).double().to(DEVICE)
    # Filter-layer params at a lower LR (fine-tune the MC init), rest normal.
    fl_params = list(net.fl.parameters())
    other = [p for p in net.parameters() if all(p is not q for q in fl_params)]
    opt = torch.optim.Adam([{"params": other, "lr": 1e-3},
                            {"params": fl_params, "lr": 3e-5}])
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS)

    best = float("inf")
    for ep in range(EPOCHS):
        net.train()
        perm = torch.randperm(len(Xtr), device=DEVICE)
        for b in range(0, len(Xtr), BATCH):
            bi = perm[b:b + BATCH]
            opt.zero_grad()
            loss = emd_loss_ring(Ytr[bi], net(Xtr[bi]), r=2)
            loss.backward(); opt.step()
        sched.step()
        net.eval()
        with torch.no_grad():
            best = min(best, emd_loss_ring(Yte, net(Xte), r=1).item())
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--filters", required=True)
    ap.add_argument("--label", default="config")
    args = ap.parse_args()

    X, ang = read_csv(args.data)
    Xn = normalize_inputs(X)
    Y = angular_target(ang)
    filt = np.load(args.filters)
    assert filt.shape[2] == X.shape[1], \
        f"filter crystals {filt.shape[2]} != data crystals {X.shape[1]}"

    print(f"[{args.label}] {len(X)} configs, {X.shape[1]} crystals, "
          f"MC-init U-Net, {len(SEEDS)} seeds (device={DEVICE})", flush=True)
    per = []
    for s in SEEDS:
        deg = train_eval_seed(Xn, Y, filt, s) * SECTOR_DEG
        per.append(deg)
        print(f"  seed {s}: best test error = {deg:.2f} deg", flush=True)
    a = np.array(per)
    print(f"  -> {args.label}: mean = {a.mean():.2f} deg (std {a.std():.2f}, "
          f"min {a.min():.2f}, max {a.max():.2f})")


if __name__ == "__main__":
    main()
