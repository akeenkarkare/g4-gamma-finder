#!/usr/bin/env python3
"""
Generate the MC directional filter templates for a detector shape, as used by
the paper's filter-layer + U-Net (Okabe et al.).

For each of `seg` angle sectors and each field channel (far, near), we simulate a
clean point source at that angle/distance and record the normalized per-crystal
readout. The stack is the (2, seg, n_crystals) filter tensor the network's first
layer is initialized with.

Usage (from project root, exe already built):
  python3 gen_filters.py --exe build/exampleB1 --shape pP --seg 64 \
      --far 400 --near 40 --events 20000 --out build/filters_pP.npy

Produces a .npy of shape (2, seg, n_crystals): channel 0 = far, 1 = near.
"""
import argparse
import os
import subprocess
import sys
import numpy as np


def read_aggregated(csv_path):
    """Return (energy_vector, n_crystals) from a single-config aggregated CSV."""
    cols, row = [], None
    for line in open(csv_path):
        line = line.rstrip("\n")
        if line.startswith("#column"):
            cols.append(line.split()[-1])
        elif line and not line.startswith("#") and not line.startswith("entries"):
            row = [float(x) for x in line.split(",")]
    e_idx = [i for i, c in enumerate(cols) if c.startswith("e") and c.endswith("_keV")]
    e = np.array([row[i] for i in e_idx])
    return e


def sim_one(exe, builddir, shape, angle, dist, events, tag):
    """Run one config (fixed angle/distance) and return the per-crystal vector."""
    mac = os.path.join(builddir, f"filt_{tag}.mac")
    base = f"filt_{tag}"
    with open(mac, "w") as f:
        f.write(f"""/det/shape {shape}
/det/padding 1.0 mm
/run/initialize
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gun/particle gamma
/gun/energy 0.5 MeV
/det/sourceAngle {angle} deg
/det/sourceDistance {dist} cm
/det/openFile {base}
/run/beamOn {events}
/det/writeFile
""")
    subprocess.run([os.path.abspath(exe), f"filt_{tag}.mac"], cwd=builddir,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    csv = os.path.join(builddir, f"{base}_nt_pixels.csv")
    e = read_aggregated(csv)
    os.remove(csv); os.remove(mac)
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="build/exampleB1")
    ap.add_argument("--builddir", default="build")
    ap.add_argument("--shape", required=True)
    ap.add_argument("--seg", type=int, default=64)
    ap.add_argument("--far", type=float, default=400.0, help="far-field distance (cm)")
    ap.add_argument("--near", type=float, default=40.0, help="near-field distance (cm)")
    ap.add_argument("--events", type=int, default=20000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    chans = [("far", args.far), ("near", args.near)]
    filt = None
    for ci, (cname, dist) in enumerate(chans):
        for k in range(args.seg):
            ang = 360.0 * k / args.seg
            e = sim_one(args.exe, args.builddir, args.shape, ang, dist,
                        args.events, f"{cname}_{k}")
            if filt is None:
                filt = np.zeros((2, args.seg, len(e)))
            # Per-sector normalize to unit sum so the template is a shape, not a scale.
            s = e.sum()
            filt[ci, k] = e / s if s > 0 else e
            print(f"  {cname} sector {k+1}/{args.seg} (ang={ang:.1f})", flush=True)

    # Drop crystal columns that are zero across ALL sectors/channels (unused by
    # the shape), so the filter width matches the active crystal count -- exactly
    # how train_compare.py trims the dataset feature columns.
    active = filt.reshape(-1, filt.shape[2]).sum(axis=0) > 0
    filt = filt[:, :, active]

    np.save(args.out, filt)
    print(f"\nSaved filter tensor {filt.shape} -> {args.out} "
          f"({filt.shape[2]} active crystals)")


if __name__ == "__main__":
    main()
