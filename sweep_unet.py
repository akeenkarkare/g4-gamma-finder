#!/usr/bin/env python3
"""
End-to-end U-Net shape sweep: for each shape, generate its dataset + MC filters,
train the MC-initialized U-Net, and print a ranking. Built for the RTX laptop
(train_unet.py auto-uses CUDA).

Two built-in modes plus a custom list:
  --mode hexominoes   all 35 hexominoes (h01..h35)        -> true best-6 under U-Net
  --mode heptominoes  compact heptominoes (footprint-controlled, via /det/cells)
  --shapes "S pP h27" explicit named shapes

Why footprint-control for N=7: a sprawling heptomino re-introduces the footprint
confound (a bigger array blurs angle). We restrict to heptominoes whose bounding
box stays within --maxspan cells of pP's footprint, so N=7 is comparable to
pP (5) and h13 (6).

Usage (Windows RTX, exe built):
  python sweep_unet.py --exe build\\exampleB1.exe --mode hexominoes
  python sweep_unet.py --exe build\\exampleB1.exe --mode heptominoes --maxspan 3
"""
import argparse
import os
import subprocess
import sys
import numpy as np

# Reuse the verified heptomino enumerator.
from heptomino_sweep import free_polyominoes


def run_geant(exe, builddir, args_list):
    subprocess.run([os.path.abspath(exe)] + args_list, cwd=builddir,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def write_mac(path, body):
    with open(path, "w") as f:
        f.write(body)


def gen_dataset(exe, builddir, shape_cmd, tag, configs, events):
    """shape_cmd is the /det line (either '/det/shape X' or '/det/cells ...')."""
    mac = f"sw_{tag}.mac"
    write_mac(os.path.join(builddir, mac), f"""{shape_cmd}
/det/padding 1.0 mm
/run/initialize
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gun/particle gamma
/gun/energy 0.5 MeV
/det/openFile ds_{tag}
/control/loop gen_config_step.mac i 1 {configs} 1
/det/writeFile
""")
    run_geant(exe, builddir, [mac])
    os.remove(os.path.join(builddir, mac))
    return os.path.join(builddir, f"ds_{tag}_nt_pixels.csv")


def gen_filters(builddir, exe, shape_cmd, tag, seg=64, far=400, near=40, events=15000):
    """Reuse gen_filters logic but allow a /det/cells shape via a small inline run.
    For named shapes we call gen_filters.py; for custom we replicate minimally."""
    out = os.path.join(builddir, f"filt_{tag}.npy")
    # gen_filters.py only supports --shape (named). For custom cells we generate
    # filters here directly using the same per-sector sim.
    if shape_cmd.startswith("/det/shape "):
        shape = shape_cmd.split()[-1]
        subprocess.run([sys.executable, "gen_filters.py", "--exe", exe,
                        "--shape", shape, "--seg", str(seg), "--far", str(far),
                        "--near", str(near), "--events", str(events), "--out", out],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        _gen_filters_custom(builddir, exe, shape_cmd, out, seg, far, near, events)
    return out


def _read_e(csv):
    cols, row = [], None
    for line in open(csv):
        line = line.rstrip("\n")
        if line.startswith("#column"):
            cols.append(line.split()[-1])
        elif line and not line.startswith("#") and not line.startswith("entries"):
            row = [float(x) for x in line.split(",")]
    e_idx = [i for i, c in enumerate(cols) if c.startswith("e") and c.endswith("_keV")]
    return np.array([row[i] for i in e_idx])


def _gen_filters_custom(builddir, exe, cells_cmd, out, seg, far, near, events):
    filt = None
    for ci, dist in enumerate([far, near]):
        for k in range(seg):
            ang = 360.0 * k / seg
            tag = f"flt_{ci}_{k}"
            write_mac(os.path.join(builddir, f"{tag}.mac"), f"""{cells_cmd}
/det/padding 1.0 mm
/run/initialize
/run/verbose 0
/gun/particle gamma
/gun/energy 0.5 MeV
/det/sourceAngle {ang} deg
/det/sourceDistance {dist} cm
/det/openFile {tag}
/run/beamOn {events}
/det/writeFile
""")
            run_geant(exe, builddir, [f"{tag}.mac"])
            e = _read_e(os.path.join(builddir, f"{tag}_nt_pixels.csv"))
            os.remove(os.path.join(builddir, f"{tag}_nt_pixels.csv"))
            os.remove(os.path.join(builddir, f"{tag}.mac"))
            if filt is None:
                filt = np.zeros((2, seg, len(e)))
            s = e.sum()
            filt[ci, k] = e / s if s > 0 else e
    active = filt.reshape(-1, filt.shape[2]).sum(axis=0) > 0
    np.save(out, filt[:, :, active])


def compact_heptominoes(maxspan):
    """Heptominoes whose bounding box fits within (maxspan+1) x (maxspan+1) cells."""
    out = []
    for cells in free_polyominoes(7):
        cols = [c for c, r in cells]; rows = [r for c, r in cells]
        if (max(cols) - min(cols)) <= maxspan and (max(rows) - min(rows)) <= maxspan:
            out.append(cells)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="build/exampleB1")
    ap.add_argument("--builddir", default="build")
    ap.add_argument("--mode", choices=["hexominoes", "heptominoes", "named"], default="named")
    ap.add_argument("--shapes", default="", help="space-separated named shapes (mode=named)")
    ap.add_argument("--maxspan", type=int, default=2,
                    help="heptomino bbox span limit in cells (2 => 3x3 box, the "
                         "footprint-controlled set ~comparable to pP/h13; 7 shapes)")
    ap.add_argument("--configs", type=int, default=1000)
    ap.add_argument("--events", type=int, default=5000)
    args = ap.parse_args()

    # Build the work list as (tag, /det command).
    jobs = []
    if args.mode == "hexominoes":
        jobs = [(f"h{n:02d}", f"/det/shape h{n:02d}") for n in range(1, 36)]
    elif args.mode == "heptominoes":
        heps = compact_heptominoes(args.maxspan)
        print(f"{len(heps)} compact heptominoes (bbox <= {args.maxspan+1} cells/side)")
        for i, cells in enumerate(heps):
            spec = ";".join(f"{c},{r}" for c, r in cells)
            jobs.append((f"hep{i:03d}", f"/det/cells {spec}"))
    else:
        jobs = [(s, f"/det/shape {s}") for s in args.shapes.split()]

    results = {}
    for tag, cmd in jobs:
        print(f"\n=== {tag} ({cmd}) ===", flush=True)
        ds = gen_dataset(args.exe, args.builddir, cmd, tag, args.configs, args.events)
        flt = gen_filters(args.builddir, args.exe, cmd, tag, events=15000)
        # Train via train_unet.py (auto-CUDA) and capture the mean.
        p = subprocess.run([sys.executable, "train_unet.py", "--data", ds,
                            "--filters", flt, "--label", tag],
                           capture_output=True, text=True)
        line = [l for l in p.stdout.splitlines() if "mean =" in l]
        if line:
            print(line[-1].strip())
            results[tag] = float(line[-1].split("mean =")[1].split("deg")[0])

    print("\n================ U-NET RANKING (lower = better) ================")
    for tag, deg in sorted(results.items(), key=lambda kv: kv[1]):
        print(f"  {tag:10s} {deg:.2f} deg")


if __name__ == "__main__":
    main()
