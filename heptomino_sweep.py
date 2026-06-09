#!/usr/bin/env python3
"""
Enumerate all free heptominoes (7-cell polyominoes) and run the full
crystal-count sweep for N=7 on this machine. Intended for the RTX 4070 Windows
laptop (GPU training is auto-used by train_compare.py).

It:
  1. enumerates the 108 free heptominoes (canonical form, dedup over the 8
     rotations/reflections),
  2. writes a Geant4 macro per heptomino using the runtime /det/cells command
     (no recompile needed),
  3. runs the simulation to produce one dataset per shape,
  4. trains + ranks them with train_compare.py to find the best 7-crystal shape.

Usage (from the project root, after building exampleB1):
  # Windows: exe is build\\exampleB1.exe ; macOS/Linux: build/exampleB1
  python3 heptomino_sweep.py --exe build/exampleB1 --configs 1000 --events 5000

Then train:
  python3 train_compare.py build/dsHep_*_nt_pixels.csv

The script prints the count of heptominoes found (should be 108) before running.
"""
import argparse
import os
import subprocess
import sys


# ---------------------------------------------------------------------------
# Free-polyomino enumeration
# ---------------------------------------------------------------------------
def normalize(cells):
    """Translate so min col/row = 0; return sorted tuple (canonical translation)."""
    mc = min(c for c, r in cells)
    mr = min(r for c, r in cells)
    return tuple(sorted((c - mc, r - mr) for c, r in cells))


def transforms(cells):
    """All 8 rotations/reflections, each normalized."""
    out = []
    pts = list(cells)
    for _ in range(4):
        pts = [(r, -c) for c, r in pts]          # rotate 90
        out.append(normalize(pts))
        out.append(normalize([(-c, r) for c, r in pts]))  # reflect
    return out


def canonical(cells):
    """Lexicographically smallest of all 8 transforms -> unique key per free shape."""
    return min(transforms(cells))


def free_polyominoes(n):
    """Generate all free n-cell polyominoes by growth from the monomino."""
    current = {((0, 0),)}
    for _ in range(n - 1):
        nxt = set()
        for shape in current:
            cellset = set(shape)
            for c, r in shape:
                for dc, dr in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    nb = (c + dc, r + dr)
                    if nb not in cellset:
                        nxt.add(canonical(cellset | {nb}))
        current = nxt
    return sorted(current)


# ---------------------------------------------------------------------------
# Macro generation + run
# ---------------------------------------------------------------------------
def cell_spec(cells):
    return ";".join(f"{c},{r}" for c, r in cells)


def write_macro(path, spec, out_base, configs, events):
    with open(path, "w") as f:
        f.write(f"""/det/cells {spec}
/det/padding 1.0 mm
/run/initialize
/run/verbose 0
/event/verbose 0
/tracking/verbose 0
/gun/particle gamma
/gun/energy 0.5 MeV
/det/openFile {out_base}
/control/loop gen_config_step.mac i 1 {configs} 1
/det/writeFile
""")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="build/exampleB1",
                    help="path to the built executable (use build\\exampleB1.exe on Windows)")
    ap.add_argument("--configs", type=int, default=1000)
    ap.add_argument("--events", type=int, default=5000)
    ap.add_argument("--builddir", default="build",
                    help="dir where exe runs and CSVs land (must contain gen_config_step.mac)")
    ap.add_argument("--limit", type=int, default=0,
                    help="only run the first N heptominoes (0 = all 108)")
    args = ap.parse_args()

    heps = free_polyominoes(7)
    print(f"Enumerated {len(heps)} free heptominoes (expected 108).")
    if args.limit:
        heps = heps[:args.limit]
        print(f"Limiting to first {len(heps)}.")

    # gen_config_step.mac must exist in the build dir.
    step = os.path.join(args.builddir, "gen_config_step.mac")
    if not os.path.exists(step):
        print(f"ERROR: {step} not found. Copy gen_config_step.mac into the build dir.")
        sys.exit(1)

    for idx, cells in enumerate(heps):
        tag = f"hep{idx:03d}"
        spec = cell_spec(cells)
        mac = os.path.join(args.builddir, f"gen_{tag}.mac")
        write_macro(mac, spec, f"dsHep_{tag}", args.configs, args.events)
        print(f"[{idx+1}/{len(heps)}] {tag}: {spec}")
        # Run from the build dir so output + gen_config_step.mac resolve.
        subprocess.run([os.path.abspath(args.exe), f"gen_{tag}.mac"],
                       cwd=args.builddir,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print("\nDone generating datasets. Now train + rank with:")
    print(f"  python3 train_compare.py {args.builddir}/dsHep_*_nt_pixels.csv")


if __name__ == "__main__":
    main()
