#!/usr/bin/env python3
"""
Aggregate the per-event Geant4 readout into a per-angle results table.

Input : a per-event ntuple CSV from the detector
        (columns: eventID, e0_keV, e1_keV, e2_keV, e3_keV, angle_deg, dist_cm)
Output: a results CSV with ONE row per angle:
        angle, p0, p1, p2, p3, p0_norm, p1_norm, p2_norm, p3_norm

  p0..p3      = total energy (keV) summed over all events at that angle
  p0_norm..   = each pixel's fraction of the 4-pixel total (sums to 1)

Usage:
  python3 aggregate_results.py [input.csv] [output.csv]
  defaults: build/dataset_nt_pixels_t0.csv -> build/results_by_angle.csv
"""
import sys
import csv
from collections import defaultdict

IN = sys.argv[1] if len(sys.argv) > 1 else "build/dataset_nt_pixels_t0.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "build/results_by_angle.csv"

# Sum per-pixel energy per (rounded) angle.
sums = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
counts = defaultdict(int)

with open(IN) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue  # skip the G4 ntuple header block
        parts = line.split(",")
        if len(parts) < 6:
            continue
        try:
            e = [float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])]
            ang = round(float(parts[5]))  # collapse float wobble (e.g. 60.0001 -> 60)
        except ValueError:
            continue
        s = sums[ang]
        for i in range(4):
            s[i] += e[i]
        counts[ang] += 1

rows = []
for ang in sorted(sums):
    p = sums[ang]
    total = sum(p)
    if total > 0:
        norm = [v / total for v in p]
    else:
        norm = [0.0, 0.0, 0.0, 0.0]
    rows.append([ang] + p + norm)

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["angle", "p0", "p1", "p2", "p3",
                "p0_norm", "p1_norm", "p2_norm", "p3_norm"])
    for r in rows:
        w.writerow([r[0]] + [f"{v:.4f}" for v in r[1:]])

print(f"Wrote {len(rows)} angle rows to {OUT}")
print(f"Total events aggregated: {sum(counts.values())}")
print("\nPreview:")
print("angle   p0       p1       p2       p3      | p0_n  p1_n  p2_n  p3_n")
for r in rows[::3]:  # every 3rd angle
    print(f"{r[0]:4d}  {r[1]:8.0f} {r[2]:8.0f} {r[3]:8.0f} {r[4]:8.0f} | "
          f"{r[5]:.2f}  {r[6]:.2f}  {r[7]:.2f}  {r[8]:.2f}")
