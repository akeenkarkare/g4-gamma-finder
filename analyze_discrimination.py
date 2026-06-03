#!/usr/bin/env python3
"""
Demonstrate background/source discrimination for the LaBr3 detector.

Reads the per-crystal energy spectra (G4 H1 CSVs) for:
  - an external Eu-152 + Ba-133 source run (directional), and
  - an internal 138La background run (intrinsic LaBr3 radioactivity),
then shows the two discrimination levers:

  1. ENERGY: the source lives in a low-energy ROI; 138La sits at 789/1436 keV.
     A low-energy window keeps the source and rejects most internal background.
  2. SYMMETRY: in-ROI, the source is directional (one pixel dominates) while the
     internal background is ~uniform across pixels -> cancels in pixel contrast.

Usage:
  python3 analyze_discrimination.py build/mix_src build/mix_int [ROI_lo ROI_hi]
  defaults: ROI = [40, 450] keV
"""
import sys


def load_spectrum(prefix, pixel):
    """Return a list counts[energy_keV] from a G4 H1 CSV."""
    path = f"{prefix}_h1_spectrum{pixel}.csv"
    counts = []
    started = False
    for line in open(path):
        if line.startswith("#"):
            continue
        if line.startswith("entries"):
            started = True
            continue
        if started:
            counts.append(float(line.split(",")[0]))
    # counts[0] = underflow; counts[i] ~ energy (i-1) keV
    return counts


def roi_sum(counts, lo, hi):
    s = 0.0
    for i in range(1, len(counts)):
        e = i - 1
        if lo <= e <= hi:
            s += counts[i]
    return s


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "build/mix_src"
    intl = sys.argv[2] if len(sys.argv) > 2 else "build/mix_int"
    lo = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
    hi = float(sys.argv[4]) if len(sys.argv) > 4 else 450.0

    print(f"ROI = [{lo:.0f}, {hi:.0f}] keV\n")

    for label, prefix in [("EXTERNAL Eu/Ba source", src),
                          ("INTERNAL 138La background", intl)]:
        per_pix_in, per_pix_all = [], []
        for p in range(4):
            c = load_spectrum(prefix, p)
            per_pix_in.append(roi_sum(c, lo, hi))
            per_pix_all.append(sum(c[1:]))
        tot_in = sum(per_pix_in)
        tot_all = sum(per_pix_all)
        frac_in = tot_in / tot_all if tot_all else 0.0

        # Directionality metric: max pixel fraction of the in-ROI total.
        if tot_in > 0:
            fr = [v / tot_in for v in per_pix_in]
            dom = max(range(4), key=lambda i: fr[i])
        else:
            fr = [0, 0, 0, 0]; dom = -1

        print(f"--- {label} ---")
        print(f"  counts in ROI / total : {tot_in:.0f} / {tot_all:.0f} "
              f"({100*frac_in:.1f}% in ROI)")
        print(f"  in-ROI pixel fractions: "
              f"p0={fr[0]:.2f} p1={fr[1]:.2f} p2={fr[2]:.2f} p3={fr[3]:.2f}")
        if dom >= 0:
            print(f"  dominant pixel        : p{dom} ({fr[dom]:.2f})")
        print()

    print("Interpretation:")
    print("  - The source puts most of its counts INSIDE the low-energy ROI and")
    print("    is strongly DIRECTIONAL (one pixel dominates).")
    print("  - The 138La background is mostly OUTSIDE the ROI (its lines are at")
    print("    789/1436 keV) and what leaks in is ~UNIFORM across pixels, so it")
    print("    contributes little to the pixel contrast the direction net uses.")


if __name__ == "__main__":
    main()
