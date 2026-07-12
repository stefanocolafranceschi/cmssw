#!/usr/bin/env python3
# plot_timing.py -- bar chart from timing_summary.csv (parse_timing.py output).
# Writes timing_comparison.png and timing_comparison.pdf (the PDF is included
# by the LaTeX report).
#
# Usage: python3 ~/ecal/step9/plot_timing.py [timing_summary.csv]

import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = sys.argv[1] if len(sys.argv) > 1 else "timing_summary.csv"

MODULE_COLS = {
    "simEcalUnsuppressedDigisSoA": "digi→SoA",
    "ecalUncalibRecHitPhase2SoA": "uncalib rechit\nproducer",
    "ecalUncalibRecHitPhase2": "SoA→legacy",
    "ecalMultiFitUncalibRecHitPh2Legacy": "CPU multifit",
}
MODE_STYLE = {  # mode -> (legend, color)
    "weightsGPU": ("weights (GPU)", "#7f7f7f"),
    "multifitGPU": ("multifit (GPU)", "#1f77b4"),
    "multifitCPU": ("multifit (CPU)", "#d62728"),
}

rows = list(csv.DictReader(open(CSV)))
if not rows:
    sys.exit("empty " + CSV)

fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [3, 1]})

# panel 1: per-module ms/event
labels = list(MODULE_COLS.values())
n_modes = len(rows)
width = 0.8 / n_modes
for i, r in enumerate(rows):
    mode = r["mode"]
    name, color = MODE_STYLE.get(mode, (mode, None))
    vals, xs = [], []
    for k, col in enumerate(MODULE_COLS):
        v = r.get(col, "")
        if v not in ("", "-", None):
            vals.append(float(v))
            xs.append(k + (i - (n_modes - 1) / 2) * width)
    ax1.bar(xs, vals, width=width, label=name, color=color)
ax1.set_xticks(range(len(labels)))
ax1.set_xticklabels(labels, fontsize=9)
ax1.set_ylabel("module real time per event [ms]")
ax1.set_yscale("log")
ax1.legend(fontsize=9)
ax1.set_title("ECAL Phase-2 local reco, per-module times", fontsize=10)

# panel 2: ECAL uncalib-rechit chain time/event (GPU chain vs CPU multifit).
# Job time is NOT plotted: full-reco + I/O noise dwarfs the ECAL difference.
vals = []
for r in rows:
    mode = r["mode"]
    name, color = MODE_STYLE.get(mode, (mode, None))
    v = r.get("ecal_gpu_chain_ms_per_event", "")
    if mode == "multifitCPU":
        v = r.get("ecalMultiFitUncalibRecHitPh2Legacy", "")
        name += "\n(module)"
    if v in ("", "-", None):
        continue
    vals.append((name, color, float(v)))
for i, (name, color, v) in enumerate(vals):
    ax2.bar(i, v, color=color)
    ax2.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
ax2.set_xticks(range(len(vals)))
ax2.set_xticklabels([n for n, _, _ in vals], fontsize=8)
ax2.set_yscale("log")
ax2.set_ylabel("ECAL uncalib-rechit chain [ms/event]")
ax2.set_title("ECAL chain only", fontsize=10)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig("timing_comparison." + ext, dpi=160)
print("wrote timing_comparison.png / .pdf")
