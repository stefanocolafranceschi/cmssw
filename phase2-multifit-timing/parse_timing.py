#!/usr/bin/env python3
# parse_timing.py -- summarize FastTimerService JSON files from step3_timing.py.
#
# Usage: python3 ~/ecal/step9/parse_timing.py resources_*.json
# Prints a table and writes timing_summary.csv (input for plot_timing.py).
#
# Reported per mode: total events, job real time per event, and real time per
# event for the ECAL-local-reco modules of interest.

import csv
import json
import sys

# module label -> short name in the table (missing modules are skipped)
MODULES = {
    "simEcalUnsuppressedDigisSoA": "digi->SoA",
    "ecalUncalibRecHitPhase2SoA": "device producer (weights or multifit)",
    "ecalUncalibRecHitPhase2": "SoA->legacy converter",
    "ecalMultiFitUncalibRecHitPh2Legacy": "CPU multifit",
}

def main(paths):
    rows = []
    for path in sorted(paths):
        with open(path) as f:
            j = json.load(f)
        mode = path.split("resources_")[-1].replace(".json", "")
        events = j.get("total", {}).get("events", 0)
        total_ms = j.get("total", {}).get("time_real", 0.0)  # ms, whole job
        per_event = total_ms / events if events else float("nan")
        row = {"mode": mode, "events": events, "job_ms_per_event": round(per_event, 3)}
        chain = 0.0
        for m in j.get("modules", []):
            label = m.get("type", "")
            # FastTimerService JSON: entries carry 'type' and 'label'
            label = m.get("label", label)
            if label in MODULES:
                t = m.get("time_real", 0.0)
                row[label] = round(t / events, 4) if events else float("nan")
                if label != "ecalMultiFitUncalibRecHitPh2Legacy":
                    chain += t / events if events else 0.0
        # apples-to-apples ECAL metric: GPU chain = digi->SoA + producer +
        # converter; CPU = the CPU multifit module itself. Job time is noisy
        # (full reco + I/O) -- do not use it as the headline number.
        row["ecal_gpu_chain_ms_per_event"] = round(chain, 4) if chain else ""
        rows.append(row)

    cols = ["mode", "events", "job_ms_per_event", "ecal_gpu_chain_ms_per_event"] + list(MODULES)
    with open("timing_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"{'mode':<14}{'events':>8}{'job ms/ev':>12}{'ecalGPUchain':>14}", end="")
    for label in MODULES:
        print(f"{label[-34:]:>36}", end="")
    print()
    for r in rows:
        print(f"{r['mode']:<14}{r['events']:>8}{r['job_ms_per_event']:>12}"
              f"{r.get('ecal_gpu_chain_ms_per_event', '') or '-':>14}", end="")
        for label in MODULES:
            print(f"{r.get(label, '-'):>36}", end="")
        print()
    print("\nwrote timing_summary.csv (feed to plot_timing.py)")
    print("NOTE: module times are host-side real times (ms/event). If a module")
    print("column is empty, check the label exists in that mode's JSON:")
    print("  python3 -c \"import json;print([m['label'] for m in json.load(open('resources_<mode>.json'))['modules']])\"")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: parse_timing.py resources_*.json")
    main(sys.argv[1:])
