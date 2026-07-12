# Phase-2 ECAL multifit — CPU-vs-GPU timing reproduction

One command reproduces the CPU-vs-GPU multifit timing on workflow **34434.612**
(ttbar, Run4 D121, no pileup). Suggested location in the fork:
`phase2-multifit-timing/` at the repo root (a plain folder, not a CMSSW package,
so `scram` ignores it).

## From scratch

```bash
cmsrel CMSSW_20_0_0_pre1 && cd CMSSW_20_0_0_pre1/src && cmsenv
git cms-merge-topic stefanocolafranceschi:from-CMSSW_20_0_0_pre1
scram b -j$(nproc)

bash phase2-multifit-timing/reproduce_timing.sh
```

The driver makes `step2.root` with one matrix pass if it is missing, times three
modes over 100 events (50 replays of the 2-event file), and prints the table.

## What is measured

| mode | reco |
|------|------|
| `weightsGPU`  | stock alpaka weights device producer (baseline) |
| `multifitGPU` | multifit device producer at the same label |
| `multifitCPU` | CPU multifit worker in its own Path |

Reported per mode: total events, job ms/event, and host-side ms/event for the
ECAL-local-reco modules (`digi→SoA`, device producer, `SoA→legacy` converter,
CPU multifit). The comparison is chain-to-chain: the GPU chain
(digi→SoA + producer + converter) vs the single CPU multifit module.

## Expected (RTX 3080, 4 threads / 4 streams, 100 events)

- CPU multifit ≈ **260 ms/event**, GPU chain ≈ **6.4 ms/event** → **×41**.
- multifit-GPU producer ≈ weights-GPU producer (≈1.16 vs 1.13 ms/eve)

## Files

- `reproduce_timing.sh` — the driver (run this).
- `step3_timing.py` — the timing job (`TIMING_MODE` env selects the mode; `NREPEAT` sets replays, default 50).
- `parse_timing.py` — prints the table and writes `timing_summary.csv`.
- `plot_timing.py` — optional per-module + chain bar chart.
