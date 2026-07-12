#!/usr/bin/env bash
# reproduce_timing.sh -- one-shot CPU-vs-GPU Phase-2 ECAL multifit timing.
#
# Prereqs: a built CMSSW_20_0_0_pre1 with branch from-CMSSW_20_0_0_pre1 merged,
#          `cmsenv` done, and a GPU visible (`nvidia-smi`).
# Run (from anywhere):
#   bash phase2-multifit-timing/reproduce_timing.sh
#
# What it does:
#   1. makes step2.root via one matrix pass for workflow 34434.612 if absent
#      (ttbar, Run4 D121, no pileup; a few minutes, needs the GPU),
#   2. times three modes over 100 events (50 replays of the 2-event file):
#        weightsGPU  -- stock alpaka weights reco            [baseline]
#        multifitGPU -- multifit device producer, same label
#        multifitCPU -- the CPU multifit worker, own Path
#   3. prints the per-module and ECAL-chain ms/event table.
# Expect: multifit-CPU ~260 ms/ev, GPU chain ~6.4 ms/ev (x41); the multifit-GPU
# producer ~= the weights-GPU producer (the NNLS upgrade is free at module level).
set -eu

HERE=$(cd "$(dirname "$0")" && pwd)
: "${CMSSW_BASE:?run cmsenv first}"
command -v cmsRun >/dev/null || { echo "cmsRun not found -- run cmsenv first"; exit 1; }
command -v nvidia-smi >/dev/null && nvidia-smi -L || echo "(no nvidia-smi; GPU modes need a visible GPU)"

cd "$CMSSW_BASE/src"
WF=$(ls -d 34434.612_*/ 2>/dev/null | head -1 || true)
if [ -z "${WF:-}" ] || [ ! -f "${WF}step2.root" ]; then
  echo "== no step2.root yet -- one matrix pass for 34434.612 (GEN..RECO) =="
  runTheMatrix.py -w upgrade -l 34434.612
  WF=$(ls -d 34434.612_*/ | head -1)
fi
cd "$WF"
echo "== workflow dir: $PWD =="
[ -f step2.root ] || { echo "step2.root still missing -- check the matrix log above"; exit 1; }

cp "$HERE/step3_timing.py" .
for M in weightsGPU multifitGPU multifitCPU; do
  echo "== timing mode: $M (100 events) =="
  TIMING_MODE="$M" cmsRun step3_timing.py > "timing_${M}.log" 2>&1 \
    || { echo "  *** $M failed -- tail timing_${M}.log:"; tail -12 "timing_${M}.log"; exit 1; }
done

echo
echo "======================== RESULTS ========================"
python3 "$HERE/parse_timing.py" resources_*.json
echo "========================================================="
echo "Chain-to-chain speedup = (multifitCPU job) / (multifitGPU ecalGPUchain)."
echo "Module/chain ratios are stable run-to-run; the job wall time is noisy on a shared box."
