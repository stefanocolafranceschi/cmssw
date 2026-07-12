#!/usr/bin/env python3
# step3_timing.py -- step-9 timing job. One wrapper, three modes, selected with
# the TIMING_MODE environment variable:
#
#   weightsGPU  : stock alpaka reco (weights device producer)          [baseline]
#   multifitGPU : multifit device producer REPLACING the weights at the
#                 same label (customisePhase2MultifitGPU, step 4)
#   multifitCPU : CPU multifit (ecalMultiFitUncalibRecHitPh2Legacy from
#                 the step-8 cff) added as its own Path
#
# Statistics: the 2-event step2.root is read NREPEAT times (default 50 ->
# 100 events; duplicate check off). VALIDATION/DQM/output paths are dropped
# from the schedule, so the job runs RAW2DIGI+RECO only; FastTimerService
# writes resources_<mode>.json.
#
# Run (from the 34434.612_* directory):
#   TIMING_MODE=weightsGPU  cmsRun step3_timing.py
#   TIMING_MODE=multifitGPU cmsRun step3_timing.py
#   TIMING_MODE=multifitCPU cmsRun step3_timing.py
# then:
#   python3 ~/ecal/step9/parse_timing.py resources_*.json
#
# Caveat: for alpaka modules FastTimerService measures host-side module time
# (acquire+produce); device kernels overlap other work. Module numbers are
# still comparable, and events/s is the honest end-to-end metric.

import os

_MODE = os.environ.get("TIMING_MODE", "")
_MODES = ("weightsGPU", "multifitGPU", "multifitCPU")
if _MODE not in _MODES:
    raise RuntimeError("set TIMING_MODE to one of %s" % (_MODES,))
_NREPEAT = int(os.environ.get("NREPEAT", "50"))

_CFG = "step3_RAW2DIGI_RECO_VALIDATION_DQM.py"
if not os.path.exists(_CFG):
    raise RuntimeError(_CFG + " not found -- run from the 34434.612_* workflow directory")
exec(open(_CFG).read())

import FWCore.ParameterSet.Config as cms

# ---- source: repeat the input file for statistics ----
process.source.fileNames = cms.untracked.vstring(*(["file:step2.root"] * _NREPEAT))
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
process.maxEvents.input = cms.untracked.int32(-1)

# ---- mode wiring ----
if _MODE == "multifitGPU":
    from RecoLocalCalo.EcalRecProducers.ecalPhase2MultifitGPUCustomise import (
        customisePhase2MultifitGPU,
    )
    process = customisePhase2MultifitGPU(process)
elif _MODE == "multifitCPU":
    process.load("RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHitPh2_cff")
    process.multifitCPUTimingPath = cms.Path(process.ecalMultiFitUncalibRecHitPh2Legacy)

# ---- schedule the ECAL alpaka chain explicitly ----
# The chain runs anyway (pulled by the reconstruction Task), but unscheduled
# modules do not get their own entries in the FastTimerService JSON summary;
# an explicit Path makes them reportable. Modules run once regardless.
_ecalChain = [m for m in ("simEcalUnsuppressedDigisSoA",
                          "ecalUncalibRecHitPhase2SoA",
                          "ecalUncalibRecHitPhase2") if hasattr(process, m)]
if _ecalChain:
    _seq = getattr(process, _ecalChain[0])
    for _m in _ecalChain[1:]:
        _seq = _seq + getattr(process, _m)
    process.ecalTimingPath = cms.Path(_seq)

# ---- schedule: RAW2DIGI + RECO + ECAL chain (+ CPU multifit path) only ----
_keep = [p for p in (getattr(process, "raw2digi_step", None),
                     getattr(process, "reconstruction_step", None),
                     getattr(process, "ecalTimingPath", None),
                     getattr(process, "multifitCPUTimingPath", None)) if p is not None]
process.schedule = cms.Schedule(*_keep)
for _name in list(process.outputModules_()):
    delattr(process, _name)

# ---- FastTimerService ----
from HLTrigger.Timer.FastTimerService_cfi import FastTimerService
process.FastTimerService = FastTimerService.clone(
    printEventSummary=False,
    printRunSummary=False,
    printJobSummary=True,
    writeJSONSummary=True,
    jsonFileName="resources_%s.json" % _MODE,
    enableDQM=False,
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50
