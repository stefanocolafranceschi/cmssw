import FWCore.ParameterSet.Config as cms

# In the alpaka-wired Phase-2 reco (ecalUncalibRecHitPhase2_cff with the alpaka
# process modifier):
#   simEcalUnsuppressedDigisSoA   : digi -> portable SoA producer
#   ecalUncalibRecHitPhase2SoA    : weights producer on device
#   ecalUncalibRecHitPhase2       : SoA -> legacy converter (consumed downstream)
#
# Since step 8 the conditions content (measured pulse template / sample
# correlation / pulse covariance + EmptyESSources) lives in
# ecalPhase2MultifitTrivialConditions_cff.py, shared with
# ecalMultiFitUncalibRecHitPh2_cff.py. This file only wires producers into an
# existing workflow config; entry points and module labels are unchanged
# since step 4 (step3_cpuVgpu.py and the step-6 comparator rely on them).

_DIGIS = "simEcalUnsuppressedDigisSoA:ebDigis"


def _loadPhase2MultifitModules(process):
    # trivial providers for the five Phase-2 records (no DB tags exist,
    # checked 150X_mcRun4_realistic_v1, Jul 2026) ...
    process.load("RecoLocalCalo.EcalRecProducers.ecalPhase2MultifitTrivialConditions_cff")
    # ... and the ES producer filling EcalMultifitConditionsPhase2Rcd from them
    process.load("RecoLocalCalo.EcalRecProducers.ecalMultifitConditionsPhase2HostESProducer_cfi")


def customisePhase2MultifitGPU(process):
    """Replace the Phase-2 weights DEVICE producer with the alpaka multifit
    producer at the same module label (ecalUncalibRecHitPhase2SoA). The
    SoA->legacy converter and everything downstream are untouched and consume
    the multifit result. Compare against a stock (weights) job to validate."""
    if not hasattr(process, "ecalUncalibRecHitPhase2SoA"):
        raise RuntimeError(
            "customisePhase2MultifitGPU: no 'ecalUncalibRecHitPhase2SoA' module in the process; "
            "this needs the alpaka-wired Phase-2 reco (ecalDevelAlpaka workflow)")

    _loadPhase2MultifitModules(process)

    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2ProducerPortable_cfi import (
        ecalUncalibRecHitPhase2ProducerPortable as _multifit,
    )
    process.globalReplace("ecalUncalibRecHitPhase2SoA", _multifit.clone(digisLabelEB=_DIGIS))
    return process


def customisePhase2MultifitGPUAlongside(process):
    """Run the alpaka multifit producer in addition to the weights producer,
    convert its output to the legacy collection format, and keep it in the
    event output for a same-file comparison:
      weights : ecalUncalibRecHitPhase2 : EcalUncalibRecHitsEB
      multifit: ecalMultiFitUncalibRecHitPhase2 : EcalUncalibRecHitsEB"""
    if not hasattr(process, "simEcalUnsuppressedDigisSoA"):
        raise RuntimeError(
            "customisePhase2MultifitGPUAlongside: no 'simEcalUnsuppressedDigisSoA' module in the "
            "process; this needs the alpaka-wired Phase-2 reco (ecalDevelAlpaka workflow)")

    _loadPhase2MultifitModules(process)

    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2ProducerPortable_cfi import (
        ecalUncalibRecHitPhase2ProducerPortable as _multifit,
    )
    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitSoAToLegacy_cfi import (
        ecalUncalibRecHitSoAToLegacy as _soaToLegacy,
    )

    process.ecalMultiFitUncalibRecHitPhase2Portable = _multifit.clone(
        digisLabelEB=_DIGIS,
        recHitsLabelEB="EcalMultiFitUncalibRecHitsEB",
    )
    process.ecalMultiFitUncalibRecHitPhase2 = _soaToLegacy.clone(
        isPhase2=True,
        inputCollectionEB="ecalMultiFitUncalibRecHitPhase2Portable:EcalMultiFitUncalibRecHitsEB",
        inputCollectionEE=None,
        outputLabelEE=None,
    )

    process.ecalMultifitPhase2Task = cms.Task(
        process.ecalMultiFitUncalibRecHitPhase2Portable,
        process.ecalMultiFitUncalibRecHitPhase2,
    )
    if hasattr(process, "reconstruction_step"):
        process.reconstruction_step.associate(process.ecalMultifitPhase2Task)
    else:
        raise RuntimeError("customisePhase2MultifitGPUAlongside: no 'reconstruction_step' in the process")

    for out in process.outputModules_().values():
        if out.type_() == "PoolOutputModule":
            out.outputCommands.append("keep *_ecalMultiFitUncalibRecHitPhase2_*_*")
    return process
