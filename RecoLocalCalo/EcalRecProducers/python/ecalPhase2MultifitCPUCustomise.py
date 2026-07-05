import FWCore.ParameterSet.Config as cms

# Same trivial Phase-2 conditions as the GPU multifit (5 records + EmptyESSources)
from RecoLocalCalo.EcalRecProducers.ecalPhase2MultifitGPUCustomise import _loadPhase2MultifitModules

# Legacy Phase-2 digis (EBDigiCollectionPh2), same default as the CPU weights
# producer and the input of the digi->SoA converter.
_LEGACY_DIGIS = "simEcalUnsuppressedDigis"


def customisePhase2MultifitCPUAlongside(process):
    """Run Thomas's CPU Phase-2 multifit worker alongside whatever else is in
    the job, on the SAME trivial conditions as the GPU multifit:

        CPU multifit: ecalMultiFitUncalibRecHitPh2CPU : EcalUncalibRecHitsEB

    Apply AFTER customisePhase2MultifitGPUAlongside for a same-file
    weights / GPU-multifit / CPU-multifit comparison.

    Config notes (validated against the branch sources, Jul 2026):
    - timealgo must be one of RatioMethod/WeightsMethod/crossCorrelationMethod
      (the ParameterSwitch has no "None" case). crossCorrelationMethod is the
      amplitude-only choice: it esConsumes NOTHING extra and its Phase-2 body
      is commented out upstream (jitter set to 0). RatioMethod would pull 4
      more records (time calib/offset/bias + sample mask) that have no source
      here; WeightsMethod needs EcalPh2TBWeights, which only Thomas's
      EcalTrivialConditionRetriever provides (clashes with our conditions).
    - useLumiInfoRunHeader=False + bunchSpacing=0 -> the activeBXs parameter
      {-2..2} is used as-is (no bunchSpacingProducer dependency); identical
      to the GPU kernel's 5 pulses.
    - ampErrorCalculation=False mirrors the GPU (no amplitude error kernel).
    - gainSwitchUseMaxSample/simplifiedNoiseModelForGainSwitch are inert:
      the Phase-2 CPU algo hardcodes hasGainSwitch=false."""
    if hasattr(process, "ecalMultiFitUncalibRecHitPh2CPU"):
        print("customisePhase2MultifitCPUAlongside: already applied, skipping")
        return process
    if not hasattr(process, "reconstruction_step"):
        raise RuntimeError(
            "customisePhase2MultifitCPUAlongside: no 'reconstruction_step' in the process")

    # conditions: reuse the GPU customise's trivial providers if already there
    if not hasattr(process, "ecalPhase2TrivialCondESProducer"):
        _loadPhase2MultifitModules(process)

    process.ecalMultiFitUncalibRecHitPh2CPU = cms.EDProducer(
        "EcalUncalibRecHitProducer",
        EBdigiCollection=cms.InputTag(_LEGACY_DIGIS),
        EEdigiCollection=cms.InputTag(""),
        EBhitCollection=cms.string("EcalUncalibRecHitsEB"),
        EEhitCollection=cms.string(""),
        IsPhase2=cms.bool(True),
        algo=cms.string("EcalUncalibRecHitWorkerMultiFitPh2"),
        algoPSet=cms.PSet(
            activeBXs=cms.vint32(-2, -1, 0, 1, 2),
            ampErrorCalculation=cms.bool(False),
            useLumiInfoRunHeader=cms.bool(False),
            bunchSpacing=cms.int32(0),
            doPrefit=cms.bool(False),
            prefitMaxChiSq=cms.double(25.0),
            dynamicPedestals=cms.bool(False),
            mitigateBadSamples=cms.bool(False),
            gainSwitchUseMaxSample=cms.bool(False),
            selectiveBadSampleCriteria=cms.bool(False),
            addPedestalUncertainty=cms.double(0.0),
            simplifiedNoiseModelForGainSwitch=cms.bool(True),
            timealgo=cms.string("crossCorrelationMethod"),
            # crossCorrelation* parameters: tracked defaults, inert for Phase-2
        ),
    )

    process.ecalMultifitPhase2CPUTask = cms.Task(process.ecalMultiFitUncalibRecHitPh2CPU)
    process.reconstruction_step.associate(process.ecalMultifitPhase2CPUTask)

    # the output-module keep is what makes the Task actually run (pull mode)
    for out in process.outputModules_().values():
        if out.type_() == "PoolOutputModule":
            out.outputCommands.append("keep *_ecalMultiFitUncalibRecHitPh2CPU_*_*")
    return process
