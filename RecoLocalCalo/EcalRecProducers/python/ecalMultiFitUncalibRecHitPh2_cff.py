import FWCore.ParameterSet.Config as cms

# Phase-2 ECAL barrel multifit amplitude reconstruction, CPU default with an
# alpaka replacement, mirroring the structure of the Run-3
# ecalMultiFitUncalibRecHit_cff.py and of the Phase-2 weights
# ecalUncalibRecHitPhase2_cff.py.
#
# This file replaces the earlier SwitchProducerCUDA-based version. Changes:
#  - SwitchProducerCUDA (removed from CMSSW_20) -> alpaka.toReplaceWith
#  - EcalTrivialConditionRetriever dropped: its Ph2Weights / Ph2SampleMask /
#    Ph2SamplesCorrelation payloads are not consumed on the amplitude-only
#    path (timealgo = crossCorrelationMethod esConsumes nothing extra), and
#    its SamplesCorrelation would clash with the trivial conditions below.
#  - the inline 2017 test-beam-sim pulse shape / covariance ESProducers
#    dropped: that template (peak index 5) is misaligned with the current
#    digitizer output; conditions come from
#    ecalPhase2MultifitTrivialConditions_cff (template measured from digis,
#    peak index 1) until DB tags exist.

# interim conditions for the five Phase-2 records (LiteDTU pedestals, CATIA
# gain ratios, pulse shape, pulse covariance, samples correlation), needed by
# both the CPU and the alpaka path.
# TODO: drop this import once CondTools writers + GlobalTag tags exist.
from RecoLocalCalo.EcalRecProducers.ecalPhase2MultifitTrivialConditions_cff import *

# ECAL Phase-2 multifit running on CPU.
# Overrides wrt the cfi (validated in the step-6 CPU-vs-GPU comparison):
#  - EBdigiCollection: Phase-2 MC digis (EBDigiCollectionPh2), not the
#    Phase-1 unpacker output
#  - timealgo = crossCorrelationMethod: the only amplitude-only choice; the
#    ParameterSwitch has no "None" case, RatioMethod would esConsume four
#    more record-less products, WeightsMethod needs EcalPh2TBWeights. Its
#    Phase-2 body is commented out upstream (jitter = 0). Because algoPSet
#    is validated through a ParameterSwitch on timealgo, the RatioMethod
#    parameters carried by the cfi are ILLEGAL (not merely inert) under any
#    other case and must be deleted in the clone (param = None deletes).
#  - useLumiInfoRunHeader = False + bunchSpacing = 0: activeBXs {-2..2} used
#    as-is, no bunchSpacingProducer dependency; same 5 pulses as the GPU.
#  - ampErrorCalculation = False mirrors the GPU (no amplitude error kernel).
#  - gainSwitchUseMaxSample = False (inert: Phase-2 algo hardcodes
#    hasGainSwitch = false).
from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHitPh2_cfi import ecalMultiFitUncalibRecHitPh2 as _ecalMultiFitUncalibRecHitPh2
ecalMultiFitUncalibRecHitPh2 = _ecalMultiFitUncalibRecHitPh2.clone(
    EBdigiCollection = 'simEcalUnsuppressedDigis',
    algoPSet = dict(
        ampErrorCalculation = False,
        useLumiInfoRunHeader = False,
        bunchSpacing = cms.int32(0),
        gainSwitchUseMaxSample = False,
        timealgo = 'crossCorrelationMethod',
        # delete the RatioMethod-case parameters (illegal under
        # crossCorrelationMethod, see above)
        timeFitParameters = None,
        amplitudeFitParameters = None,
        timeFitLimits_Lower = None,
        timeFitLimits_Upper = None,
        timeConstantTerm = None,
        timeNconst = None,
        outOfTimeThresholdGain10p = None,
        outOfTimeThresholdGain10m = None,
        outOfTimeThresholdGain1p = None,
        outOfTimeThresholdGain1m = None,
        amplitudeThreshold = None,
    )
)
ecalMultiFitUncalibRecHitPh2Legacy = ecalMultiFitUncalibRecHitPh2.clone()

ecalMultiFitUncalibRecHitPh2Task = cms.Task(
    # ECAL Phase-2 multifit running on CPU
    ecalMultiFitUncalibRecHitPh2
)

# process modifier to run the alpaka implementation
from Configuration.ProcessModifiers.alpaka_cff import alpaka

# ECAL multifit conditions on the device, filled from the five Phase-2
# records. Enclosed in a Task to prevent the construction of the ESProducer
# in the default (CPU) configuration.
from RecoLocalCalo.EcalRecProducers.ecalMultifitConditionsPhase2HostESProducer_cfi import ecalMultifitConditionsPhase2HostESProducer
ecalMultiFitUncalibRecHitPh2PortableConditions = cms.Task(ecalMultifitConditionsPhase2HostESProducer)

# The digi -> portable collection producer (simEcalUnsuppressedDigisSoA) is
# deliberately NOT defined or scheduled here; it is referenced by label only,
# mirroring the Run-3 multifit cff, which likewise consumes
# 'ecalDigisPortable:ebDigis' without owning that module. In the Phase-2
# workflows the module is provided by ecalUncalibRecHitPhase2_cff (weights)
# via ecalLocalRecoSequence. It cannot be defined here:
#  - an identical clone fails at process.load ("Trying to override definition
#    of simEcalUnsuppressedDigisSoA while it is used by the task
#    calolocalrecoTask") -- task-referenced labels are compared by object
#    identity, not parameter equality;
#  - importing the object from the weights cff is a circular import, because
#    on this branch ecalUncalibRecHitPhase2_cfi imports THIS cff to define
#    its CPU placeholder.
# A standalone (weights-free) multifit config must schedule the digi->SoA
# producer itself.

# ECAL Phase-2 multifit running on the accelerator
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2ProducerPortable_cfi import ecalUncalibRecHitPhase2ProducerPortable as _ecalUncalibRecHitPhase2ProducerPortable
ecalMultiFitUncalibRecHitPh2Portable = _ecalUncalibRecHitPhase2ProducerPortable.clone(
    digisLabelEB = 'simEcalUnsuppressedDigisSoA:ebDigis'
)

# convert the uncalibrated rechits from SoA to legacy format at the same
# module label consumed downstream
from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitSoAToLegacy_cfi import ecalUncalibRecHitSoAToLegacy as _ecalUncalibRecHitSoAToLegacy
alpaka.toReplaceWith(ecalMultiFitUncalibRecHitPh2, _ecalUncalibRecHitSoAToLegacy.clone(
    isPhase2 = True,
    inputCollectionEB = 'ecalMultiFitUncalibRecHitPh2Portable:EcalUncalibRecHitsEB',
    inputCollectionEE = None,
    outputLabelEE = None
))

alpaka.toReplaceWith(ecalMultiFitUncalibRecHitPh2Task, cms.Task(
    # ECAL multifit conditions on the device
    ecalMultiFitUncalibRecHitPh2PortableConditions,
    # ECAL Phase-2 multifit running on the device
    # (simEcalUnsuppressedDigisSoA is provided by the weights cff, see above)
    ecalMultiFitUncalibRecHitPh2Portable,
    # convert the uncalibrated rechits from SoA to legacy format
    ecalMultiFitUncalibRecHitPh2,
))

# for GPU validation run the CPU multifit alongside the alpaka modules
from Configuration.ProcessModifiers.gpuValidationEcal_cff import gpuValidationEcal
_ecalMultiFitUncalibRecHitPh2TaskValidation = ecalMultiFitUncalibRecHitPh2Task.copy()
_ecalMultiFitUncalibRecHitPh2TaskValidation.add(ecalMultiFitUncalibRecHitPh2Legacy)
gpuValidationEcal.toReplaceWith(ecalMultiFitUncalibRecHitPh2Task, _ecalMultiFitUncalibRecHitPh2TaskValidation)

# for alpaka validation compare alpaka serial with alpaka
from Configuration.ProcessModifiers.alpakaValidationEcal_cff import alpakaValidationEcal
from HeterogeneousCore.AlpakaCore.functions import makeSerialClone
ecalMultiFitUncalibRecHitPh2PortableSerialSync = makeSerialClone(ecalMultiFitUncalibRecHitPh2Portable)
ecalMultiFitUncalibRecHitPh2SerialSync = _ecalUncalibRecHitSoAToLegacy.clone(
    isPhase2 = True,
    inputCollectionEB = 'ecalMultiFitUncalibRecHitPh2PortableSerialSync:EcalUncalibRecHitsEB',
    inputCollectionEE = None,
    outputLabelEE = None
)
_ecalMultiFitUncalibRecHitPh2TaskValidation = ecalMultiFitUncalibRecHitPh2Task.copy()
_ecalMultiFitUncalibRecHitPh2TaskValidation.add(ecalMultiFitUncalibRecHitPh2PortableSerialSync)
_ecalMultiFitUncalibRecHitPh2TaskValidation.add(ecalMultiFitUncalibRecHitPh2SerialSync)
alpakaValidationEcal.toReplaceWith(ecalMultiFitUncalibRecHitPh2Task, _ecalMultiFitUncalibRecHitPh2TaskValidation)
