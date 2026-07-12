import FWCore.ParameterSet.Config as cms

# CPU definition of the Phase-2 uncalibrated-rechit producer label.
# Under the alpaka process modifier this module is replaced by the SoA->legacy
# converter fed by the device weights producer (see ecalUncalibRecHitPhase2_cff),
# so this is effectively a placeholder for non-alpaka configurations.
#
# It clones the CPU multifit *cfi* (a leaf, with its own defaults). Earlier
# this file imported the multifit *cff* -- a cfi importing a cff is inverted
# layering and created a circular-import chain
# (weights cff -> this cfi -> multifit cff -> ...); do not reintroduce it.

from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHitPh2_cfi import ecalMultiFitUncalibRecHitPh2 as _ecalMultiFitUncalibRecHitPh2

ecalUncalibRecHitPhase2 = _ecalMultiFitUncalibRecHitPh2.clone(
    EBdigiCollection = cms.InputTag("simEcalUnsuppressedDigis")
)
