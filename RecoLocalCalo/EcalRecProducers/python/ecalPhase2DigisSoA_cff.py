import FWCore.ParameterSet.Config as cms

# The ECAL Phase-2 digi -> portable (SoA) conversion, in its own leaf cff so
# that BOTH local-reco cffs (weights: ecalUncalibRecHitPhase2_cff, multifit:
# ecalMultiFitUncalibRecHitPh2_cff) can import the SAME module OBJECT.
#
# Why object sharing matters: process.load compares task-referenced labels by
# object identity, not parameter equality -- two cffs may only define the
# same label if they share the Python object (step-8 lesson 1). A leaf cff
# with no further reco imports also keeps the import graph acyclic
# (step-8 lesson 2).

from RecoLocalCalo.EcalRecProducers.ecalPhase2DigiToPortableProducer_cfi import ecalPhase2DigiToPortableProducer as _ecalPhase2DigiToPortableProducer

simEcalUnsuppressedDigisSoA = _ecalPhase2DigiToPortableProducer.clone()
