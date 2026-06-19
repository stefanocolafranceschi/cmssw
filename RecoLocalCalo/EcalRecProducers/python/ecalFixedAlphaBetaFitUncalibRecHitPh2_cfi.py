import FWCore.ParameterSet.Config as cms

# producer of rechits starting from digis
ecalFixedAlphaBetaFitUncalibRecHitPh2 = cms.EDProducer("EcalUncalibRecHitProducer",
    EEdigiCollection = cms.InputTag(""),
    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    EEhitCollection = cms.string(""),
    EBhitCollection = cms.string("EcalUncalibRecHitsEB"),
    IsPhase2 = cms.bool(True),
    algo = cms.string("EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2"),
    algoPSet = cms.PSet(
      alpha = cms.double(1.138),
      beta = cms.double(1.655),
      alphaBetaFilename = cms.untracked.string(""),
      minAmpl = cms.double(8.0),
      useDynamicPedestal = cms.bool(True)
    )
)
