import FWCore.ParameterSet.Config as cms

clusterSplitterProducer = cms.EDProducer("HelperSplitter",
    Candidate = cms.InputTag("ak4CaloJets", "", "RECO"),
    siPixelClusters = cms.InputTag("siPixelClustersPreSplitting", "", "RECO"),
    ptMin = cms.double(70),
    pixelCPE = cms.string("PixelCPEGeneric"),
    tanLorentzAngle = cms.double(0.0),
    tanLorentzAngleBarrelLayer1 = cms.double(0.0),
    verbose = cms.bool(False),
)

clusterSplitterProducerTask = cms.Task(clusterSplitterProducer)
