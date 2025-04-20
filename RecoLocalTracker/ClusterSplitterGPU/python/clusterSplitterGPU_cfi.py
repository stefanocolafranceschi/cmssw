import FWCore.ParameterSet.Config as cms

clusterSplitterGPU = cms.EDProducer("trial",
    nHits = cms.uint32(100),
    offset = cms.int32(10),
    ptMin = cms.double(100),
    deltaR = cms.double(0.05),
    chargeFracMin = cms.double(2.0),
    tanLorentzAngle = cms.double(0.001),
    tanLorentzAngleBarrelLayer1 = cms.double(0.001),
    expSizeXAtLorentzAngleIncidence = cms.double(1.5),
    expSizeXDeltaPerTanAlpha = cms.double(0.0),
    expSizeYAtNormalIncidence = cms.double(1.3),
    centralMIPCharge = cms.double(26000),
    chargePerUnit = cms.double(2000),
    forceXError = cms.double(100),
    forceYError = cms.double(150),
    fractionalWidth = cms.double(0.4),
    siPixelClusters = cms.InputTag("siPixelClustersPreSplittingAlpaka", "", "RECO"),
    siPixelDigis = cms.InputTag("candidateDataSoA"),
    trackingRecHits = cms.InputTag("siPixelRecHitsPreSplittingAlpaka"),
    candidateInput = cms.InputTag("candidateDataSoA"),
    geometryInput = cms.InputTag("candidateDataSoA"),
    verbose = cms.bool(False),
    debugMode = cms.bool(False),
    targetDetId = cms.int32(304181256),
    targetClusterOffset = cms.int32(2),
    targetEvent = cms.int32(1),
    vertices = cms.InputTag("offlinePrimaryVertices")
)

clusterSplitterGPUTask = cms.Task(clusterSplitterGPU)
