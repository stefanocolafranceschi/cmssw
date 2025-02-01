import FWCore.ParameterSet.Config as cms

trialProducer = cms.EDProducer(
    "trial",
    Candidate = cms.InputTag("Candidate"),
    SiPixelClusters = cms.InputTag("SiPixelClusters"),
    ptMin = cms.double(0.5),            # Default value
    tanLorentzAngle = cms.double(0.1),  # Default value
    tanLorentzAngleBarrelLayer1 = cms.double(0.2),  # Default value
    verbose = cms.bool(False)
)
