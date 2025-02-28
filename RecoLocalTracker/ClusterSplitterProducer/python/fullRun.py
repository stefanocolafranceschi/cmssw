import FWCore.ParameterSet.Config as cms

from Geometry.CommonTopologies.globalTrackingGeometry_cfi import *
from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
from Configuration.ProcessModifiers.alpaka_cff import alpaka
process = cms.Process("RECOCC",alpaka)

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoLocalTracker.SiPixelClusterizer.siPixelClustersPreSplitting_cff")
process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
#process.load("RecoTracker.Configuration.RecoPixelVertexing_cff")
process.load('Configuration.EventContent.EventContent_cff')
#process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')
process.load('RecoTracker.PixelTrackFitting.PixelTracks_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# HelperSplitter producer
process.candidateDataSoA = cms.EDProducer("HelperSplitter",
    Candidate = cms.InputTag("ak4CaloJets", "", "RECO"),
    siPixelClusters = cms.InputTag("siPixelClusters","","RECO"),
    ptMin = cms.double(50),
    tanLorentzAngle = cms.double(0.1),
    tanLorentzAngleBarrelLayer1 = cms.double(0.2),
    verbose = cms.bool(True)
)

# trial producer (which uses the output from HelperSplitter)
process.trial = cms.EDProducer(
    "trial",
    nHits=cms.uint32(100),
    offset=cms.int32(10),
    ptMin=cms.double(2),
    deltaR=cms.double(100.05),
    chargeFracMin=cms.double(2.0),
    tanLorentzAngle=cms.double(0.02),
    tanLorentzAngleBarrelLayer1=cms.double(0.015),
    expSizeXAtLorentzAngleIncidence=cms.double(0.1),
    expSizeXDeltaPerTanAlpha=cms.double(0.02),
    expSizeYAtNormalIncidence=cms.double(0.1),
    centralMIPCharge=cms.double(26000),
    chargePerUnit=cms.double(2000),
    forceXError=cms.double(100),
    forceYError=cms.double(150),
    fractionalWidth=cms.double(0.4),
    siPixelClusters=cms.InputTag("siPixelClustersPreSplittingAlpaka"),
    siPixelDigis=cms.InputTag("siPixelClustersPreSplittingAlpaka"),
    trackingRecHits = cms.InputTag("siPixelRecHitsPreSplittingAlpaka"),
    candidateInput=cms.InputTag("candidateDataSoA"),
    zVertex=cms.InputTag("pixelVerticesAlpaka"),
    geometryInput=cms.InputTag("candidateDataSoA"),
    verbose=cms.bool(True),
)

process.HelperSplitter_step = cms.Path(process.candidateDataSoA)
process.offlineBeamSpotDevice_step = cms.Path(process.offlineBeamSpotDevice)
process.siPixelClustersPreSplitting_step = cms.Path(process.siPixelClustersPreSplittingAlpaka)
process.siPixelRecHitsPreSplitting_step = cms.Path(process.siPixelRecHitsPreSplittingAlpaka)
#process.pixelVertexing_step = cms.Path(process.recopixelvertexing)
process.reconstruction_step1 = cms.Path(process.reconstruction_pixelTrackingOnly)
process.trial_step = cms.Path(process.trial)

# Set the schedule so that HelperSplitter runs before trial
process.schedule = cms.Schedule(
    process.HelperSplitter_step,
    process.siPixelClustersPreSplitting_step,
    process.siPixelRecHitsPreSplitting_step,
    #process.pixelVertexing_step,  
    process.reconstruction_step1,  
    process.trial_step
)

#-----------------------------------------------------------------
# Input and Output Configuration
#-----------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

#process.output = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('file:step_output.root'),
#    outputCommands = cms.untracked.vstring("keep *_*_*_*")
#)

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step4.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.endpath = cms.EndPath(process.RECOSIMoutput)
process.schedule.append(process.endpath)
