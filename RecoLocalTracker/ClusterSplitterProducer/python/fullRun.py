import FWCore.ParameterSet.Config as cms

from Geometry.CommonTopologies.globalTrackingGeometry_cfi import *
from Geometry.TrackerGeometryBuilder.trackerGeometry_cfi import *
from Configuration.ProcessModifiers.alpaka_cff import alpaka
process = cms.Process("RECOCC",alpaka)

process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')

process.load("RecoLocalTracker.SiPixelClusterizer.siPixelClustersPreSplitting_cff")
process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi")
process.load('RecoTracker.PixelTrackFitting.PixelTracks_cff')

process.load("RecoVertex.Configuration.RecoPixelVertexing_cff")
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.DQMFileSaver_cfi")
process.load("DQMServices.Components.DQMStoreStats_cfi")

process.dqmSaver.workflow = cms.untracked.string('/MyTest/JetCoreClusterSplitter/Timing')
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# HelperSplitter producer
process.candidateDataSoA = cms.EDProducer("HelperSplitter",
    Candidate = cms.InputTag("ak4CaloJets", "", "RECO"),
    siPixelClusters = cms.InputTag('siPixelClustersPreSplitting', '', 'RECO'),
    #siPixelClustersSoA=cms.InputTag("siPixelClustersPreSplittingAlpaka"),
    ptMin = cms.double(100),
    pixelCPE = cms.string("PixelCPEGeneric"),
    tanLorentzAngle = cms.double(0.0),
    tanLorentzAngleBarrelLayer1 = cms.double(0.0),
    verbose = cms.bool(False),
)

# trial producer (which uses the output from HelperSplitter)
process.trial = cms.EDProducer(
    "trial",
    #nHits=cms.uint32(100),
    #offset=cms.int32(10),
    ptMin=cms.double(100),
    deltaR=cms.double(0.05),
    chargeFracMin=cms.double(2.0),
    tanLorentzAngle=cms.double(0.0),
    tanLorentzAngleBarrelLayer1=cms.double(0.0),
    expSizeXAtLorentzAngleIncidence=cms.double(1.5),
    expSizeXDeltaPerTanAlpha=cms.double(0.0),
    expSizeYAtNormalIncidence=cms.double(1.3),
    centralMIPCharge=cms.double(26000),
    chargePerUnit=cms.double(2000),
    forceXError=cms.double(100),
    forceYError=cms.double(150),
    fractionalWidth=cms.double(0.4),
    #siPixelClusters=cms.InputTag("siPixelClustersPreSplittingAlpaka"),
    siPixelClusters = cms.InputTag("siPixelClustersPreSplittingAlpaka", "", "RECO"),
    siPixelDigis=cms.InputTag("candidateDataSoA"),              
    #siPixelDigis=cms.InputTag("siPixelClustersPreSplittingAlpaka"),              
    #siPixelDigis =cms.InputTag("siPixelClustersPreSplittingAlpaka", "", "RECO"), #assuming in the file
    ##trackingRecHits = cms.InputTag("siPixelRecHitsPreSplittingAlpaka"),
    #trackingRecHits = cms.InputTag("siPixelRecHitsPreSplittingAlpaka", "", "RECO"),
    candidateInput=cms.InputTag("candidateDataSoA"),
    #zVertex=cms.InputTag("pixelVerticesAlpaka"),
    geometryInput=cms.InputTag("candidateDataSoA"),
    verbose=cms.bool(False),
    debugMode = cms.bool(False),             #is True, only one cluster will be analyzed
    targetDetId = cms.int32(304152592),
    targetClusterOffset = cms.int32(0),
    targetEvent = cms.int32(1),    
    vertices = cms.InputTag('offlinePrimaryVertices'),
)

process.FastTimerService = cms.Service("FastTimerService",
    printEventSummary        = cms.untracked.bool(True),  # Print summary at the end
    printRunSummary          = cms.untracked.bool(False),
    printJobSummary          = cms.untracked.bool(True),  # Print total job performance
    enableDQM                = cms.untracked.bool(True),  # Enable DQM monitoring
    enableDQMbyModule        = cms.untracked.bool(True),  # Track time per module
    enableDQMbyPathActive    = cms.untracked.bool(True),  # Time per active path
    enableDQMbyPathTotal     = cms.untracked.bool(True),  # Total time per path
    enableDQMbyProcesses     = cms.untracked.bool(False), # If using subprocesses
    writeJSONSummary = cms.untracked.bool(True),
    jsonFileName = cms.untracked.string('resources'),
)

# DQM File Saver (Saves monitoring histograms)
process.dqmSaver.workflow = cms.untracked.string('/JetCoreClusterSplitter/Reco/DQMTest')
process.dqmSaver.convention = cms.untracked.string('Offline')
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)
process.dqm_step = cms.Path(process.dqmSaver)  # DQM Step


process.offlineBeamSpotDevice_step = cms.Path(process.offlineBeamSpotDevice)
process.siPixelClustersPreSplitting_step = cms.Path(process.siPixelClustersPreSplittingAlpaka)
process.beamSpotProducer_step = cms.Path(process.offlineBeamSpotDevice)
process.HelperSplitter_step = cms.Path(process.candidateDataSoA)
process.siPixelRecHitsPreSplitting_step = cms.Path(process.siPixelRecHitsPreSplittingAlpaka)
process.pixelVertexing_step = cms.Path(process.recopixelvertexing)
process.reconstruction_step1 = cms.Path(process.reconstruction_pixelTrackingOnly)
process.trial_step = cms.Path(process.trial)

# Set the schedule so that HelperSplitter runs before trial
process.schedule = cms.Schedule(
    process.reconstruction_step1,  
    #process.siPixelClustersPreSplitting_step,
    #process.siPixelRecHitsPreSplitting_step,
    #process.beamSpotProducer_step,
    process.HelperSplitter_step,
    #process.pixelVertexing_step,  
    process.trial_step,
    process.dqm_step    
)

#-----------------------------------------------------------------
# Input and Output Configuration
#-----------------------------------------------------------------
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:step3my.root')
    fileNames = cms.untracked.vstring('file:largestep3.root')
    #fileNames = cms.untracked.vstring('file:/gpu_data/store/relval/CMSSW_15_0_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_142X_mcRun3_2025_realistic_v7_STD_2025_PU-v3/2580000/1c2caeef-e246-4b6d-bebc-4fb6df4f9bbd.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

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
