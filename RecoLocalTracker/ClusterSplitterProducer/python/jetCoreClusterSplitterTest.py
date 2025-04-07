import FWCore.ParameterSet.Config as cms

process = cms.Process("RECOOOOOO")

# Standard services, geometry, magnetic field, and GlobalTag
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoLocalTracker.SiStripClusterizer.SiStripClusterChargeCut_cfi")
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi")
process.load("HLTrigger.Timer.FastTimerService_cfi")
process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.DQMFileSaver_cfi")
process.load("DQMServices.Components.DQMStoreStats_cfi")

process.dqmSaver.workflow = cms.untracked.string('/MyTest/JetCoreClusterSplitter/Timing')
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

# Load Pixel Clusterizer
process.load("RecoLocalTracker.SiPixelClusterizer.siPixelClustersPreSplitting_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

# Define the JetCoreClusterSplitter EDProducer
process.jetCoreClusterSplitterTest = cms.EDProducer("JetCoreClusterSplitter",
    pixelClusters = cms.InputTag('siPixelClustersPreSplitting', '', 'RECO'),
    vertices              = cms.InputTag('offlinePrimaryVertices'),
    pixelCPE              = cms.string("PixelCPEGeneric"),
    verbose               = cms.bool(False),
    debugMode             = cms.bool(False),         #is True, only one cluster will be analyzed
    targetDetId           = cms.int32(304181256),
    targetClusterOffset   = cms.int32(2),
    targetEvent           = cms.int32(1),
    ptMin                 = cms.double(100),
    cores                 = cms.InputTag("ak4CaloJets", "", "RECO"),
    chargeFractionMin     = cms.double(2.0),
    deltaRmax             = cms.double(0.05),
    forceXError           = cms.double(100),
    forceYError           = cms.double(150),
    fractionalWidth       = cms.double(0.4),
    chargePerUnit         = cms.double(2000),
    centralMIPCharge      = cms.double(26000)
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



# Define the process path
process.raw2digi_step = cms.Path(process.RawToDigi_pixelOnly)
process.siPixelClustersPreSplitting_step = cms.Path(process.siPixelClustersPreSplitting)
process.jetCoreClusterSplitter_step = cms.Path(process.jetCoreClusterSplitterTest)
process.reconstruction_step1 = cms.Path(process.reconstruction_pixelTrackingOnly)
process.dqm_step = cms.Path(process.dqmSaver)  # DQM Step

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:step3my.root')
    fileNames = cms.untracked.vstring('file:largestep3.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(999))

# Output module
process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier   = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:jetCoreClusterSplitter_output.root'),
    outputCommands = cms.untracked.vstring("keep *_*_*_*"),
    splitLevel = cms.untracked.int32(0)
)
process.out = cms.EndPath(process.RECOSIMoutput)


# Set the schedule
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.siPixelClustersPreSplitting_step,
    #process.reconstruction_step1,
    process.jetCoreClusterSplitter_step,
    process.out,
    process.dqm_step
)


