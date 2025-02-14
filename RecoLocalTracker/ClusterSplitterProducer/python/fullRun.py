import FWCore.ParameterSet.Config as cms

#-----------------------------------------------------------------
# Process and Geometry Setup
#-----------------------------------------------------------------
# Use a unique process name (here "RECOCC")
process = cms.Process("RECOCC")

# Load geometry from the DB, which also produces the GlobalTrackingGeometryRecord,
# TrackerTopologyRcd, and related records.
process.load("Configuration.Geometry.GeometryRecoDB_cff")
# Optionally, load the global tracking geometry if needed
process.load("Geometry.CommonTopologies.globalTrackingGeometry_cfi")

# Load the GlobalTag configuration
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# Adjust the GlobalTag as needed for your release/conditions
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2024_realistic', '')

#-----------------------------------------------------------------
# Define Your Custom Producers
#-----------------------------------------------------------------
# HelperSplitter producer
process.HelperSplitter = cms.EDProducer("HelperSplitter",
    Candidate = cms.InputTag("ak4CaloJets", "", "RECO"),
    SiPixelClusters = cms.InputTag("SiPixelClusters", "", "RECO"),
    ptMin = cms.double(0.5),
    tanLorentzAngle = cms.double(0.1),
    tanLorentzAngleBarrelLayer1 = cms.double(0.2),
    verbose = cms.bool(True)
)

# trial producer (which uses the output from HelperSplitter)
process.trial = cms.EDProducer("trial",
    nHits = cms.uint32(100),
    offset = cms.int32(10),
    ptMin = cms.double(200),
    deltaR = cms.double(0.05),
    chargeFracMin = cms.double(2.0),
    tanLorentzAngle = cms.double(0.02),
    tanLorentzAngleBarrelLayer1 = cms.double(0.015),
    expSizeXAtLorentzAngleIncidence = cms.double(0.1),
    expSizeXDeltaPerTanAlpha = cms.double(0.02),
    expSizeYAtNormalIncidence = cms.double(0.1),
    centralMIPCharge = cms.double(26000),
    chargePerUnit = cms.double(2000),
    forceXError = cms.double(100),
    forceYError = cms.double(150),
    fractionalWidth = cms.double(0.4),
    siPixelClusters = cms.InputTag("SiPixelClustersSoACollection"),
    siPixelDigis = cms.InputTag("SiPixelDigisSoACollection"),
    trackingRecHits = cms.InputTag("trackingRecHitsSoACollection"),
    candidateInput = cms.InputTag("candidateDataSoA"),
    zVertex = cms.InputTag("zVertex"),
    geometryInput = cms.InputTag("ClusterGeometrySoA"),
    verbose = cms.bool(False)
)

#-----------------------------------------------------------------
# Define Execution Paths and Schedule
#-----------------------------------------------------------------
# Here we define one Path for each producer
process.HelperSplitter_step = cms.Path(process.HelperSplitter)
process.trial_step = cms.Path(process.trial)

# Set the schedule so that HelperSplitter runs before trial
process.schedule = cms.Schedule(
    process.HelperSplitter_step,
    process.trial_step
)

#-----------------------------------------------------------------
# Input and Output Configuration
#-----------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:step_output.root'),
    outputCommands = cms.untracked.vstring("keep *_*_*_*")
)
process.endpath = cms.EndPath(process.output)
process.schedule.append(process.endpath)