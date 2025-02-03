import FWCore.ParameterSet.Config as cms

def customizeForAlpakaSplitter(process):

    # Load necessary producers
    process.load("RecoLocalCalo.CaloTowersCreator.calotowermaker_cfi")  # Load Calo Towers producer
    process.load("RecoJets.JetProducers.ak4CaloJets_cfi")  # Load the ak4CaloJets producer

    # Define a Task to ensure dependencies are managed explicitly
    process.towerMakerTask = cms.Task(process.calotowermaker)
    process.ak4CaloJetsTask = cms.Task(process.ak4CaloJets)

    # Add Tasks to Paths
    process.towerMaker_step = cms.Path(process.towerMakerTask)
    process.ak4CaloJets_step = cms.Path(process.ak4CaloJetsTask)

    # Define the HelperSplitter producer
    process.HelperSplitter = cms.EDProducer("HelperSplitter",
        Candidate = cms.InputTag("ak4CaloJets"),          # Input for HelperSplitter
        SiPixelClusters = cms.InputTag("SiPixelClusters"),
        ptMin = cms.double(0.5),                          # Default value
        tanLorentzAngle = cms.double(0.1),                # Default value
        tanLorentzAngleBarrelLayer1 = cms.double(0.2),    # Default value
        verbose = cms.bool(False)
    )

    # Create a Task for HelperSplitter
    process.HelperSplitterTask = cms.Task(process.HelperSplitter)
    process.HelperSplitter_step = cms.Path(process.HelperSplitterTask)

    # Enforce execution order in the schedule
    process.schedule.extend([process.towerMaker_step, process.ak4CaloJets_step, process.HelperSplitter_step])

    return process
