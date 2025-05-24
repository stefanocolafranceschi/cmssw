import FWCore.ParameterSet.Config as cms

def customize_step3(process):
    # -------------------------------------------------------
    # Load producers
    process.load("RecoLocalTracker.ClusterSplitterProducer.clusterSplitterProducer_cfi")
    process.candidateDataSoA = process.clusterSplitterProducer.clone()

    process.load("RecoLocalTracker.ClusterSplitterGPU.clusterSplitterGPU_cfi")
    process.trial = process.clusterSplitterGPU.clone(
        candidateInput = cms.InputTag("candidateDataSoA"),
        geometryInput = cms.InputTag("candidateDataSoA"),
        siPixelDigis = cms.InputTag("candidateDataSoA"),
        clusterPixelCounts = cms.InputTag("candidateDataSoA")
    )

    if not hasattr(process, "clusterSplitterSequence"):
        process.clusterSplitterSequence = cms.Sequence(
            process.candidateDataSoA +
            process.trial
        )
        process.pathClusterSplitter = cms.Path(process.clusterSplitterSequence)
        process.schedule.append(process.pathClusterSplitter)

    # -------------------------------------------------------
    # FastTimerService with JSON and DQM
    process.load("DQMServices.Core.DQMStore_cfi")
    process.load("DQMServices.Components.DQMFileSaver_cfi")
    process.load("DQMServices.Components.DQMStoreStats_cfi")

    process.dqmSaver.workflow = cms.untracked.string('/JetCoreClusterSplitter/Reco/DQMTest')
    process.dqmSaver.convention = cms.untracked.string('Offline')
    process.dqmSaver.saveByRun = cms.untracked.int32(-1)
    process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True)

    process.FastTimerService = cms.Service("FastTimerService",
        printEventSummary        = cms.untracked.bool(True),
        printRunSummary          = cms.untracked.bool(False),
        printJobSummary          = cms.untracked.bool(True),
        enableDQM                = cms.untracked.bool(True),
        enableDQMbyModule        = cms.untracked.bool(True),
        enableDQMbyPathActive    = cms.untracked.bool(True),
        enableDQMbyPathTotal     = cms.untracked.bool(True),
        enableDQMbyProcesses     = cms.untracked.bool(False),
        writeJSONSummary         = cms.untracked.bool(True),
        jsonFileName             = cms.untracked.string('step3_timing.json'),
    )

#    process.options = cms.untracked.PSet(
#        IgnoreCompletely = cms.untracked.vstring(),
#        Rethrow = cms.untracked.vstring(),
#        TryToContinue = cms.untracked.vstring(),
#        accelerators = cms.untracked.vstring('cpu'),
#    )


    if not hasattr(process, "dqm_step"):
        process.dqm_step = cms.Path(process.dqmSaver)
        process.schedule.append(process.dqm_step)

    return process
