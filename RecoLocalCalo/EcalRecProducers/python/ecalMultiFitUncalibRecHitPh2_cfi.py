import FWCore.ParameterSet.Config as cms

ecalMultiFitUncalibRecHitPh2 = cms.EDProducer("EcalUncalibRecHitProducer",
    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
    EEdigiCollection = cms.InputTag(""),
    EBhitCollection = cms.string("EcalUncalibRecHitsEB"),
    EEhitCollection = cms.string(""),
    IsPhase2 = cms.bool(True),
    algo = cms.string("EcalUncalibRecHitWorkerMultiFitPh2"),
    algoPSet = cms.PSet(
      # for multifit method
      activeBXs = cms.vint32(-2,-1,0,1,2),
      ampErrorCalculation = cms.bool(True),
      useLumiInfoRunHeader = cms.bool(True),
  
      doPrefit = cms.bool(False),
      prefitMaxChiSq = cms.double(25.),
      
      dynamicPedestals = cms.bool(False),
      mitigateBadSamples = cms.bool(False),
      gainSwitchUseMaxSample = cms.bool(True),
      selectiveBadSampleCriteria = cms.bool(False),
      simplifiedNoiseModelForGainSwitch = cms.bool(True),
      addPedestalUncertainty = cms.double(0.),
  
      # decide which algorithm to be use to calculate the jitter
      timealgo = cms.string("RatioMethod"),
  
      # for ratio method
      timeFitParameters = cms.vdouble(-2.015452e+00, 3.130702e+00, -1.234730e+01, 4.188921e+01, -8.283944e+01, 9.101147e+01, -5.035761e+01, 1.105621e+01),
      amplitudeFitParameters = cms.vdouble(1.138, 1.652),
      timeFitLimits_Lower = cms.double(0.2),
      timeFitLimits_Upper = cms.double(1.4),
      # for time error
      timeConstantTerm= cms.double(.6),
      # for kOutOfTime flag
      timeNconst = cms.double(28.5),
      outOfTimeThresholdGain10p = cms.double(5),
      outOfTimeThresholdGain10m = cms.double(5),
      outOfTimeThresholdGain1p  = cms.double(5),
      outOfTimeThresholdGain1m  = cms.double(5),
      amplitudeThreshold = cms.double(10),
   )                                           
)
