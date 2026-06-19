#import  RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2WeightsProducer_cfi as _mod
#ecalUncalibRecHitPhase2 = _mod.ecalUncalibRecHitPhase2WeightsProducer.clone()

from SimCalorimetry.EcalSimProducers.esCATIAGainProducer_cfi import *
from SimCalorimetry.EcalSimProducers.esEcalLiteDTUPedestalsProducer_cfi import *

#from RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHitPh2_cfi import *
#ecalUncalibRecHitPhase2 = ecalFixedAlphaBetaFitUncalibRecHitPh2.clone(
#    EBdigiCollection = cms.InputTag("simEcalUnsuppressedDigis")
#)

# use only the cpu branch of the switchProducer for now
from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHitPh2_cff import *
ecalUncalibRecHitPhase2 = ecalMultiFitUncalibRecHitPh2.cpu.clone(
    EBdigiCollection = cms.InputTag("simEcalUnsuppressedDigis")
)

