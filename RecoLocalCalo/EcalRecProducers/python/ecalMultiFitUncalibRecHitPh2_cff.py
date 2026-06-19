import FWCore.ParameterSet.Config as cms
from HeterogeneousCore.CUDACore.SwitchProducerCUDA import SwitchProducerCUDA

# provides the Phase 2 SampleMask and the Phase 2 SamplesCorrelation
# all other payloades are not produced
from CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi import EcalTrivialConditionRetriever as _EcalTrivialConditionRetriever
EcalTrivialConditionRetriever = _EcalTrivialConditionRetriever.clone(
  producedEcalAlignmentEB = cms.untracked.bool(False),
  producedEcalAlignmentEE = cms.untracked.bool(False),
  producedEcalAlignmentES = cms.untracked.bool(False),
  producedEcalSimPulseShape = cms.untracked.bool(False),
  producedEcalChannelStatus = cms.untracked.bool(False),
  producedEcalDQMTowerStatus = cms.untracked.bool(False),
  producedEcalDQMChannelStatus = cms.untracked.bool(False),
  producedEcalDCSTowerStatus = cms.untracked.bool(False),
  producedEcalDAQTowerStatus = cms.untracked.bool(False),
  producedEcalTrgChannelStatus = cms.untracked.bool(False),
  producedEcalPedestals = cms.untracked.bool(False),
  producedEcalWeights = cms.untracked.bool(False),
  producedEcalLinearCorrections = cms.untracked.bool(False),
  producedEcalIntercalibConstants = cms.untracked.bool(False),
  producedEcalIntercalibConstantsMC = cms.untracked.bool(False),
  producedEcalIntercalibErrors = cms.untracked.bool(False),
  producedEcalTimeCalibConstants = cms.untracked.bool(False),
  producedEcalTimeCalibErrors = cms.untracked.bool(False),
  producedEcalTimeOffsetConstant = cms.untracked.bool(False),
  producedEcalLaserCorrection = cms.untracked.bool(False),
  producedEcalGainRatios = cms.untracked.bool(False),
  producedEcalADCToGeVConstant = cms.untracked.bool(False),
  producedEcalMappingElectronics = cms.untracked.bool(False),
  producedEcalClusterLocalContCorrParameters = cms.untracked.bool(False),
  producedEcalClusterCrackCorrParameters = cms.untracked.bool(False),
  producedEcalClusterEnergyUncertaintyParameters = cms.untracked.bool(False),
  producedEcalClusterEnergyCorrectionParameters = cms.untracked.bool(False),
  producedEcalClusterEnergyCorrectionObjectSpecificParameters = cms.untracked.bool(False),
  producedEcalSampleMask = cms.untracked.bool(False),
  producedEcalTimeBiasCorrections = cms.untracked.bool(False),
  producedEcalSamplesCorrelation = cms.untracked.bool(False),

  getWeightsFromFile = cms.untracked.bool(False),
  producedEcalPh2Weights = cms.untracked.bool(True),
  producedEcalPh2SampleMask = cms.untracked.bool(True),
  producedEcalPh2SamplesCorrelation = cms.untracked.bool(True)
)

# provides the Phase 2 pulse shape (same for all channels)
ecalPh2PulseShapesRcd = cms.ESSource("EmptyESSource",
    recordName = cms.string("EcalPh2PulseShapesRcd"),
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True)
)
ecalPh2PulseShapesESProducer = cms.ESProducer("EcalPh2PulseShapesESProducer",
    pulseShapes = cms.vdouble(0., 0., 0., 0.0113979, 0.758151, 1.0, 0.887744, 0.673548, 0.474332, 0.319561, 0.215144, 0.147464, 0.101087, 0.0693181, 0.0475044, 0.)
)

# provides the Phase 2 pulse covariances (same for all channels)
ecalPh2PulseCovariancesRcd = cms.ESSource("EmptyESSource",
    recordName = cms.string("EcalPh2PulseCovariancesRcd"),
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True)
)
ecalPh2PulseCovariancesESProducer = cms.ESProducer("EcalPh2PulseCovariancesESProducer",
    pulseCovariances = cms.vdouble(
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  3.001e-06,  1.233e-05,  0.000e+00, -4.416e-06, -4.571e-06, -3.614e-06, -2.636e-06, -1.286e-06, -8.410e-07, -5.296e-07,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  1.233e-05,  6.154e-05,  0.000e+00, -2.200e-05, -2.309e-05, -1.838e-05, -1.373e-05, -7.334e-06, -5.088e-06, -3.745e-06, -2.428e-06,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -4.416e-06, -2.200e-05,  0.000e+00,  8.319e-06,  8.545e-06,  6.792e-06,  5.059e-06,  2.678e-06,  1.816e-06,  1.223e-06,  8.245e-07,  5.589e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -4.571e-06, -2.309e-05,  0.000e+00,  8.545e-06,  9.182e-06,  7.219e-06,  5.388e-06,  2.853e-06,  1.944e-06,  1.324e-06,  9.083e-07,  6.335e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -3.614e-06, -1.838e-05,  0.000e+00,  6.792e-06,  7.219e-06,  6.016e-06,  4.437e-06,  2.385e-06,  1.636e-06,  1.118e-06,  7.754e-07,  5.556e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -2.636e-06, -1.373e-05,  0.000e+00,  5.059e-06,  5.388e-06,  4.437e-06,  3.602e-06,  1.917e-06,  1.322e-06,  9.079e-07,  6.529e-07,  4.752e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -1.286e-06, -7.334e-06,  0.000e+00,  2.678e-06,  2.853e-06,  2.385e-06,  1.917e-06,  1.375e-06,  9.100e-07,  6.455e-07,  4.693e-07,  3.657e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -8.410e-07, -5.088e-06,  0.000e+00,  1.816e-06,  1.944e-06,  1.636e-06,  1.322e-06,  9.100e-07,  9.115e-07,  6.062e-07,  4.436e-07,  3.422e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00, -5.296e-07, -3.745e-06,  0.000e+00,  1.223e-06,  1.324e-06,  1.118e-06,  9.079e-07,  6.455e-07,  6.062e-07,  7.217e-07,  4.862e-07,  3.768e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00, -2.428e-06,  0.000e+00,  8.245e-07,  9.083e-07,  7.754e-07,  6.529e-07,  4.693e-07,  4.436e-07,  4.862e-07,  6.509e-07,  4.418e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  5.589e-07,  6.335e-07,  5.556e-07,  4.752e-07,  3.657e-07,  3.422e-07,  3.768e-07,  4.418e-07,  6.142e-07,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
        0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,  0.000e+00,
    )
)

# ECAL multifit running on CPU
from RecoLocalCalo.EcalRecProducers.ecalMultiFitUncalibRecHitPh2_cfi import ecalMultiFitUncalibRecHitPh2 as _ecalMultiFitUncalibRecHitPh2
ecalMultiFitUncalibRecHitPh2 = SwitchProducerCUDA(
  cpu = _ecalMultiFitUncalibRecHitPh2.clone()
)

ecalMultiFitUncalibRecHitPh2Task = cms.Task(
  # ECAL multifit running on CPU
  ecalMultiFitUncalibRecHitPh2
)

