import FWCore.ParameterSet.Config as cms

# In the alpaka-wired Phase-2 reco (ecalUncalibRecHitPhase2_cff with the alpaka
# process modifier):
#   simEcalUnsuppressedDigisSoA   : digi -> portable SoA producer
#   ecalUncalibRecHitPhase2SoA    : weights producer on device
#   ecalUncalibRecHitPhase2       : SoA -> legacy converter (consumed downstream)

_DIGIS = "simEcalUnsuppressedDigisSoA:ebDigis"

# Measured from step2.root digis with dump_phase2_pulse_template.py
# (21 signal frames, 116069 noise frames; modal peak sample 5 = kMaxSampleIdx;
# true-pedestal (13.0) subtraction, no per-frame mean).
# 16-sample pulse template, template[0] at digi sample 4 (peak = template[1] at
# sample 5); small negative tail values are measurement fluctuations around 0.
_PULSE_SHAPE = [0.646882, 1.0, 0.763388, 0.384452, 0.145384, 0.0385834,
                0.0115881, -0.00364541, -0.0147806, -0.00316459, -0.00744192,
                -0.00618433, 0.0, 0.0, 0.0, 0.0]

# noise sample correlation vs |i-j| from pedestal-only frames; measured in
# gain 10, reused for gain 1 (not separately measurable from these digis).
# Toeplitz matrix is positive definite (min eigenvalue 0.154, no shrinkage).
_SAMPLE_CORR = [1.0, 0.714871, 0.572326, 0.4849, 0.434294, 0.394461,
                0.377368, 0.364767, 0.356238, 0.343968, 0.43524, 0.434521,
                0.434619, 0.436439, 0.437694, 0.43699]

# pulse-shape covariance (row-major 16x16, template indexing) from
# normalized-pulse residuals; row/col 1 vanish by peak normalization,
# rows/cols 12-15 are the extrapolated tail (no data)
_PULSE_COV = [
    0.000429833, 0, -8.96947e-05, -0.000104676, 0.000234094, 0.000217236, 0.000213305, 8.44262e-05, 2.53562e-05, 0.000183617, 8.54719e-05, -3.95251e-06, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -8.96947e-05, 0, 0.00043817, 4.46164e-05, -8.2173e-05, -8.09094e-05, -1.82909e-06, -1.49307e-05, -4.19575e-05, -0.000125516, -0.000180117, -6.32121e-05, 0, 0, 0, 0,
    -0.000104676, 0, 4.46164e-05, 0.000198505, 8.41339e-05, -2.51144e-06, -6.37521e-05, 3.86776e-05, 4.40365e-05, -4.26698e-05, -9.86186e-05, -1.06118e-06, 0, 0, 0, 0,
    0.000234094, 0, -8.2173e-05, 8.41339e-05, 0.000494388, 0.000347399, 0.000165319, 0.00017497, 0.000270604, 0.000278583, 5.26999e-05, 6.3043e-05, 0, 0, 0, 0,
    0.000217236, 0, -8.09094e-05, -2.51144e-06, 0.000347399, 0.000585171, 0.000452957, 0.000433725, 0.000493059, 0.000528928, 0.000151409, 0.000243072, 0, 0, 0, 0,
    0.000213305, 0, -1.82909e-06, -6.37521e-05, 0.000165319, 0.000452957, 0.000620015, 0.000382771, 0.000278701, 0.000414631, 0.000144976, 0.000182151, 0, 0, 0, 0,
    8.44262e-05, 0, -1.49307e-05, 3.86776e-05, 0.00017497, 0.000433725, 0.000382771, 0.000632292, 0.000562114, 0.000415537, 0.00010883, 0.000259739, 0, 0, 0, 0,
    2.53562e-05, 0, -4.19575e-05, 4.40365e-05, 0.000270604, 0.000493059, 0.000278701, 0.000562114, 0.000872625, 0.000574681, 0.000258053, 0.000348987, 0, 0, 0, 0,
    0.000183617, 0, -0.000125516, -4.26698e-05, 0.000278583, 0.000528928, 0.000414631, 0.000415537, 0.000574681, 0.000765147, 0.000265969, 0.000206068, 0, 0, 0, 0,
    8.54719e-05, 0, -0.000180117, -9.86186e-05, 5.26999e-05, 0.000151409, 0.000144976, 0.00010883, 0.000258053, 0.000265969, 0.000428097, 0.000193513, 0, 0, 0, 0,
    -3.95251e-06, 0, -6.32121e-05, -1.06118e-06, 6.3043e-05, 0.000243072, 0.000182151, 0.000259739, 0.000348987, 0.000206068, 0.000193513, 0.000319126, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]

_PHASE2_RECORDS = (
    "EcalLiteDTUPedestalsRcd",
    "EcalCATIAGainRatiosRcd",
    "EcalPh2PulseShapesRcd",
    "EcalPh2PulseCovariancesRcd",
    "EcalPh2SamplesCorrelationRcd",
)


def _loadPhase2MultifitModules(process):
    # No tags exist in the conditions DB for the five Phase-2 records (checked
    # 150X_mcRun4_realistic_v1, Jul 2026): serve them from the trivial
    # ESProducer, with EmptyESSources providing the IOVs.
    from RecoLocalCalo.EcalRecProducers.ecalPhase2TrivialCondESProducer_cfi import (
        ecalPhase2TrivialCondESProducer as _trivialCond,
    )
    process.ecalPhase2TrivialCondESProducer = _trivialCond.clone(
        pulseShape=_PULSE_SHAPE,
        sampleCorrelationG10=_SAMPLE_CORR,
        sampleCorrelationG1=_SAMPLE_CORR,
        pulseCovariance=_PULSE_COV,
    )
    for rec in _PHASE2_RECORDS:
        setattr(
            process,
            "empty" + rec + "Source",
            cms.ESSource(
                "EmptyESSource",
                recordName=cms.string(rec),
                firstValid=cms.vuint32(1),
                iovIsRunNotTime=cms.bool(True),
            ),
        )

    # ES producer filling EcalMultifitConditionsPhase2Rcd from the
    # EcalPh2* / LiteDTU / CATIA condition records
    process.load("RecoLocalCalo.EcalRecProducers.ecalMultifitConditionsPhase2HostESProducer_cfi")


def customisePhase2MultifitGPU(process):
    """Replace the Phase-2 weights DEVICE producer with the alpaka multifit
    producer at the same module label (ecalUncalibRecHitPhase2SoA). The
    SoA->legacy converter and everything downstream are untouched and consume
    the multifit result. Compare against a stock (weights) job to validate."""
    if not hasattr(process, "ecalUncalibRecHitPhase2SoA"):
        raise RuntimeError(
            "customisePhase2MultifitGPU: no 'ecalUncalibRecHitPhase2SoA' module in the process; "
            "this needs the alpaka-wired Phase-2 reco (ecalDevelAlpaka workflow)")

    _loadPhase2MultifitModules(process)

    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2ProducerPortable_cfi import (
        ecalUncalibRecHitPhase2ProducerPortable as _multifit,
    )
    process.globalReplace("ecalUncalibRecHitPhase2SoA", _multifit.clone(digisLabelEB=_DIGIS))
    return process


def customisePhase2MultifitGPUAlongside(process):
    """Run the alpaka multifit producer in addition to the weights producer,
    convert its output to the legacy collection format, and keep it in the
    event output for a same-file comparison:
      weights : ecalUncalibRecHitPhase2 : EcalUncalibRecHitsEB
      multifit: ecalMultiFitUncalibRecHitPhase2 : EcalUncalibRecHitsEB"""
    if not hasattr(process, "simEcalUnsuppressedDigisSoA"):
        raise RuntimeError(
            "customisePhase2MultifitGPUAlongside: no 'simEcalUnsuppressedDigisSoA' module in the "
            "process; this needs the alpaka-wired Phase-2 reco (ecalDevelAlpaka workflow)")

    _loadPhase2MultifitModules(process)

    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitPhase2ProducerPortable_cfi import (
        ecalUncalibRecHitPhase2ProducerPortable as _multifit,
    )
    from RecoLocalCalo.EcalRecProducers.ecalUncalibRecHitSoAToLegacy_cfi import (
        ecalUncalibRecHitSoAToLegacy as _soaToLegacy,
    )

    process.ecalMultiFitUncalibRecHitPhase2Portable = _multifit.clone(
        digisLabelEB=_DIGIS,
        recHitsLabelEB="EcalMultiFitUncalibRecHitsEB",
    )
    process.ecalMultiFitUncalibRecHitPhase2 = _soaToLegacy.clone(
        isPhase2=True,
        inputCollectionEB="ecalMultiFitUncalibRecHitPhase2Portable:EcalMultiFitUncalibRecHitsEB",
        inputCollectionEE=None,
        outputLabelEE=None,
    )

    process.ecalMultifitPhase2Task = cms.Task(
        process.ecalMultiFitUncalibRecHitPhase2Portable,
        process.ecalMultiFitUncalibRecHitPhase2,
    )
    if hasattr(process, "reconstruction_step"):
        process.reconstruction_step.associate(process.ecalMultifitPhase2Task)
    else:
        raise RuntimeError("customisePhase2MultifitGPUAlongside: no 'reconstruction_step' in the process")

    for out in process.outputModules_().values():
        if out.type_() == "PoolOutputModule":
            out.outputCommands.append("keep *_ecalMultiFitUncalibRecHitPhase2_*_*")
    return process
