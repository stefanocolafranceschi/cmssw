#ifndef CondFormats_EcalObjects_EcalMultifitConditionsPhase2SoA_h
#define CondFormats_EcalObjects_EcalMultifitConditionsPhase2SoA_h

#include <array>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"

// Multifit conditions for the Phase-2 ECAL barrel (CATIA amplifier with
// 2 gains x10/x1, LiteDTU ADC with 16 samples per frame). Phase-2 ECAL is
// barrel-only, so unlike EcalMultifitConditionsSoA there is no EE section
// and no offsetEE. Timing-related conditions are not included yet; they will
// be added when the Phase-2 GPU timing computation is addressed.
using PulseShapePhase2Array = std::array<float, EcalPh2PulseShape::TEMPLATESAMPLES>;
using SampleCorrelationPhase2Array = std::array<double, ecalPh2::sampleSize>;

// Pulse-shape covariance stored ROW-MAJOR and CONTIGUOUS per channel
// (flat [tr*TEMPLATESAMPLES + tc], layout-compatible with
// EcalPh2PulseCovariance::covval), deliberately as a plain SOA_COLUMN and
// NOT as a SOA_EIGEN_COLUMN: SoA Eigen columns store each matrix component
// in its own stride-separated column, which is incompatible with the
// kernels' per-channel EcalPh2PulseCovariance aliasing. (Run 3 has exactly
// this Eigen-column fill / contiguous reinterpret_cast read mismatch in
// EcalMultifitConditionsSoA + AmplitudeComputationKernels.dev.cc -- flagged
// upstream; do not copy that pairing here.)
using PulseCovariancePhase2Array =
    std::array<float, EcalPh2PulseShape::TEMPLATESAMPLES * EcalPh2PulseShape::TEMPLATESAMPLES>;

GENERATE_SOA_LAYOUT(EcalMultifitConditionsPhase2SoALayout,
                    SOA_COLUMN(uint32_t, rawid),
                    // pedestal mean and rms per CATIA gain (indices: ecalPh2::gainId10, ecalPh2::gainId1)
                    SOA_COLUMN(float, pedestals_mean_g10),
                    SOA_COLUMN(float, pedestals_mean_g1),
                    SOA_COLUMN(float, pedestals_rms_g10),
                    SOA_COLUMN(float, pedestals_rms_g1),
                    // CATIA gain ratio g10/g1 (EcalCATIAGainRatios), nominally 10
                    SOA_COLUMN(float, gain10Over1),
                    // 16-sample pulse shape template (2017 test-beam simulation)
                    SOA_COLUMN(PulseShapePhase2Array, pulseShapes),
                    // NxN N=TEMPLATESAMPLES(16) for each xtal, row-major flat
                    SOA_COLUMN(PulseCovariancePhase2Array, pulseCovariance),
                    // Sample correlation scalars: array of 16 values per CATIA gain (barrel only)
                    SOA_SCALAR(SampleCorrelationPhase2Array, sampleCorrelation_g10),
                    SOA_SCALAR(SampleCorrelationPhase2Array, sampleCorrelation_g1))

using EcalMultifitConditionsPhase2SoA = EcalMultifitConditionsPhase2SoALayout<>;

#endif
