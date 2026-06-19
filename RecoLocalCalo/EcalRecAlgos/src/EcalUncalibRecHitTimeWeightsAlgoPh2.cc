#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitTimeWeightsAlgoPh2.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"

#include "TVectorD.h"

/// Compute time
double EcalUncalibRecHitTimeWeightsAlgoPh2::time(const EcalDataFrame_Ph2 &dataFrame,
                                                 const std::vector<double> &amplitudes,
                                                 const EcalLiteDTUPedestalsMap::Item *aped,
                                                 const EcalCATIAGainRatio *aGain,
                                                 const FullSampleVector &fullpulse,
                                                 const EcalPh2WeightSet::EcalWeightMatrix **weights) {
  const unsigned int nsample = EcalDataFrame_Ph2::MAXSAMPLES;

  double maxamplitude = -std::numeric_limits<double>::max();

  double pulsenorm = 0.;
  int iGainSwitch = 0;

  ROOT::Math::SVector<double, nsample> pedSubSamples;
  for (unsigned int iSample = 0; iSample < nsample; ++iSample) {
    const EcalLiteDTUSample &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    const int gainId = sample.gainId();

    double pedestal = aped->mean(gainId);
    double gainratio = *aGain;

    if (gainId == ecalPh2::gainId1) {
      iGainSwitch = 1;
    } else {
      iGainSwitch = 0;
    }

    amplitude = (static_cast<double>(sample.adc()) - pedestal) * gainratio;

    pedSubSamples(iSample) = amplitude;

    if (amplitude > maxamplitude) {
      maxamplitude = amplitude;
    }
    pulsenorm += fullpulse(iSample);
  }

  int ipulse = 0;
  for (const auto &amplit : amplitudes) {
    int bx = ipulse / ecalPh2::kNSamplesPerLHCPeriod - ecalPh2::kNActiveLHCPeriods / 2;
    int firstsamplet = std::max(0, bx + 3);
    int offset = ecalPh2::kNSamplesPerLHCPeriod * (ecalPh2::kNActiveLHCPeriods / 2);

    for (unsigned int isample = firstsamplet; isample < nsample; ++isample) {
      const auto pulse = fullpulse(isample + offset);
      pedSubSamples(isample) = std::max(0., pedSubSamples(isample) - amplit * pulse / pulsenorm);
    }
    ++ipulse;
  }

  // Compute parameters
  double amplitude(-1.), jitter(-1.);
  ROOT::Math::SVector<double, EcalUncalibRecHitRecAbsAlgo<EcalDataFrame_Ph2>::nWeightsRows> param =
      (*weights[iGainSwitch]) * pedSubSamples;
  amplitude = param(EcalUncalibRecHitRecAbsAlgo<EcalDataFrame_Ph2>::iAmplitude);
  if (amplitude) {
    jitter = -param(EcalUncalibRecHitRecAbsAlgo<EcalDataFrame_Ph2>::iTime) / amplitude;
  } else {
    jitter = 0.;
  }

  return jitter;
}
