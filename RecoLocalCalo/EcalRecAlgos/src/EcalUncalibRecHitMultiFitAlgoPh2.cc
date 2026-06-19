#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgoPh2.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

EcalUncalibRecHitMultiFitAlgoPh2::EcalUncalibRecHitMultiFitAlgoPh2()
    : computeErrors_(true),
      doPrefit_(false),
      prefitMaxChiSq_(1.),
      dynamicPedestals_(false),
      mitigateBadSamples_(false),
      selectiveBadSampleCriteria_(false),
      addPedestalUncertainty_(0.),
      simplifiedNoiseModelForGainSwitch_(true),
      gainSwitchUseMaxSample_(false) {
  singlebx_.resize(1);
  singlebx_ << 0;

  pulsefuncSingle_.disableErrorCalculation();
  pulsefuncSingle_.setMaxIters(1);
  pulsefuncSingle_.setMaxIterWarnings(false);
}

/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgoPh2::makeRecHit(const EcalDataFrame_Ph2 &dataFrame,
                                                                    const EcalLiteDTUPedestalsMap::Item *aped,
                                                                    const EcalCATIAGainRatio *aGain,
                                                                    const SampleMatrixGainArray &noisecors,
                                                                    const FullSampleVector &fullpulse,
                                                                    const FullSampleMatrix &fullpulsecov,
                                                                    const BXVector &activeBX) {
  const uint32_t flags = 0;

  constexpr unsigned int nsample = EcalDataFrame_Ph2::MAXSAMPLES;

  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;

  double pedval = 0.;

  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  const bool hasGainSwitch = false;

  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = dynamicPedestals_ && !hasGainSwitch;

  for (unsigned int iSample = 0; iSample < nsample; ++iSample) {
    const auto &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    const int gainId = sample.gainId();

    const double pedestal = aped->mean(gainId);
    const double gainratio = *aGain;

    if (gainId == 1) {
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    } else {
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1;  //-1 for static pedestal
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc()) * gainratio;
    } else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }

    amplitudes[iSample] = amplitude;

    if (iSample == iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }
  }

  double amplitude, amperr, chisq;
  bool status = false;

  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && gainSwitchUseMaxSample_) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh(dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags);
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse = 0; ipulse < pulsefunc_.BXs().rows(); ++ipulse) {
      int bx = pulsefunc_.BXs().coeff(ipulse);
      if (bx != 0) {
        rh.setOutOfTimeAmplitude(bx + 5, 0.0);
      }
    }
    return rh;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = mitigateBadSamples_ && hasGainSwitch && iSampleMax > 0;
  mitigateBadSample &=
      (!selectiveBadSampleCriteria_ || (gainsNoise.coeff(iSampleMax - 1) != gainsNoise.coeff(iSampleMax)));
  if (mitigateBadSample) {
    badSamples[iSampleMax - 1] = 1;
  }

  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    std::array<double, ecalPh2::NGAINS> pedrmss;
    std::array<double, ecalPh2::NGAINS> gainratios;
    for (unsigned int i; i < ecalPh2::NGAINS; ++i) {
      pedrmss[i] = aped->rms(i);
      gainratios[i] = ecalPh2::gains[i];
    }
    if (simplifiedNoiseModelForGainSwitch_) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax] * gainratios[gainidxmax] * pedrmss[gainidxmax] * pedrmss[gainidxmax] *
                 noisecors[gainidxmax];
      if (!dynamicPedestal && addPedestalUncertainty_ > 0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += addPedestalUncertainty_ * addPedestalUncertainty_ * SampleMatrix::Ones();
      }
    } else {
      noisecov = SampleMatrix::Zero();
      for (unsigned int gainidx = 0; gainidx < noisecors.size(); ++gainidx) {
        SampleGainVector mask = gainidx * SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array() == mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff() > 0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx] * gainratios[gainidx] * pedrmss[gainidx] * pedrmss[gainidx] *
                      pedestal.asDiagonal() * noisecors[gainidx] * pedestal.asDiagonal();
          if (!dynamicPedestal && addPedestalUncertainty_ > 0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx] * gainratios[gainidx] * addPedestalUncertainty_ * addPedestalUncertainty_ *
                        pedestal.asDiagonal() * SampleMatrix::Ones() * pedestal.asDiagonal();
          }
        }
      }
    }
  } else {
    noisecov = aped->rms(ecalPh2::gainId10) * aped->rms(ecalPh2::gainId10) * noisecors[0];
    if (!dynamicPedestal && addPedestalUncertainty_ > 0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += addPedestalUncertainty_ * addPedestalUncertainty_ * SampleMatrix::Ones();
    }
  }

  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (doPrefit_) {
    status =
        pulsefuncSingle_.DoFit(amplitudes, noisecov, singlebx_, fullpulse, fullpulsecov, gainsPedestal, badSamples);
    amplitude = status ? pulsefuncSingle_.X()[0] : 0.;
    amperr = status ? pulsefuncSingle_.Errors()[0] : 0.;
    chisq = pulsefuncSingle_.ChiSq();

    if (chisq < prefitMaxChiSq_) {
      usePrefit = true;
    }
  }

  if (!usePrefit) {
    if (!computeErrors_)
      pulsefunc_.disableErrorCalculation();
    status = pulsefunc_.DoFit(amplitudes, noisecov, activeBX, fullpulse, fullpulsecov, gainsPedestal, badSamples);
    chisq = pulsefunc_.ChiSq();

    if (!status) {
      edm::LogWarning("EcalUncalibRecHitMultiFitAlgoPh2::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse = 0; ipulse < pulsefunc_.BXs().rows(); ++ipulse) {
      if (pulsefunc_.BXs().coeff(ipulse) == 0) {
        ipulseintime = ipulse;
        break;
      }
    }

    amplitude = status ? pulsefunc_.X()[ipulseintime] : 0.;
    amperr = status ? pulsefunc_.Errors()[ipulseintime] : 0.;
  }

  double jitter = 0.;

  EcalUncalibratedRecHit rh(dataFrame.id(), amplitude, pedval, jitter, chisq, flags);
  rh.setAmplitudeError(amperr);

  if (!usePrefit) {
    for (unsigned int ipulse = 0; ipulse < pulsefunc_.BXs().rows(); ++ipulse) {
      int bx = pulsefunc_.BXs().coeff(ipulse);
      if (bx != 0 && std::abs(bx) < 100) {
        rh.setOutOfTimeAmplitude(bx + 5, status ? pulsefunc_.X().coeff(ipulse) : 0.);
      } else if (bx == (100 + gainsPedestal[iSampleMax])) {
        rh.setPedestal(status ? pulsefunc_.X().coeff(ipulse) : 0.);
      }
    }
  }

  return rh;
}
