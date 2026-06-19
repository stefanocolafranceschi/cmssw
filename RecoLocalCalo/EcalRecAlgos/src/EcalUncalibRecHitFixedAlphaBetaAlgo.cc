/** \class EcalUncalibRecHitFixedAlphaBetaAlgo
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using an analytical function fit, with the pulse parameters alpha and beta fixed.
  *  It follows a fast fit algorithms devolped on test beam data by P. Jarry
  *  If the pedestal is not given, it is calculated from the first 2 pre-samples; 
  *
  *  \author A.Ghezzi
  */

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitFixedAlphaBetaAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

template <class C, class P>
EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::EcalUncalibRecHitFixedAlphaBetaAlgo()
    : minAmpl_(16.),
      dyn_pedestal_(true),
      fAlpha_(0.),
      fBeta_(0.),
      fAmp_max_(-1.),
      fTim_max_(-1.),
      fPed_max_(0.),
      alfabeta_(0.),
      fNb_iter_(4),
      fNum_samp_bef_max_(1),
      fNum_samp_after_max_(3),
      fSigma_ped_(1.1),
      un_sur_sigma_(1. / double(fSigma_ped_)),
      doFit_(false),
      DM1_(3),
      temp_(3) {
  for (int i = 0; i < 36; ++i) {
    for (int j = 0; j < 1701; ++j) {
      alpha_table_[i][j] = 1.2;
      beta_table_[i][j] = 1.7;
    }
  }
}

template <class C, class P>
EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::EcalUncalibRecHitFixedAlphaBetaAlgo(int n_iter,
                                                                               int n_bef_max,
                                                                               int n_aft_max,
                                                                               float sigma_ped)
    : EcalUncalibRecHitFixedAlphaBetaAlgo() {
  fNb_iter_ = n_iter;
  fNum_samp_bef_max_ = n_bef_max;
  fNum_samp_after_max_ = n_aft_max;
  fSigma_ped_ = sigma_ped;
  un_sur_sigma_ = 1. / double(sigma_ped);
}

/// Compute parameters
template <class C, class P>
EcalUncalibratedRecHit EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::makeRecHit(
    const C& dataFrame,
    const double* pedestals,
    const double* gainRatios,
    const EcalWeightSet::EcalWeightMatrix** weights,
    const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix) {
  double chi2_(-1.);

  double frame[C::MAXSAMPLES];  // will contain the ADC values
  double pedestal = 0;          // carries pedestal for highest gain i.e. gainId == 1 (Phase 1) gainId == 0 (Phase 2)

  constexpr int gainId0 =
      (P::NGAINS > 2) ? 1 : 0;  // expected gainId at the beginning of dataFrame. 1 for Phase 1 and 0 for Phase 2
  bool iGainSwitch = false;     // flags whether there's any gainId other than gainId0
  int gainId = 0;               // stores gainId at every sample
  double maxsample(-1);         // ADC value of maximal ped-subtracted sample
  int imax(-1);                 // sample number of maximal ped-subtracted sample
  bool external_pede = false;
  bool isSaturated = false;  // flag reporting whether gain0 has been found

  // Set the pedestal
  if (pedestals) {
    external_pede = true;
    if (dyn_pedestal_) {
      pedestal = (double(dataFrame.sample(0).adc() + dataFrame.sample(1).adc())) / 2.;
    } else {
      pedestal = pedestals[0];
    }
  } else {  // pedestal from pre-sample
    external_pede = false;
    pedestal = double(dataFrame.sample(0).adc() + dataFrame.sample(1).adc()) / 2.;
  }

  // Get time samples checking for Gain Switch and pedestals
  for (int iSample = 0; iSample < C::MAXSAMPLES; ++iSample) {
    //create frame in adc max gain equivalent
    gainId = dataFrame.sample(iSample).gainId();

    // warning: For Phase 1 the vector pedestal is supposed to have in the order G12, G6 and G1
    // if gainId is zero treat it as 3 temporarily to protect against undefined
    // frame will be set to ~max of gain1
    if (P::NGAINS > 2 && gainId == 0) {
      gainId = 3;
      isSaturated = true;
    }

    if (pedestals && gainId != gainId0) {
      frame[iSample] = (double(dataFrame.sample(iSample).adc()) - pedestals[gainId - 1]) * gainRatios[gainId - 1];
    } else {
      frame[iSample] = double(dataFrame.sample(iSample).adc()) - pedestal;
    }

    iGainSwitch = gainId != gainId0;

    if (frame[iSample] > maxsample) {
      maxsample = frame[iSample];
      imax = iSample;
    }
  }

  if ((iGainSwitch && !external_pede) ||  // ... thus you return dummy rechit
      imax == -1) {                       // protect against all frames being <-1
    return EcalUncalibratedRecHit(dataFrame.id(), -1., -100., -1., -1.);
  }

  InitFitParameters(frame, imax);
  chi2_ = PerformAnalyticFit(frame, imax);
  uint32_t flags = 0;
  if (isSaturated)
    flags = EcalUncalibratedRecHit::kSaturated;

  return EcalUncalibratedRecHit(dataFrame.id(), fAmp_max_, pedestal + fPed_max_, fTim_max_ - 5, chi2_, flags);
}

template <class C, class P>
double EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::pulseShapeFunction(const double t) const {
  if (alfabeta_ <= 0.)
    return 0.;
  double dtsbeta, variable, puiss;
  const double dt = t - fTim_max_;
  if (dt > -alfabeta_) {
    dtsbeta = dt / fBeta_;
    variable = 1. + dt / alfabeta_;
    puiss = pow(variable, fAlpha_);
    return fAmp_max_ * puiss * exp(-dtsbeta) + fPed_max_;
  }
  return fPed_max_;
}

template <class C, class P>
void EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::InitFitParameters(const double* samples, const int max_sample) {
  // in a first attempt just use the value of the maximum sample
  fAmp_max_ = samples[max_sample];
  fTim_max_ = max_sample;
  fPed_max_ = 0;

  // amplitude too low for fit to converge
  // timing set correctly is assumed
  if (fAmp_max_ < minAmpl_) {
    fAmp_max_ = samples[5];
    double sumA = samples[5] + samples[4] + samples[6];
    if (sumA != 0) {
      fTim_max_ = 5 + (samples[6] - samples[4]) / sumA;
    } else {
      fTim_max_ = -993;
    }  //-999+6
    doFit_ = false;
  }
  // if timing very badly off, that just use max sample
  else if (max_sample < 1 || max_sample > 7) {
    doFit_ = false;
  } else {
    //y=a*(x-xM)^2+b*(x-xM)+c
    doFit_ = true;
    float a = float(samples[max_sample - 1] + samples[max_sample + 1] - 2 * samples[max_sample]) / 2.;
    if (a == 0) {
      doFit_ = false;
      return;
    }
    float b = float(samples[max_sample + 1] - samples[max_sample - 1]) / 2.;

    fTim_max_ = max_sample - b / (2 * a);
    fAmp_max_ = samples[max_sample] - b * b / (4 * a);
  }
}

template <class C, class P>
float EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::PerformAnalyticFit(const double* samples, const int max_sample) {
  //! fit electronic function from simulation
  //! parameters fAlpha_ and fBeta_ are fixed and fit is providing the 3 following parameters
  //! the maximum amplitude ( fAmp_max_ )
  //! the time of the maximum  ( fTim_max_)
  //| the pedestal (fPed_max_)

  double chi2 = -1, db[3];

  int num_fit_min = (int)(max_sample - fNum_samp_bef_max_);
  int num_fit_max = (int)(max_sample + fNum_samp_after_max_);

  if (num_fit_min < 0)
    num_fit_min = 0;
  if (num_fit_max >= C::MAXSAMPLES) {
    num_fit_max = C::MAXSAMPLES - 1;
  }

  if (!doFit_) {
    LogDebug("EcalUncalibRecHitFixedAlphaBetaAlgo") << "No fit performed. The amplitude of sample 5 will be used";
    return -1;
  }

  double func, delta;
  double variation_func_max = 0.;
  double variation_tim_max = 0.;
  double variation_ped_max = 0.;
  //!          Loop on iterations
  for (int iter = 0; iter < fNb_iter_; ++iter) {
    //!          initialization inside iteration loop !
    chi2 = 0.;

    for (int i1 = 0; i1 < 3; ++i1) {
      temp_[i1] = 0;
      for (int j1 = i1; j1 < 3; ++j1) {
        DM1_.fast(j1 + 1, i1 + 1) = 0;
      }
    }

    fAmp_max_ += variation_func_max;
    fTim_max_ += variation_tim_max;
    fPed_max_ += variation_ped_max;

    //! Then we loop on samples to be fitted
    for (int i = num_fit_min; i <= num_fit_max; ++i) {
      //! calculate function to be fitted
      func = pulseShapeFunction((double)i);
      //! then calculate derivatives of function to be fitted
      const double dt = (double)i - fTim_max_;
      if (dt > -alfabeta_) {
        const double dt_sur_beta = dt / fBeta_;
        const double variable = 1. + dt / alfabeta_;
        const double expo = exp(-dt_sur_beta);
        const double puissance = pow(variable, fAlpha_);

        db[0] = un_sur_sigma_ * puissance * expo;
        db[1] = fAmp_max_ * db[0] * dt_sur_beta / (alfabeta_ * variable);
      } else {
        db[0] = 0.;
        db[1] = 0.;
      }
      db[2] = un_sur_sigma_;
      //! compute matrix elements DM1
      for (int i1 = 0; i1 < 3; ++i1) {
        for (int j1 = i1; j1 < 3; ++j1) {
          //double & fast(int row, int col);
          DM1_.fast(j1 + 1, i1 + 1) += db[i1] * db[j1];
        }
      }
      //! compute delta
      delta = (samples[i] - func) * un_sur_sigma_;
      //! compute vector elements PROD
      for (int ii = 0; ii < 3; ++ii) {
        temp_[ii] += delta * db[ii];
      }
      chi2 += delta * delta;
    }  //! end of loop on samples

    int fail = 0;
    DM1_.invert(fail);
    if (fail != 0.) {
      //just a guess from the value of the parameters in the previous interaction;
      InitFitParameters(samples, max_sample);
      return -101;
    }
    //! compute variations of parameters fAmp_max and fTim_max
    CLHEP::HepVector PROD = DM1_ * temp_;

    // Probably the fastest way to protect against
    // +-inf value in the matrix DM1_ after inversion
    // (which is nevertheless flagged as successfull...)
    if (edm::isNotFinite(PROD[0])) {
      InitFitParameters(samples, max_sample);
      return -103;
    }

    variation_func_max = PROD[0];
    variation_tim_max = PROD[1];
    variation_ped_max = PROD[2];
  }  //!end of loop on iterations

  //!   protection again diverging/unstable fit
  if (variation_func_max > 2000. || variation_func_max < -1000.) {
    InitFitParameters(samples, max_sample);
    return -102;
  }

  //!      results of the fit are calculated
  fAmp_max_ += variation_func_max;
  fTim_max_ += variation_tim_max;
  fPed_max_ += variation_ped_max;

  // protection against unphysical results:
  // ampli mismatched to MaxSample, ampli largely negative, time off preselected range
  if (fAmp_max_ > 2 * samples[max_sample] || fAmp_max_ < -100 || fTim_max_ < 0 || 9 < fTim_max_) {
    InitFitParameters(samples, max_sample);
    return -104;
  }

  return chi2;
}

template <class C, class P>
void EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::SetAlphaBeta(const double alpha, const double beta) {
  fAlpha_ = alpha;
  fBeta_ = beta;
  alfabeta_ = fAlpha_ * fBeta_;
}

template <class C, class P>
void EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::SetMinAmpl(const double ampl) {
  minAmpl_ = ampl;
}
template <class C, class P>
void EcalUncalibRecHitFixedAlphaBetaAlgo<C, P>::SetDynamicPedestal(const bool p) {
  dyn_pedestal_ = p;
}

// Define type combinations that the algorithm can be used with
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
template class EcalUncalibRecHitFixedAlphaBetaAlgo<EBDataFrame, ecalPh1>;
template class EcalUncalibRecHitFixedAlphaBetaAlgo<EEDataFrame, ecalPh1>;
template class EcalUncalibRecHitFixedAlphaBetaAlgo<EcalDataFrame_Ph2, ecalPh2>;
