#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitFixedAlphaBetaAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitFixedAlphaBetaAlgo_HH

/** \class EcalUncalibRecHitFixedAlphaBetaAlgo
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using an analytical function fit, with the pulse parameters alpha and beta fixed.
  *  It follows a fast fit algorithms devolped on test beam data by P. Jarry
  *  If the pedestal is not given, it is calculated from the first 2 pre-samples; 
  *
  *  \author A.Ghezzi
  */

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"

template <class C, class P>
class EcalUncalibRecHitFixedAlphaBetaAlgo : public EcalUncalibRecHitRecAbsAlgo<C> {
public:
  EcalUncalibRecHitFixedAlphaBetaAlgo();
  EcalUncalibRecHitFixedAlphaBetaAlgo(int n_iter, int n_bef_max = 1, int n_aft_max = 3, float sigma_ped = 1.1);
  ~EcalUncalibRecHitFixedAlphaBetaAlgo() override = default;
  EcalUncalibratedRecHit makeRecHit(const C& dataFrame,
                                    const double* pedestals,
                                    const double* gainRatios,
                                    const EcalWeightSet::EcalWeightMatrix** weights,
                                    const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix) override;
  void SetAlphaBeta(const double alpha, const double beta);
  void SetMinAmpl(const double ampl);
  void SetDynamicPedestal(const bool dyn_pede);

private:
  double minAmpl_;
  bool dyn_pedestal_;
  double fAlpha_;    //parameter of the shape
  double fBeta_;     //parameter of the shape
  double fAmp_max_;  // peak amplitude
  double fTim_max_;  // time of the peak (in 25ns units)
  double fPed_max_;  // pedestal value
  double alfabeta_;

  int fNb_iter_;
  int fNum_samp_bef_max_;
  int fNum_samp_after_max_;

  float fSigma_ped_;
  double un_sur_sigma_;
  //temporary solution for different alpha and beta
  float alpha_table_[36][1701];
  float beta_table_[36][1701];

  bool doFit_;

  CLHEP::HepSymMatrix DM1_;
  CLHEP::HepVector temp_;

  double pulseShapeFunction(const double t) const;
  float PerformAnalyticFit(const double* samples, const int max_sample);
  void InitFitParameters(const double* samples, const int max_sample);
};

#endif
