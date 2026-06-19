#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitTimeWeightsAlgoPh2_h
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitTimeWeightsAlgoPh2_h

/** \class EcalUncalibRecHitTimeWeightsAlgoPh2
  *  Template used to compute amplitude, pedestal, time jitter, chi2 of a pulse
  *  using a weights method
  *
  *  The chi2 computation with matrix is replaced by the chi2express which is  moved outside the weight algo
  *  (need to clean up the interface in next iteration so that we do not pass-by useless arrays)
  */

#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSetT.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include <vector>

class EcalUncalibRecHitTimeWeightsAlgoPh2 {
public:
  using FullSampleVector = typename EigenMatrixTypes<ecalPh2>::FullSampleVector;

  EcalUncalibRecHitTimeWeightsAlgoPh2(){};
  ~EcalUncalibRecHitTimeWeightsAlgoPh2() = default;

  /// Compute time
  double time(const EcalDataFrame_Ph2 &dataFrame,
              const std::vector<double> &amplitudes,
              const EcalLiteDTUPedestalsMap::Item *aped,
              const EcalCATIAGainRatio *aGain,
              const FullSampleVector &fullpulse,
              const EcalPh2WeightSet::EcalWeightMatrix **weights);
};
#endif
