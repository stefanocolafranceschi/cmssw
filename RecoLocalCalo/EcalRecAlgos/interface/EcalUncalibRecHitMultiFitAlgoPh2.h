#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoPh2_h
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMultiFitAlgoPh2_h

/** \class EcalUncalibRecHitMultiFitAlgoPh2
  *  Amplitude reconstucted from Phase 2 digis by the multi-template fit
  */

#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"

class EcalUncalibRecHitMultiFitAlgoPh2 {
public:
  using SampleVector = typename EigenMatrixTypes<ecalPh2>::SampleVector;
  using FullSampleVector = typename EigenMatrixTypes<ecalPh2>::FullSampleVector;
  using BXVector = typename EigenMatrixTypes<ecalPh2>::BXVector;
  using SampleGainVector = typename EigenMatrixTypes<ecalPh2>::SampleGainVector;
  using SampleMatrix = typename EigenMatrixTypes<ecalPh2>::SampleMatrix;
  using FullSampleMatrix = typename EigenMatrixTypes<ecalPh2>::FullSampleMatrix;
  using SampleMatrixGainArray = typename EigenMatrixTypes<ecalPh2>::SampleMatrixGainArray;

  EcalUncalibRecHitMultiFitAlgoPh2();
  ~EcalUncalibRecHitMultiFitAlgoPh2(){};
  EcalUncalibratedRecHit makeRecHit(const EcalDataFrame_Ph2 &dataFrame,
                                    const EcalLiteDTUPedestalsMap::Item *aped,
                                    const EcalCATIAGainRatio *aGain,
                                    const SampleMatrixGainArray &noisecors,
                                    const FullSampleVector &fullpulse,
                                    const FullSampleMatrix &fullpulsecov,
                                    const BXVector &activeBX);
  void disableErrorCalculation() { computeErrors_ = false; }
  void setDoPrefit(const bool b) { doPrefit_ = b; }
  void setPrefitMaxChiSq(const double x) { prefitMaxChiSq_ = x; }
  void setDynamicPedestals(const bool b) { dynamicPedestals_ = b; }
  void setMitigateBadSamples(const bool b) { mitigateBadSamples_ = b; }
  void setSelectiveBadSampleCriteria(const bool b) { selectiveBadSampleCriteria_ = b; }
  void setAddPedestalUncertainty(const double x) { addPedestalUncertainty_ = x; }
  void setSimplifiedNoiseModelForGainSwitch(const bool b) { simplifiedNoiseModelForGainSwitch_ = b; }
  void setGainSwitchUseMaxSample(const bool b) { gainSwitchUseMaxSample_ = b; }

private:
  PulseChiSqSNNLS<ecalPh2> pulsefunc_;
  PulseChiSqSNNLS<ecalPh2> pulsefuncSingle_;
  bool computeErrors_;
  bool doPrefit_;
  double prefitMaxChiSq_;
  bool dynamicPedestals_;
  bool mitigateBadSamples_;
  bool selectiveBadSampleCriteria_;
  double addPedestalUncertainty_;
  bool simplifiedNoiseModelForGainSwitch_;
  bool gainSwitchUseMaxSample_;
  BXVector singlebx_;
};

#endif
