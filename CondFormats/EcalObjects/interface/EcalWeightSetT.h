#ifndef CondFormats_EcalObjects_EcalWeightSetT_h
#define CondFormats_EcalObjects_EcalWeightSetT_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "DataFormats/Math/interface/Matrix.h"
#include <iostream>

template <class P>
class EcalWeightSetT {
public:
  using EcalWeightMatrix =
      typename math::Matrix<3, P::sampleSize>::type;  // 3 different types of weights (amplitude, pedestal, time)
  using EcalChi2WeightMatrix = typename math::Matrix<P::sampleSize, P::sampleSize>::type;

  EcalWeightSetT();
  EcalWeightSetT(const EcalWeightSetT<P>& aset);
  ~EcalWeightSetT() = default;

  EcalWeightMatrix& getWeightsBeforeGainSwitch() { return wgtBeforeSwitch_; }
  EcalWeightMatrix& getWeightsAfterGainSwitch() { return wgtAfterSwitch_; }
  EcalChi2WeightMatrix& getChi2WeightsBeforeGainSwitch() { return wgtChi2BeforeSwitch_; }
  EcalChi2WeightMatrix& getChi2WeightsAfterGainSwitch() { return wgtChi2AfterSwitch_; }

  const EcalWeightMatrix& getWeightsBeforeGainSwitch() const { return wgtBeforeSwitch_; }
  const EcalWeightMatrix& getWeightsAfterGainSwitch() const { return wgtAfterSwitch_; }
  const EcalChi2WeightMatrix& getChi2WeightsBeforeGainSwitch() const { return wgtChi2BeforeSwitch_; }
  const EcalChi2WeightMatrix& getChi2WeightsAfterGainSwitch() const { return wgtChi2AfterSwitch_; }

  EcalWeightSetT& operator=(const EcalWeightSetT<P>& rhs);

  void print(std::ostream& o) const {
    using namespace std;
    o << "wgtBeforeSwitch_.: " << wgtBeforeSwitch_ << " wgtAfterSwitch_.: " << wgtAfterSwitch_
      << " wgtChi2BeforeSwitch_.: " << wgtChi2BeforeSwitch_ << " wgtChi2AfterSwitch_.: " << wgtChi2AfterSwitch_ << endl;
  }

private:
  EcalWeightMatrix wgtBeforeSwitch_;
  EcalWeightMatrix wgtAfterSwitch_;
  EcalChi2WeightMatrix wgtChi2BeforeSwitch_;
  EcalChi2WeightMatrix wgtChi2AfterSwitch_;

  COND_SERIALIZABLE;
};

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
using EcalPh2WeightSet = EcalWeightSetT<ecalPh2>;

#endif
