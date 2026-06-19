#include "CondFormats/EcalObjects/interface/EcalWeightSetT.h"

template <class P>
EcalWeightSetT<P>::EcalWeightSetT() {}

template <class P>
EcalWeightSetT<P>::EcalWeightSetT(const EcalWeightSetT<P>& rhs) {
  wgtBeforeSwitch_ = rhs.wgtBeforeSwitch_;
  wgtAfterSwitch_ = rhs.wgtAfterSwitch_;
  wgtChi2BeforeSwitch_ = rhs.wgtChi2BeforeSwitch_;
  wgtChi2AfterSwitch_ = rhs.wgtChi2AfterSwitch_;
}

template <class P>
EcalWeightSetT<P>& EcalWeightSetT<P>::operator=(const EcalWeightSetT<P>& rhs) {
  wgtBeforeSwitch_ = rhs.wgtBeforeSwitch_;
  wgtAfterSwitch_ = rhs.wgtAfterSwitch_;
  wgtChi2BeforeSwitch_ = rhs.wgtChi2BeforeSwitch_;
  wgtChi2AfterSwitch_ = rhs.wgtChi2AfterSwitch_;
  return *this;
}

template class EcalWeightSetT<ecalPh2>;
