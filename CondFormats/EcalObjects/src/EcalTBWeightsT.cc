#include "CondFormats/EcalObjects/interface/EcalTBWeightsT.h"

template <class P>
EcalTBWeightsT<P>::EcalTBWeightsT() {}

template <class P>
void EcalTBWeightsT<P>::setValue(const EcalXtalGroupId& groupId,
                                 const EcalTDCId& tdcId,
                                 const EcalWeightSetT<P>& weight) {
  setValue(std::make_pair(groupId, tdcId), weight);
}

template <class P>
void EcalTBWeightsT<P>::setValue(const std::pair<EcalXtalGroupId, EcalTDCId>& keyPair,
                                 const EcalWeightSetT<P>& weight) {
  map_.insert(std::make_pair(keyPair, weight));
}

template class EcalTBWeightsT<ecalPh2>;
