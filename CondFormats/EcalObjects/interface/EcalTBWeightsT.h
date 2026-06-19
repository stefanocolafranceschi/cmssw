#ifndef CondFormats_EcalObjects_EcalTBWeightsT_H
#define CondFormats_EcalObjects_EcalTBWeightsT_H

#include "CondFormats/Serialization/interface/Serializable.h"

#include <map>
#include "CondFormats/EcalObjects/interface/EcalXtalGroupId.h"
#include "CondFormats/EcalObjects/interface/EcalWeightSetT.h"

template <class P>
class EcalTBWeightsT {
public:
  using EcalTDCId = int;
  using EcalTBWeightMap = std::map<std::pair<EcalXtalGroupId, EcalTDCId>, EcalWeightSetT<P>>;

  EcalTBWeightsT();
  ~EcalTBWeightsT() = default;

  // modifiers
  void setValue(const EcalXtalGroupId& groupId, const EcalTDCId& tdcId, const EcalWeightSetT<P>& weight);
  void setValue(const std::pair<EcalXtalGroupId, EcalTDCId>& keyPair, const EcalWeightSetT<P>& weight);

  // accessors
  const EcalTBWeightMap& getMap() const { return map_; }

protected:
  EcalTBWeightMap map_;

  COND_SERIALIZABLE;
};

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
using EcalPh2TBWeights = EcalTBWeightsT<ecalPh2>;

#endif
