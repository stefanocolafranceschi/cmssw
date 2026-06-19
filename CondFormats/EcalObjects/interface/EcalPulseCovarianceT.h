#ifndef CondFormats_EcalObjects_EcalPulseCovarianceT_h
#define CondFormats_EcalObjects_EcalPulseCovarianceT_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"

template <typename T>
struct EcalPulseCovarianceT {
public:
  static const int TEMPLATESAMPLES = T::TEMPLATESAMPLES;

  EcalPulseCovarianceT();

  float covval[TEMPLATESAMPLES][TEMPLATESAMPLES];

  float val(int i, int j) const { return covval[i][j]; }

  COND_SERIALIZABLE;
};

using EcalPh2PulseCovariance = EcalPulseCovarianceT<EcalPh2PulseShape>;

typedef EcalCondObjectContainer<EcalPh2PulseCovariance> EcalPh2PulseCovariancesMap;
typedef EcalPh2PulseCovariancesMap::const_iterator EcalPh2PulseCovariancesMapIterator;
typedef EcalPh2PulseCovariancesMap EcalPh2PulseCovariances;

#endif
