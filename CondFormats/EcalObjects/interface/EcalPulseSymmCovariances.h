#ifndef EcalPulseSymmCovariances_h
#define EcalPulseSymmCovariances_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"

#include <algorithm>

struct EcalPulseSymmCovariance {
public:
  static const int TEMPLATESAMPLES = EcalPulseShape::TEMPLATESAMPLES;

  EcalPulseSymmCovariance();

  float covval[TEMPLATESAMPLES * (TEMPLATESAMPLES + 1) / 2];

  int indexFor(int i, int j) const {
    int m = std::min(i, j);
    int n = std::max(i, j);
    return n + TEMPLATESAMPLES * m - m * (m + 1) / 2;
  }

  float val(int i, int j) const { return covval[indexFor(i, j)]; }
  float& val(int i, int j) { return covval[indexFor(i, j)]; }

  COND_SERIALIZABLE;
};

typedef EcalCondObjectContainer<EcalPulseSymmCovariance> EcalPulseSymmCovariancesMap;
typedef EcalPulseSymmCovariancesMap::const_iterator EcalPulseSymmCovariancesMapIterator;
typedef EcalPulseSymmCovariancesMap EcalPulseSymmCovariances;

struct EcalPh2PulseSymmCovariance : public EcalPulseSymmCovariance {
public:
  static const int TEMPLATESAMPLES = EcalPh2PulseShape::TEMPLATESAMPLES;

  COND_SERIALIZABLE;
};

typedef EcalCondObjectContainer<EcalPulseSymmCovariance> EcalPh2PulseSymmCovariancesMap;
typedef EcalPh2PulseSymmCovariancesMap::const_iterator EcalPh2PulseSymmCovariancesMapIterator;
typedef EcalPh2PulseSymmCovariancesMap EcalPh2PulseSymmCovariances;

#endif
