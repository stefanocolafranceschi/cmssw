#ifndef CondFormats_EcalObjects_EcalPh2SamplesCorrelation_h
#define CondFormats_EcalObjects_EcalPh2SamplesCorrelation_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include <iostream>
#include <vector>

class EcalPh2SamplesCorrelation {
public:
  EcalPh2SamplesCorrelation();
  EcalPh2SamplesCorrelation(const EcalPh2SamplesCorrelation& aset);
  ~EcalPh2SamplesCorrelation() = default;

  std::vector<double> g10SamplesCorrelation;
  std::vector<double> g1SamplesCorrelation;

  void print(std::ostream& o) const;

  COND_SERIALIZABLE;
};

#endif
