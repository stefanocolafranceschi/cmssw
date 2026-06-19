#include "CondFormats/EcalObjects/interface/EcalPh2SamplesCorrelation.h"

EcalPh2SamplesCorrelation::EcalPh2SamplesCorrelation() {}

EcalPh2SamplesCorrelation::EcalPh2SamplesCorrelation(const EcalPh2SamplesCorrelation& aset)
    : g10SamplesCorrelation(aset.g10SamplesCorrelation), g1SamplesCorrelation(aset.g1SamplesCorrelation) {}

template <typename T>
static inline void print_vector(std::ostream& o, const std::vector<T>& vect) {
  o << "[";
  for (const auto& corr : vect) {
    o << corr << ", ";
  }
  o << "]";
}

void EcalPh2SamplesCorrelation::print(std::ostream& o) const {
  o << "Gain 10 correlations:";
  print_vector<double>(o, g10SamplesCorrelation);
  o << std::endl;
  o << "Gain 1 correlations:";
  print_vector<double>(o, g1SamplesCorrelation);
  o << std::endl;
}
