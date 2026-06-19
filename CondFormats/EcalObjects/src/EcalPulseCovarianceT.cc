#include "CondFormats/EcalObjects/interface/EcalPulseCovarianceT.h"

template <typename T>
EcalPulseCovarianceT<T>::EcalPulseCovarianceT() {
  for (int k = 0; k < std::pow(TEMPLATESAMPLES, 2); ++k) {
    int i = k / TEMPLATESAMPLES;
    int j = k % TEMPLATESAMPLES;
    covval[i][j] = 0.;
  }
}

template struct EcalPulseCovarianceT<EcalPh2PulseShape>;
