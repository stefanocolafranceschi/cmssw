#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"

template <class P>
EcalPulseShapeT<P>::EcalPulseShapeT() {
  for (int s = 0; s < TEMPLATESAMPLES; ++s)
    pdfval[s] = 0.;
}

template struct EcalPulseShapeT<ecalPh2>;
