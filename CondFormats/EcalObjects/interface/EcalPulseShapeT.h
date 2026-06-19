#ifndef CondFormats_EcalObjects_EcalPulseShapeT_h
#define CondFormats_EcalObjects_EcalPulseShapeT_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/EcalObjects/interface/EcalCondObjectContainer.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

template <class P>
struct EcalPulseShapeT {
public:
  static constexpr int TEMPLATESAMPLES = static_cast<int>(P::kPulseShapeTemplateSampleSize);

  EcalPulseShapeT();

  float pdfval[TEMPLATESAMPLES];

  float val(int isample) const { return pdfval[isample]; }

  COND_SERIALIZABLE;
};

using EcalPh2PulseShape = EcalPulseShapeT<ecalPh2>;

typedef EcalCondObjectContainer<EcalPh2PulseShape> EcalPh2PulseShapesMap;
typedef EcalPh2PulseShapesMap::const_iterator EcalPh2PulseShapesMapIterator;
typedef EcalPh2PulseShapesMap EcalPh2PulseShapes;

#endif
