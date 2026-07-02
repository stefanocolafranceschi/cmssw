#ifndef CondFormats_EcalObjects_interface_EcalMultifitConditionsPhase2Host_h
#define CondFormats_EcalObjects_interface_EcalMultifitConditionsPhase2Host_h

#include "CondFormats/EcalObjects/interface/EcalMultifitConditionsPhase2SoA.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

using EcalMultifitConditionsPhase2Host = PortableHostCollection<EcalMultifitConditionsPhase2SoA>;

#endif
