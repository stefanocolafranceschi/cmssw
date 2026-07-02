#ifndef CondFormats_EcalObjects_interface_alpaka_EcalMultifitConditionsPhase2Device_h
#define CondFormats_EcalObjects_interface_alpaka_EcalMultifitConditionsPhase2Device_h

#include "CondFormats/EcalObjects/interface/EcalMultifitConditionsPhase2Host.h"
#include "CondFormats/EcalObjects/interface/EcalMultifitConditionsPhase2SoA.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using ::EcalMultifitConditionsPhase2Host;
  using EcalMultifitConditionsPhase2Device = PortableCollection<EcalMultifitConditionsPhase2SoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
