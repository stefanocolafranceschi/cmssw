#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoPhase2Portable_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoPhase2Portable_h

#include <vector>

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsDevice.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "DeclsForKernels.h"
#include "EcalMultifitParameters.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit {

  using InputProduct = EcalDigiPhase2DeviceCollection;
  using OutputProduct = EcalUncalibratedRecHitDeviceCollection;

  void launchKernels(Queue& queue,
                     InputProduct const& digisDevEB,
                     OutputProduct& uncalibRecHitsDevEB,
                     EcalMultifitConditionsDevice const& conditionsDev,
                     EcalMultifitParametersPhase2 const* paramsDev,
                     ConfigurationParametersPhase2 const& configParams);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit

#endif  // RecoLocalCalo_EcalRecProducers_plugins_alpaka_EcalUncalibRecHitMultiFitAlgoPhase2Portable_h
