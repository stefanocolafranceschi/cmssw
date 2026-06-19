#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernelsPhase2_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationKernelsPhase2_h

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsDevice.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "DeclsForKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit {

  using InputProduct = EcalDigiPhase2DeviceCollection;
  using OutputProduct = EcalUncalibratedRecHitDeviceCollection;

  void minimization_procedure(Queue& queue,
                              InputProduct const& digisDevEB,
                              OutputProduct& uncalibRecHitsDevEB,
                              EventDataForScratchDevicePhase2& scratch,
                              EcalMultifitConditionsDevice const& conditionsDev,
                              ConfigurationParametersPhase2 const& configParams,
                              uint32_t const totalChannels);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationKernelsPhase2_h
