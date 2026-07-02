#ifndef CondFormats_DataRecord_EcalMultifitConditionsPhase2Rcd_h
#define CondFormats_DataRecord_EcalMultifitConditionsPhase2Rcd_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"

#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLiteDTUPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2SamplesCorrelationRcd.h"

// Record for the Phase-2 ECAL barrel multifit conditions (CATIA 2 gains,
// 16 samples). Analogous to EcalMultifitConditionsRcd, but depending on the
// Phase-2 condition records. Timing-related records will be added here once
// the GPU timing computation for Phase 2 is addressed.
class EcalMultifitConditionsPhase2Rcd
    : public edm::eventsetup::DependentRecordImplementation<EcalMultifitConditionsPhase2Rcd,
                                                            edm::mpl::Vector<EcalCATIAGainRatiosRcd,
                                                                             EcalLiteDTUPedestalsRcd,
                                                                             EcalPh2PulseCovariancesRcd,
                                                                             EcalPh2PulseShapesRcd,
                                                                             EcalPh2SamplesCorrelationRcd>> {};
#endif
