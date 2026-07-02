#include <algorithm>
#include <cassert>
#include <cstring>

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLiteDTUPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseCovariancesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2SamplesCorrelationRcd.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalPh2SamplesCorrelation.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovarianceT.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"

#include "CondFormats/DataRecord/interface/EcalMultifitConditionsPhase2Rcd.h"
#include "CondFormats/EcalObjects/interface/EcalMultifitConditionsPhase2SoA.h"
#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPhase2Device.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ESProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/ModuleFactory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  // ESProducer for the Phase-2 ECAL barrel multifit conditions.
  // Fills the 16-sample, 2-CATIA-gain conditions SoA from the Phase-2
  // condition records. Mirrors EcalMultifitConditionsHostESProducer (Run 3).
  // Timing-related conditions are not filled yet (Phase-2 GPU timing deferred).
  class EcalMultifitConditionsPhase2HostESProducer : public ESProducer {
  public:
    EcalMultifitConditionsPhase2HostESProducer(edm::ParameterSet const& iConfig) : ESProducer(iConfig) {
      auto cc = setWhatProduced(this);
      pedestalsToken_ = cc.consumes();
      gainRatiosToken_ = cc.consumes();
      pulseShapesToken_ = cc.consumes();
      pulseCovariancesToken_ = cc.consumes();
      samplesCorrelationToken_ = cc.consumes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      descriptions.addWithDefaultLabel(desc);
    }

    std::unique_ptr<EcalMultifitConditionsPhase2Host> produce(EcalMultifitConditionsPhase2Rcd const& iRecord) {
      auto const& pedestalsData = iRecord.get(pedestalsToken_);
      auto const& gainRatiosData = iRecord.get(gainRatiosToken_);
      auto const& pulseShapesData = iRecord.get(pulseShapesToken_);
      auto const& pulseCovariancesData = iRecord.get(pulseCovariancesToken_);
      auto const& samplesCorrelationData = iRecord.get(samplesCorrelationToken_);

      // Phase-2 ECAL is barrel-only
      auto const& pedestalsEB = pedestalsData.barrelItems();
      auto const& gainRatiosEB = gainRatiosData.barrelItems();
      auto const& pulseShapesEB = pulseShapesData.barrelItems();
      auto const& pulseCovariancesEB = pulseCovariancesData.barrelItems();

      size_t const numberOfXtals = pedestalsEB.size();  // ecalPh2::kEBChannels
      assert(numberOfXtals > 0);
      assert(gainRatiosEB.size() == numberOfXtals);
      assert(pulseShapesEB.size() == numberOfXtals);
      assert(pulseCovariancesEB.size() == numberOfXtals);

      auto product = std::make_unique<EcalMultifitConditionsPhase2Host>(cms::alpakatools::host(), numberOfXtals);
      auto view = product->view();

      for (size_t i = 0; i < numberOfXtals; ++i) {
        auto vi = view[i];

        vi.rawid() = EBDetId::unhashIndex(i).rawId();

        vi.pedestals_mean_g10() = pedestalsEB[i].mean(ecalPh2::gainId10);
        vi.pedestals_rms_g10() = pedestalsEB[i].rms(ecalPh2::gainId10);
        vi.pedestals_mean_g1() = pedestalsEB[i].mean(ecalPh2::gainId1);
        vi.pedestals_rms_g1() = pedestalsEB[i].rms(ecalPh2::gainId1);

        vi.gain10Over1() = gainRatiosEB[i];

        std::memcpy(
            vi.pulseShapes().data(), pulseShapesEB[i].pdfval, sizeof(float) * EcalPh2PulseShape::TEMPLATESAMPLES);

        for (int j = 0; j < EcalPh2PulseCovariance::TEMPLATESAMPLES; ++j) {
          for (int k = 0; k < EcalPh2PulseCovariance::TEMPLATESAMPLES; ++k) {
            vi.pulseCovariance()(j, k) = pulseCovariancesEB[i].val(j, k);
          }
        }
      }  // end barrel loop

      // === Scalar data (not per xtal)
      // Sample correlations per CATIA gain
      assert(samplesCorrelationData.g10SamplesCorrelation.size() >= ecalPh2::sampleSize);
      assert(samplesCorrelationData.g1SamplesCorrelation.size() >= ecalPh2::sampleSize);
      std::memcpy(view.sampleCorrelation_g10().data(),
                  samplesCorrelationData.g10SamplesCorrelation.data(),
                  sizeof(double) * ecalPh2::sampleSize);
      std::memcpy(view.sampleCorrelation_g1().data(),
                  samplesCorrelationData.g1SamplesCorrelation.data(),
                  sizeof(double) * ecalPh2::sampleSize);

      return product;
    }

  private:
    edm::ESGetToken<EcalLiteDTUPedestalsMap, EcalLiteDTUPedestalsRcd> pedestalsToken_;
    edm::ESGetToken<EcalCATIAGainRatios, EcalCATIAGainRatiosRcd> gainRatiosToken_;
    edm::ESGetToken<EcalPh2PulseShapes, EcalPh2PulseShapesRcd> pulseShapesToken_;
    edm::ESGetToken<EcalPh2PulseCovariances, EcalPh2PulseCovariancesRcd> pulseCovariancesToken_;
    edm::ESGetToken<EcalPh2SamplesCorrelation, EcalPh2SamplesCorrelationRcd> samplesCorrelationToken_;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_EVENTSETUP_ALPAKA_MODULE(EcalMultifitConditionsPhase2HostESProducer);
