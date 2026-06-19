#include "CondFormats/DataRecord/interface/EcalMultifitConditionsRcd.h"
#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsDevice.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "DataFormats/Portable/interface/PortableObject.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/MoveToDeviceCache.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDGetToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"

#include "DeclsForKernels.h"
#include "EcalMultifitParameters.h"
#include "EcalUncalibRecHitMultiFitAlgoPhase2Portable.h"

#include <algorithm>

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace {
    using EcalMultifitParametersCache =
        cms::alpakatools::MoveToDeviceCache<Device, PortableHostObject<EcalMultifitParametersPhase2>>;
  }

  class EcalUncalibRecHitPhase2ProducerPortable
      : public stream::SynchronizingEDProducer<edm::GlobalCache<EcalMultifitParametersCache>> {
  public:
    explicit EcalUncalibRecHitPhase2ProducerPortable(edm::ParameterSet const& ps, EcalMultifitParametersCache const*);
    ~EcalUncalibRecHitPhase2ProducerPortable() override = default;
    static void fillDescriptions(edm::ConfigurationDescriptions&);
    static std::unique_ptr<EcalMultifitParametersCache> initializeGlobalCache(edm::ParameterSet const& ps);

    void acquire(device::Event const&, device::EventSetup const&) override;
    void produce(device::Event&, device::EventSetup const&) override;

    static void globalEndJob(EcalMultifitParametersCache*) {}

  private:
    using InputProduct = EcalDigiPhase2DeviceCollection;
    const device::EDGetToken<InputProduct> digisTokenEB_;
    using OutputProduct = EcalUncalibratedRecHitDeviceCollection;
    const device::EDPutToken<OutputProduct> uncalibRecHitsTokenEB_;

    // conditions tokens
    const device::ESGetToken<EcalMultifitConditionsDevice, EcalMultifitConditionsRcd> multifitConditionsToken_;

    // configuration parameters
    ecal::multifit::ConfigurationParametersPhase2 configParameters_;

    cms::alpakatools::host_buffer<uint32_t> ebDigisSizeHostBuf_;
  };

  void EcalUncalibRecHitPhase2ProducerPortable::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
    edm::ParameterSetDescription desc;

    desc.add<edm::InputTag>("digisLabelEB", edm::InputTag("ecalPhase2DigiToPortableProducer", "ebDigis"));

    desc.add<std::string>("recHitsLabelEB", "EcalUncalibRecHitsEB");

    desc.add<double>("EBtimeFitLimits_Lower", 0.2);
    desc.add<double>("EBtimeFitLimits_Upper", 1.4);
    desc.add<double>("EBtimeConstantTerm", .6);
    desc.add<double>("EBtimeNconst", 28.5);
    desc.add<double>("outOfTimeThresholdGain12pEB", 5);
    desc.add<double>("outOfTimeThresholdGain12mEB", 5);
    desc.add<double>("outOfTimeThresholdGain61pEB", 5);
    desc.add<double>("outOfTimeThresholdGain61mEB", 5);
    desc.add<double>("amplitudeThresholdEB", 10);

    desc.add<std::vector<double>>("EBtimeFitParameters",
                                  {-2.015452e+00,
                                   3.130702e+00,
                                   -1.234730e+01,
                                   4.188921e+01,
                                   -8.283944e+01,
                                   9.101147e+01,
                                   -5.035761e+01,
                                   1.105621e+01});
    desc.add<std::vector<double>>("EBamplitudeFitParameters", {1.138, 1.652});

    desc.addUntracked<std::vector<uint32_t>>("kernelMinimizeThreads", {32, 1, 1});
    desc.add<bool>("shouldRunTimingComputation", true);

    confDesc.addWithDefaultLabel(desc);
  }

  std::unique_ptr<EcalMultifitParametersCache> EcalUncalibRecHitPhase2ProducerPortable::initializeGlobalCache(
      edm::ParameterSet const& ps) {
    PortableHostObject<EcalMultifitParametersPhase2> params(cms::alpakatools::host());

    auto const ebTimeFitParamsFromPSet = ps.getParameter<std::vector<double>>("EBtimeFitParameters");
    // Assert that there are as many parameters as the EcalMultiFitParametersSoA expects
    assert(ebTimeFitParamsFromPSet.size() == EcalMultifitParametersPhase2::kNTimeFitParams);
    std::ranges::copy(ebTimeFitParamsFromPSet, params->timeFitParamsEB.begin());

    std::vector<float> ebAmplitudeFitParameters_;
    auto const ebAmplFitParamsFromPSet = ps.getParameter<std::vector<double>>("EBamplitudeFitParameters");
    // Assert that there are as many parameters as the EcalMultiFitParametersSoA expects
    assert(ebAmplFitParamsFromPSet.size() == EcalMultifitParametersPhase2::kNAmplitudeFitParams);
    std::ranges::copy(ebAmplFitParamsFromPSet, params->amplitudeFitParamsEB.begin());

    return std::make_unique<EcalMultifitParametersCache>(std::move(params));
  }

  EcalUncalibRecHitPhase2ProducerPortable::EcalUncalibRecHitPhase2ProducerPortable(const edm::ParameterSet& ps,
                                                                                   EcalMultifitParametersCache const*)
      : SynchronizingEDProducer(ps),
        digisTokenEB_{consumes(ps.getParameter<edm::InputTag>("digisLabelEB"))},
        uncalibRecHitsTokenEB_{produces(ps.getParameter<std::string>("recHitsLabelEB"))},
        multifitConditionsToken_{esConsumes()},
        ebDigisSizeHostBuf_{cms::alpakatools::make_host_buffer<uint32_t>()} {
    std::pair<double, double> EBtimeFitLimits, EEtimeFitLimits;
    EBtimeFitLimits.first = ps.getParameter<double>("EBtimeFitLimits_Lower");
    EBtimeFitLimits.second = ps.getParameter<double>("EBtimeFitLimits_Upper");

    auto EBtimeConstantTerm = ps.getParameter<double>("EBtimeConstantTerm");
    auto EBtimeNconst = ps.getParameter<double>("EBtimeNconst");

    auto outOfTimeThreshG12pEB = ps.getParameter<double>("outOfTimeThresholdGain12pEB");
    auto outOfTimeThreshG12mEB = ps.getParameter<double>("outOfTimeThresholdGain12mEB");
    auto outOfTimeThreshG61pEB = ps.getParameter<double>("outOfTimeThresholdGain61pEB");
    auto outOfTimeThreshG61mEB = ps.getParameter<double>("outOfTimeThresholdGain61mEB");
    auto amplitudeThreshEB = ps.getParameter<double>("amplitudeThresholdEB");

    // switch to run timing computation kernels
    configParameters_.shouldRunTimingComputation = ps.getParameter<bool>("shouldRunTimingComputation");

    // minimize kernel launch conf
    auto threadsMinimize = ps.getUntrackedParameter<std::vector<uint32_t>>("kernelMinimizeThreads");
    configParameters_.kernelMinimizeThreads[0] = threadsMinimize[0];
    configParameters_.kernelMinimizeThreads[1] = threadsMinimize[1];
    configParameters_.kernelMinimizeThreads[2] = threadsMinimize[2];

    //
    // configuration and physics parameters: done once
    // assume there is a single device
    // use sync copying
    //

    // time fit parameters and limits
    configParameters_.timeFitLimitsFirstEB = EBtimeFitLimits.first;
    configParameters_.timeFitLimitsSecondEB = EBtimeFitLimits.second;

    // time constant terms
    configParameters_.timeConstantTermEB = EBtimeConstantTerm;

    // time N const
    configParameters_.timeNconstEB = EBtimeNconst;

    // amplitude threshold for time flags
    configParameters_.amplitudeThreshEB = amplitudeThreshEB;

    // out of time thresholds gain-dependent
    configParameters_.outOfTimeThreshG12pEB = outOfTimeThreshG12pEB;
    configParameters_.outOfTimeThreshG61pEB = outOfTimeThreshG61pEB;
    configParameters_.outOfTimeThreshG12mEB = outOfTimeThreshG12mEB;
    configParameters_.outOfTimeThreshG61mEB = outOfTimeThreshG61mEB;
  }

  void EcalUncalibRecHitPhase2ProducerPortable::acquire(device::Event const& event, device::EventSetup const& setup) {
    auto& queue = event.queue();

    // get device collections from event
    auto const& ebDigisDev = event.get(digisTokenEB_);

    // copy the actual numbers of digis in the collections to host
    auto ebDigisSizeDevConstView =
        cms::alpakatools::make_device_view<const uint32_t>(queue, ebDigisDev.const_view().size());
    alpaka::memcpy(queue, ebDigisSizeHostBuf_, ebDigisSizeDevConstView);
  }

  void EcalUncalibRecHitPhase2ProducerPortable::produce(device::Event& event, device::EventSetup const& setup) {
    auto& queue = event.queue();

    // get device collections from event
    auto const& ebDigisDev = event.get(digisTokenEB_);

    // get the actual numbers of digis in the collections
    auto const ebDigisSize = static_cast<int>(*ebDigisSizeHostBuf_.data());

    // output device collections
    OutputProduct uncalibRecHitsDevEB{queue, ebDigisSize};
    // reset the size scalar of the SoA
    // memset takes an alpaka view that is created from the scalar in a view to the portable device collection
    auto uncalibRecHitSizeViewEB =
        cms::alpakatools::make_device_view<uint32_t>(queue, uncalibRecHitsDevEB.view().size());
    alpaka::memset(queue, uncalibRecHitSizeViewEB, 0);

    // stop here if there are no digis
    if (ebDigisSize > 0) {
      // conditions
      auto const& multifitConditionsDev = setup.getData(multifitConditionsToken_);
      auto const* multifitParametersDev = globalCache()->get(queue).const_data();

      //
      // schedule algorithms
      //
      ecal::multifit::launchKernels(
          queue, ebDigisDev, uncalibRecHitsDevEB, multifitConditionsDev, multifitParametersDev, configParameters_);
    }

    // put into the event
    event.emplace(uncalibRecHitsTokenEB_, std::move(uncalibRecHitsDevEB));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(EcalUncalibRecHitPhase2ProducerPortable);
