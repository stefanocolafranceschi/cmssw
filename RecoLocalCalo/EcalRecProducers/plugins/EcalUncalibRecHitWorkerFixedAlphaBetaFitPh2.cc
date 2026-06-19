#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <FWCore/ParameterSet/interface/ConfigurationDescriptions.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLiteDTUPedestalsRcd.h"

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitFixedAlphaBetaAlgo.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerRunOneDigiBase.h"

#include <cmath>
#include <fstream>

namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

class EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2 : public EcalUncalibRecHitWorkerRunOneDigiBase {
public:
  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2();
  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2(const edm::ParameterSet& ps, edm::ConsumesCollector&);
  ~EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2() override = default;

  void set(const edm::EventSetup& es) override;
  bool run(const edm::Event& evt,
           const edm::DataFrameContainer::const_iterator& digi,
           EcalUncalibratedRecHitCollection& result) override;

  edm::ParameterSetDescription getAlgoDescription() override;

private:
  EcalUncalibRecHitFixedAlphaBetaAlgo<EcalDataFrame_Ph2, ecalPh2> algo_;

  const double alpha_;
  const double beta_;
  const std::string alphaBetaFilename_;
  std::vector<std::vector<std::pair<double, double> > >
      alphaBetaValues_;  // List of alpha and Beta values [SM#][CRY#](alpha, beta)
  bool useAlphaBetaArray_;

  bool setAlphaBeta();  // Sets the alphaBetaValues_ vectors by the values provided in alphaBetaFilename_

  edm::ESHandle<EcalCATIAGainRatios> pRatio;
  edm::ESHandle<EcalLiteDTUPedestalsMap> pedHandle;
  edm::ESGetToken<EcalCATIAGainRatios, EcalCATIAGainRatiosRcd> ratiosToken_;
  edm::ESGetToken<EcalLiteDTUPedestalsMap, EcalLiteDTUPedestalsRcd> pedestalsToken_;
};

EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2() : alpha_(0.), beta_(0.) {}

EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2(const edm::ParameterSet& ps,
                                                                                         edm::ConsumesCollector& c)
    : EcalUncalibRecHitWorkerRunOneDigiBase(ps, c),
      alpha_(ps.getParameter<double>("alpha")),
      beta_(ps.getParameter<double>("beta")),
      alphaBetaFilename_(ps.getUntrackedParameter<std::string>("alphaBetaFilename")),
      ratiosToken_(c.esConsumes<EcalCATIAGainRatios, EcalCATIAGainRatiosRcd>()),
      pedestalsToken_(c.esConsumes<EcalLiteDTUPedestalsMap, EcalLiteDTUPedestalsRcd>()) {
  useAlphaBetaArray_ = setAlphaBeta();  // set crystalwise values of alpha and beta
  if (!useAlphaBetaArray_) {
    edm::LogInfo("EcalUncalibRecHitError") << " No alpha-beta file found. Using the default values.";
  }

  algo_.SetMinAmpl(ps.getParameter<double>("minAmpl"));
  algo_.SetDynamicPedestal(ps.getParameter<bool>("useDynamicPedestal"));
}

void EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::set(const edm::EventSetup& es) {
  // Gain Ratios
  LogDebug("EcalUncalibRecHitDebug") << "fetching gainRatios....";
  pRatio = es.getHandle(ratiosToken_);

  // fetch the pedestals from the cond DB via EventSetup
  LogDebug("EcalUncalibRecHitDebug") << "fetching pedestals....";
  pedHandle = es.getHandle(pedestalsToken_);
}

//Sets the alphaBetaValues_ vectors by the values provided in alphaBetaFilename_
bool EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::setAlphaBeta() {
  std::ifstream file(alphaBetaFilename_.c_str());
  if (!file.is_open())
    return false;

  alphaBetaValues_.resize(36);

  char buffer[100];
  int sm, cry, ret;
  float a, b;
  std::pair<double, double> p(-1, -1);

  while (!file.getline(buffer, 100).eof()) {
    ret = sscanf(buffer, "%d %d %f %f", &sm, &cry, &a, &b);
    if ((ret != 4) || (sm <= 0) || (sm > 36) || (cry <= 0) || (cry > 1700)) {
      // send warning
      edm::LogWarning("EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2")
          << "Incomplete input line with " << ret << "/4 values read or the channel in SM " << sm << ", crystal " << cry
          << " does not exist. Input line is ignored.";
      continue;
    }

    if (alphaBetaValues_[sm - 1].empty()) {
      alphaBetaValues_[sm - 1].resize(1700, p);
    }
    alphaBetaValues_[sm - 1][cry - 1].first = a;
    alphaBetaValues_[sm - 1][cry - 1].second = b;
  }

  file.close();
  return true;
}

bool EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::run(const edm::Event& evt,
                                                      const edm::DataFrameContainer::const_iterator& itdg,
                                                      EcalUncalibratedRecHitCollection& result) {
  const DetId detId(itdg->id());

  const EcalCATIAGainRatioMap& gainMap = pRatio.product()->getMap();      // map of gain ratios
  const EcalLiteDTUPedestalsMap& pedMap = pedHandle.product()->getMap();  // map of pedestals

  std::stringstream ss;

  // find pedestals for this channel
  const auto pedIter = pedMap.find(detId);
  if (pedIter == pedMap.end()) {
    ss << "Could not find pedestals for channel: ";
    detId.subdetId() == EcalBarrel ? ss << EBDetId(detId) : ss << EEDetId(detId);
    edm::LogWarning("EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2")
        << ss.str() << "\nNo uncalib rechit will be made for this digi.";
    return false;
  }
  double peds[ecalPh2::NGAINS];
  for (unsigned int i = 0; i < ecalPh2::NGAINS; ++i) {
    peds[i] = pedIter->mean(i);
  }

  // find gain ratios
  ss.str("");
  const auto gainIter = gainMap.find(detId);
  if (gainIter == gainMap.end()) {
    ss << "Could not find gain ratios for channel: ";
    detId.subdetId() == EcalBarrel ? ss << EBDetId(detId) : ss << EEDetId(detId);
    edm::LogWarning("EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2")
        << ss.str() << "\nNo uncalib rechit will be made for this digi.";
    return false;
  }
  double gainRatios[] = {1., *gainIter};  // The CATIA ASIC has two gains

  if (detId.subdetId() == EcalBarrel) {
    // Define alpha and beta either by stored values or by default universal values
    const EBDetId ebDetId(detId);
    double a = alpha_;
    double b = beta_;
    if (useAlphaBetaArray_) {
      if (!alphaBetaValues_[ebDetId.ism() - 1].empty()) {
        a = alphaBetaValues_[ebDetId.ism() - 1][ebDetId.ic() - 1].first;
        b = alphaBetaValues_[ebDetId.ism() - 1][ebDetId.ic() - 1].second;
        if (a == -1. && b == -1.) {
          a = alpha_;
          b = beta_;
        }
      }
    }
    algo_.SetAlphaBeta(a, b);
    result.push_back(algo_.makeRecHit(*itdg, peds, gainRatios, nullptr, nullptr));
  } else {
    edm::LogWarning("EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2")
        << "Channel " << EEDetId(detId) << " is not in the ECAL barrel. No uncalib rechit will be made for this digi.";
    return false;
  }
  return true;
}

edm::ParameterSetDescription EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2::getAlgoDescription() {
  edm::ParameterSetDescription psd;

  psd.addNode(edm::ParameterDescription<double>("alpha", 1.138, true) and
              edm::ParameterDescription<std::string>("alphaBetaFilename", "", false) and
              edm::ParameterDescription<double>("beta", 1.655, true) and
              edm::ParameterDescription<double>("minAmpl", 8.0, true) and
              edm::ParameterDescription<bool>("useDynamicPedestal", true, true));

  return psd;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitWorkerFactory,
                  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2,
                  "EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2");
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitFillDescriptionWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitFillDescriptionWorkerFactory,
                  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2,
                  "EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2");
