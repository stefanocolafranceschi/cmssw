// -*- C++ -*-
//
// Class:      EcalPh2PulseCovariancesESProducer
//
// Original Author:  Thomas Reis
//         Created:  Fri, 25 Feb 2022 11:37:30 GMT
//

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseCovariancesRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovarianceT.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

class EcalPh2PulseCovariancesESProducer : public edm::ESProducer {
public:
  EcalPh2PulseCovariancesESProducer(const edm::ParameterSet&);
  ~EcalPh2PulseCovariancesESProducer() override = default;

  using ReturnType = std::unique_ptr<EcalPh2PulseCovariances>;

  ReturnType produce(const EcalPh2PulseCovariancesRcd&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  // stores the covarinace matrix as a vector
  std::vector<double> pulseCov_;
};

EcalPh2PulseCovariancesESProducer::EcalPh2PulseCovariancesESProducer(const edm::ParameterSet& iConfig)
    : pulseCov_(iConfig.getParameter<std::vector<double>>("pulseCovariances")) {
  setWhatProduced(this);
}

EcalPh2PulseCovariancesESProducer::ReturnType EcalPh2PulseCovariancesESProducer::produce(
    const EcalPh2PulseCovariancesRcd& iRecord) {
  auto product = std::make_unique<EcalPh2PulseCovariances>();

  const auto pulseCovSize = pulseCov_.size();
  const auto nTempSamples = EcalPh2PulseCovariance::TEMPLATESAMPLES;
  const auto covMatrixElemSize = nTempSamples * nTempSamples;

  if (pulseCovSize != covMatrixElemSize) {
    edm::LogWarning("EcalPh2PulseCovariancesESProducer")
        << "Input vector of pulse covariances has " << pulseCovSize << " elements but the conditions format requires "
        << covMatrixElemSize << " elements to form a " << nTempSamples << "x" << nTempSamples
        << "matrix. Excess elements will be truncated and missing elements will be set to zero";
  }

  // set the same pulse covariances for all channels
  for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
    EcalPh2PulseCovariance pc;
    const EBDetId ebDetId = EBDetId::unhashIndex(i);

    size_t col = 0;
    size_t row = 0;
    for (unsigned int j = 0; j < covMatrixElemSize; ++j) {
      pc.covval[row][col] = (j < pulseCovSize) ? static_cast<float>(pulseCov_[j]) : 0.;
      product->insert(std::make_pair(ebDetId.rawId(), pc));

      // move to next element in matrix
      if (j != nTempSamples - 1) {
        ++col;
      } else {
        col = 0;
        ++row;
      }
    }
  }

  return product;
}

void EcalPh2PulseCovariancesESProducer::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<double>>("pulseCovariances", {0.});
  confDesc.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(EcalPh2PulseCovariancesESProducer);
