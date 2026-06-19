// -*- C++ -*-
//
// Class:      EcalPh2PulseShapesESProducer
//
// Original Author:  Thomas Reis
//         Created:  Thu, 24 Feb 2022 17:29:40 GMT
//

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseShapesRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

class EcalPh2PulseShapesESProducer : public edm::ESProducer {
public:
  EcalPh2PulseShapesESProducer(const edm::ParameterSet&);
  ~EcalPh2PulseShapesESProducer() override = default;

  using ReturnType = std::unique_ptr<EcalPh2PulseShapes>;

  ReturnType produce(const EcalPh2PulseShapesRcd&);

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  std::vector<double> pulseShapes_;
};

EcalPh2PulseShapesESProducer::EcalPh2PulseShapesESProducer(const edm::ParameterSet& iConfig)
    : pulseShapes_(iConfig.getParameter<std::vector<double>>("pulseShapes")) {
  setWhatProduced(this);
}

EcalPh2PulseShapesESProducer::ReturnType EcalPh2PulseShapesESProducer::produce(const EcalPh2PulseShapesRcd& iRecord) {
  auto product = std::make_unique<EcalPh2PulseShapes>();

  const auto pulseShapesSize = pulseShapes_.size();
  if (pulseShapesSize != EcalPh2PulseShape::TEMPLATESAMPLES) {
    edm::LogWarning("EcalPh2PulseShapesESProducer")
        << "Input vector of pulse shapes has " << pulseShapesSize << " elements but the conditions format requires "
        << EcalPh2PulseShape::TEMPLATESAMPLES
        << ". Excess elements will be truncated and missing elements will be set to zero";
  }

  // set the same pulse shape for all channels
  for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
    EcalPh2PulseShape ps;
    const EBDetId ebDetId = EBDetId::unhashIndex(i);

    for (unsigned int j = 0; j < EcalPh2PulseShape::TEMPLATESAMPLES; ++j) {
      ps.pdfval[j] = (j < pulseShapesSize) ? static_cast<float>(pulseShapes_[j]) : 0.;
      product->insert(std::make_pair(ebDetId.rawId(), ps));
    }
  }

  return product;
}

void EcalPh2PulseShapesESProducer::fillDescriptions(edm::ConfigurationDescriptions& confDesc) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<double>>("pulseShapes", {0.});
  confDesc.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(EcalPh2PulseShapesESProducer);
