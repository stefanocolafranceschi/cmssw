// Trivial ESProducer for the Phase-2 ECAL barrel multifit condition records.
//
// None of these records has a tag in the production conditions DB yet (checked
// 150X_mcRun4_realistic_v1 and a DB-wide search, Jul 2026), so this module
// serves all five from configuration values, uniformly for all EB channels:
//   EcalLiteDTUPedestalsRcd      (pedestal mean/rms per CATIA gain)
//   EcalCATIAGainRatiosRcd       (gain10/gain1 ratio)
//   EcalPh2PulseShapesRcd        (16-sample pulse template)
//   EcalPh2PulseCovariancesRcd   (16x16 pulse-shape covariance)
//   EcalPh2SamplesCorrelationRcd (noise sample correlation per gain)
//
// The IOVs must be provided by EmptyESSource instances in the configuration
// (see ecalPhase2MultifitGPUCustomise.py). Remove once real tags exist.

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

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

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

class EcalPhase2TrivialCondESProducer : public edm::ESProducer {
public:
  explicit EcalPhase2TrivialCondESProducer(edm::ParameterSet const& ps) {
    pedestalsMean_ = ps.getParameter<std::vector<double>>("pedestalsMean");
    pedestalsRMS_ = ps.getParameter<std::vector<double>>("pedestalsRMS");
    gainRatio_ = ps.getParameter<double>("gain10Over1");
    pulseShape_ = ps.getParameter<std::vector<double>>("pulseShape");
    pulseCovariance_ = ps.getParameter<std::vector<double>>("pulseCovariance");
    sampleCorrelationG10_ = ps.getParameter<std::vector<double>>("sampleCorrelationG10");
    sampleCorrelationG1_ = ps.getParameter<std::vector<double>>("sampleCorrelationG1");

    constexpr size_t nGains = ecalPh2::NGAINS;
    constexpr size_t nT = EcalPh2PulseShape::TEMPLATESAMPLES;
    if (pedestalsMean_.size() != nGains || pedestalsRMS_.size() != nGains) {
      throw cms::Exception("Configuration") << "pedestalsMean/RMS must have " << nGains << " entries";
    }
    if (pulseShape_.size() != nT) {
      throw cms::Exception("Configuration") << "pulseShape must have " << nT << " entries";
    }
    if (pulseCovariance_.size() != nT * nT) {
      throw cms::Exception("Configuration") << "pulseCovariance must have " << nT * nT << " entries (row-major)";
    }
    if (sampleCorrelationG10_.size() != ecalPh2::sampleSize || sampleCorrelationG1_.size() != ecalPh2::sampleSize) {
      throw cms::Exception("Configuration") << "sampleCorrelation vectors must have " << ecalPh2::sampleSize
                                            << " entries";
    }

    setWhatProduced(this, &EcalPhase2TrivialCondESProducer::producePedestals);
    setWhatProduced(this, &EcalPhase2TrivialCondESProducer::produceGainRatios);
    setWhatProduced(this, &EcalPhase2TrivialCondESProducer::producePulseShapes);
    setWhatProduced(this, &EcalPhase2TrivialCondESProducer::producePulseCovariances);
    setWhatProduced(this, &EcalPhase2TrivialCondESProducer::produceSamplesCorrelation);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    // defaults: EcalLiteDTUPedestals class defaults and nominal CATIA ratio;
    // the pulse shape has NO meaningful default -- it must be configured
    // (e.g. with the template measured from the digis).
    desc.add<std::vector<double>>("pedestalsMean", {13., 8.});
    desc.add<std::vector<double>>("pedestalsRMS", {2.8, 1.2});
    desc.add<double>("gain10Over1", 10.);
    desc.add<std::vector<double>>("pulseShape", std::vector<double>(16, 0.));
    desc.add<std::vector<double>>("pulseCovariance", std::vector<double>(256, 0.));
    // delta correlation (white noise) placeholders; override with measured values
    auto delta = std::vector<double>(16, 0.);
    delta[0] = 1.;
    desc.add<std::vector<double>>("sampleCorrelationG10", delta);
    desc.add<std::vector<double>>("sampleCorrelationG1", delta);
    descriptions.addWithDefaultLabel(desc);
  }

  std::unique_ptr<EcalLiteDTUPedestalsMap> producePedestals(EcalLiteDTUPedestalsRcd const&) {
    auto product = std::make_unique<EcalLiteDTUPedestalsMap>();
    EcalLiteDTUPedestals item;
    for (unsigned int g = 0; g < ecalPh2::NGAINS; ++g) {
      item.setMean(g, pedestalsMean_[g]);
      item.setRMS(g, pedestalsRMS_[g]);
    }
    for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
      product->insert(std::make_pair(EBDetId::unhashIndex(i).rawId(), item));
    }
    return product;
  }

  std::unique_ptr<EcalCATIAGainRatios> produceGainRatios(EcalCATIAGainRatiosRcd const&) {
    auto product = std::make_unique<EcalCATIAGainRatios>();
    for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
      product->insert(std::make_pair(EBDetId::unhashIndex(i).rawId(), static_cast<EcalCATIAGainRatio>(gainRatio_)));
    }
    return product;
  }

  std::unique_ptr<EcalPh2PulseShapes> producePulseShapes(EcalPh2PulseShapesRcd const&) {
    auto product = std::make_unique<EcalPh2PulseShapes>();
    EcalPh2PulseShape shape;
    for (int t = 0; t < EcalPh2PulseShape::TEMPLATESAMPLES; ++t) {
      shape.pdfval[t] = pulseShape_[t];
    }
    for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
      product->insert(std::make_pair(EBDetId::unhashIndex(i).rawId(), shape));
    }
    return product;
  }

  std::unique_ptr<EcalPh2PulseCovariances> producePulseCovariances(EcalPh2PulseCovariancesRcd const&) {
    auto product = std::make_unique<EcalPh2PulseCovariances>();
    EcalPh2PulseCovariance cov;
    constexpr int nT = EcalPh2PulseCovariance::TEMPLATESAMPLES;
    for (int j = 0; j < nT; ++j) {
      for (int k = 0; k < nT; ++k) {
        cov.covval[j][k] = pulseCovariance_[j * nT + k];
      }
    }
    for (unsigned int i = 0; i < ecalPh2::kEBChannels; ++i) {
      product->insert(std::make_pair(EBDetId::unhashIndex(i).rawId(), cov));
    }
    return product;
  }

  std::unique_ptr<EcalPh2SamplesCorrelation> produceSamplesCorrelation(EcalPh2SamplesCorrelationRcd const&) {
    auto product = std::make_unique<EcalPh2SamplesCorrelation>();
    product->g10SamplesCorrelation = sampleCorrelationG10_;
    product->g1SamplesCorrelation = sampleCorrelationG1_;
    return product;
  }

private:
  std::vector<double> pedestalsMean_;
  std::vector<double> pedestalsRMS_;
  double gainRatio_;
  std::vector<double> pulseShape_;
  std::vector<double> pulseCovariance_;
  std::vector<double> sampleCorrelationG10_;
  std::vector<double> sampleCorrelationG1_;
};

DEFINE_FWK_EVENTSETUP_MODULE(EcalPhase2TrivialCondESProducer);
