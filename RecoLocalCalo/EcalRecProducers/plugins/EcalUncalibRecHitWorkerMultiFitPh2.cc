#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/EmptyGroupDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSwitch.h"

#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"
#include "CondFormats/EcalObjects/interface/EcalLiteDTUPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalCATIAGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalWeightXtalGroups.h"
#include "CondFormats/EcalObjects/interface/EcalTBWeightsT.h"
#include "CondFormats/EcalObjects/interface/EcalSampleMask.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"
#include "CondFormats/EcalObjects/interface/EcalPh2SamplesCorrelation.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapeT.h"
#include "CondFormats/EcalObjects/interface/EcalPulseCovarianceT.h"
#include "CondFormats/DataRecord/interface/EcalCATIAGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalLiteDTUPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/EcalWeightXtalGroupsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2TBWeightsRcd.h"
#include "CondFormats/DataRecord/interface/EcalSampleMaskRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2SamplesCorrelationRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseShapesRcd.h"
#include "CondFormats/DataRecord/interface/EcalPh2PulseCovariancesRcd.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitTimingCCAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgoPh2.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitTimeWeightsAlgoPh2.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecChi2Algo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRatioMethodAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerBaseClass.h"

#include <cmath>

class EcalUncalibRecHitWorkerMultiFitPh2 final : public EcalUncalibRecHitWorkerBaseClass {
public:
  EcalUncalibRecHitWorkerMultiFitPh2();
  EcalUncalibRecHitWorkerMultiFitPh2(const edm::ParameterSet&, edm::ConsumesCollector& c);
  ~EcalUncalibRecHitWorkerMultiFitPh2() override{};

  edm::ParameterSetDescription getAlgoDescription() override;

private:
  using SampleVector = typename EigenMatrixTypes<ecalPh2>::SampleVector;
  using FullSampleVector = typename EigenMatrixTypes<ecalPh2>::FullSampleVector;
  using SampleMatrix = typename EigenMatrixTypes<ecalPh2>::SampleMatrix;
  using FullSampleMatrix = typename EigenMatrixTypes<ecalPh2>::FullSampleMatrix;
  using SampleMatrixGainArray = typename EigenMatrixTypes<ecalPh2>::SampleMatrixGainArray;
  using BXVector = typename EigenMatrixTypes<ecalPh2>::BXVector;

  void set(const edm::EventSetup& es) override;
  void set(const edm::Event& evt) override;
  void run(const edm::Event& evt,
           const edm::DataFrameContainer& digis,
           EcalUncalibratedRecHitCollection& result) override;
  double timeCorrection(const float ampli,
                        const std::vector<float>& amplitudeBins,
                        const std::vector<float>& shiftBins) const;

  edm::ESHandle<EcalLiteDTUPedestalsMap> peds_;
  edm::ESGetToken<EcalLiteDTUPedestalsMap, EcalLiteDTUPedestalsRcd> pedsToken_;
  edm::ESHandle<EcalCATIAGainRatios> gains_;
  edm::ESGetToken<EcalCATIAGainRatios, EcalCATIAGainRatiosRcd> gainsToken_;
  edm::ESGetToken<EcalPh2SamplesCorrelation, EcalPh2SamplesCorrelationRcd> noiseConvariancesToken_;
  edm::ESHandle<EcalPh2PulseShapes> pulseshapes_;
  edm::ESGetToken<EcalPh2PulseShapes, EcalPh2PulseShapesRcd> pulseShapesToken_;
  edm::ESHandle<EcalPh2PulseCovariances> pulsecovariances_;
  edm::ESGetToken<EcalPh2PulseCovariances, EcalPh2PulseCovariancesRcd> pulseConvariancesToken_;

  // multifit method
  EcalUncalibRecHitMultiFitAlgoPh2 multiFitMethod_;
  SampleMatrixGainArray noisecors_;
  BXVector activeBX_;
  // uncertainty calculation (CPU intensive)
  const bool ampErrorCalculation_;
  const bool useLumiInfoRunHeader_;

  int bunchSpacingManual_;
  edm::EDGetTokenT<unsigned int> bunchSpacing_;

  const bool doPrefit_;
  const double prefitMaxChiSq_;
  const bool dynamicPedestals_;
  const bool mitigateBadSamples_;
  const bool gainSwitchUseMaxSample_;
  const bool selectiveBadSampleCriteria_;
  const double addPedestalUncertainty_;
  const bool simplifiedNoiseModelForGainSwitch_;

  // time algorithm to be used to set the jitter and its uncertainty
  enum TimeAlgo { noMethod, ratioMethod, weightsMethod, crossCorrelationMethod };
  TimeAlgo timealgo_ = noMethod;

  // ratio method method for timing
  std::unique_ptr<EcalUncalibRecHitRatioMethodAlgo<EcalDataFrame_Ph2, EcalPh2SampleMask>> ratioMethodAlgo_;
  edm::ESHandle<EcalTimeBiasCorrections> timeCorrBias_;
  edm::ESGetToken<EcalTimeBiasCorrections, EcalTimeBiasCorrectionsRcd> timeCorrBiasToken_;
  edm::ESHandle<EcalTimeCalibConstants> itime_;
  edm::ESGetToken<EcalTimeCalibConstants, EcalTimeCalibConstantsRcd> itimeToken_;
  edm::ESHandle<EcalTimeOffsetConstant> offtime_;
  edm::ESGetToken<EcalTimeOffsetConstant, EcalTimeOffsetConstantRcd> offtimeToken_;
  edm::ESHandle<EcalPh2SampleMask> sampleMaskHand_;
  edm::ESGetToken<EcalPh2SampleMask, EcalSampleMaskRcd> sampleMaskToken_;
  std::vector<double> timeFitParameters_;
  std::vector<double> amplitudeFitParameters_;
  std::pair<double, double> timeFitLimits_;
  double timeConstantTerm_;
  double timeNconst_;
  double outOfTimeThreshG10p_;
  double outOfTimeThreshG10m_;
  double outOfTimeThreshG1p_;
  double outOfTimeThreshG1m_;
  double amplitudeThresh_;

  // time weights method
  std::unique_ptr<EcalUncalibRecHitTimeWeightsAlgoPh2> weightsMethodAlgo_;
  edm::ESHandle<EcalWeightXtalGroups> grps_;
  edm::ESGetToken<EcalWeightXtalGroups, EcalWeightXtalGroupsRcd> grpsToken_;
  edm::ESHandle<EcalPh2TBWeights> wgts_;
  edm::ESGetToken<EcalPh2TBWeights, EcalPh2TBWeightsRcd> wgtsToken_;
  const EcalPh2WeightSet::EcalWeightMatrix* weights_[2];

  //Timing Cross Correlation Algo
  std::unique_ptr<EcalUncalibRecHitTimingCCAlgo> computeCC_;
  double CCminTimeToBeLateMin_;
  double CCminTimeToBeLateMax_;
  double CCTimeShiftWrtRations_;
  double CCtargetTimePrecision_;
  double CCtargetTimePrecisionForDelayedPulses_;
};

EcalUncalibRecHitWorkerMultiFitPh2::EcalUncalibRecHitWorkerMultiFitPh2()
    : EcalUncalibRecHitWorkerBaseClass(),
      ampErrorCalculation_(false),
      useLumiInfoRunHeader_(false),
      doPrefit_(false),
      prefitMaxChiSq_(0.),
      dynamicPedestals_(false),
      mitigateBadSamples_(false),
      gainSwitchUseMaxSample_(false),
      selectiveBadSampleCriteria_(false),
      addPedestalUncertainty_(0.),
      simplifiedNoiseModelForGainSwitch_(false) {}

EcalUncalibRecHitWorkerMultiFitPh2::EcalUncalibRecHitWorkerMultiFitPh2(const edm::ParameterSet& ps,
                                                                       edm::ConsumesCollector& c)
    : EcalUncalibRecHitWorkerBaseClass(ps, c),
      ampErrorCalculation_(ps.getParameter<bool>("ampErrorCalculation")),
      useLumiInfoRunHeader_(ps.getParameter<bool>("useLumiInfoRunHeader")),
      doPrefit_(ps.getParameter<bool>("doPrefit")),
      prefitMaxChiSq_(ps.getParameter<double>("prefitMaxChiSq")),
      dynamicPedestals_(ps.getParameter<bool>("dynamicPedestals")),
      mitigateBadSamples_(ps.getParameter<bool>("mitigateBadSamples")),
      gainSwitchUseMaxSample_(ps.getParameter<bool>("gainSwitchUseMaxSample")),
      selectiveBadSampleCriteria_(ps.getParameter<bool>("selectiveBadSampleCriteria")),
      addPedestalUncertainty_(ps.getParameter<double>("addPedestalUncertainty")),
      simplifiedNoiseModelForGainSwitch_(ps.getParameter<bool>("simplifiedNoiseModelForGainSwitch")) {
  // get the BX for the pulses to be activated
  std::vector<int32_t> activeBXs = ps.getParameter<std::vector<int32_t>>("activeBXs");
  activeBX_.resize(activeBXs.size());
  for (unsigned int ibx = 0; ibx < activeBXs.size(); ++ibx) {
    activeBX_.coeffRef(ibx) = activeBXs[ibx];
  }

  if (useLumiInfoRunHeader_) {
    bunchSpacing_ = c.consumes<unsigned int>(edm::InputTag("bunchSpacingProducer"));
    bunchSpacingManual_ = 0;
  } else {
    bunchSpacingManual_ = ps.getParameter<int>("bunchSpacing");
  }

  pedsToken_ = c.esConsumes<EcalLiteDTUPedestalsMap, EcalLiteDTUPedestalsRcd>();
  gainsToken_ = c.esConsumes<EcalCATIAGainRatios, EcalCATIAGainRatiosRcd>();
  noiseConvariancesToken_ = c.esConsumes<EcalPh2SamplesCorrelation, EcalPh2SamplesCorrelationRcd>();
  pulseShapesToken_ = c.esConsumes<EcalPh2PulseShapes, EcalPh2PulseShapesRcd>();
  pulseConvariancesToken_ = c.esConsumes<EcalPh2PulseCovariances, EcalPh2PulseCovariancesRcd>();

  // algorithm to be used for timing and required configurations
  auto const& timeAlgoName = ps.getParameter<std::string>("timealgo");
  if (timeAlgoName == "RatioMethod") {
    timealgo_ = ratioMethod;
    itimeToken_ = c.esConsumes<EcalTimeCalibConstants, EcalTimeCalibConstantsRcd>();
    offtimeToken_ = c.esConsumes<EcalTimeOffsetConstant, EcalTimeOffsetConstantRcd>();
    sampleMaskToken_ = c.esConsumes<EcalPh2SampleMask, EcalSampleMaskRcd>();
    timeCorrBiasToken_ = c.esConsumes<EcalTimeBiasCorrections, EcalTimeBiasCorrectionsRcd>();
    ratioMethodAlgo_ = std::make_unique<EcalUncalibRecHitRatioMethodAlgo<EcalDataFrame_Ph2, EcalPh2SampleMask>>();
    timeFitParameters_ = ps.getParameter<std::vector<double>>("timeFitParameters");
    amplitudeFitParameters_ = ps.getParameter<std::vector<double>>("amplitudeFitParameters");
    timeFitLimits_.first = ps.getParameter<double>("timeFitLimits_Lower");
    timeFitLimits_.second = ps.getParameter<double>("timeFitLimits_Upper");
    timeConstantTerm_ = ps.getParameter<double>("timeConstantTerm");
    timeNconst_ = ps.getParameter<double>("timeNconst");
    outOfTimeThreshG10p_ = ps.getParameter<double>("outOfTimeThresholdGain10p");
    outOfTimeThreshG10m_ = ps.getParameter<double>("outOfTimeThresholdGain10m");
    outOfTimeThreshG1p_ = ps.getParameter<double>("outOfTimeThresholdGain1p");
    outOfTimeThreshG1m_ = ps.getParameter<double>("outOfTimeThresholdGain1m");
    amplitudeThresh_ = ps.getParameter<double>("amplitudeThreshold");
  } else if (timeAlgoName == "WeightsMethod") {
    timealgo_ = weightsMethod;
    grpsToken_ = c.esConsumes<EcalWeightXtalGroups, EcalWeightXtalGroupsRcd>();
    wgtsToken_ = c.esConsumes<EcalPh2TBWeights, EcalPh2TBWeightsRcd>();
    weightsMethodAlgo_ = std::make_unique<EcalUncalibRecHitTimeWeightsAlgoPh2>();
  } else if (timeAlgoName == "crossCorrelationMethod") {
    timealgo_ = crossCorrelationMethod;
    const auto startTime = ps.getParameter<double>("crossCorrelationStartTime");
    const auto stopTime = ps.getParameter<double>("crossCorrelationStopTime");
    CCtargetTimePrecision_ = ps.getParameter<double>("crossCorrelationTargetTimePrecision");
    CCtargetTimePrecisionForDelayedPulses_ =
        ps.getParameter<double>("crossCorrelationTargetTimePrecisionForDelayedPulses");
    CCminTimeToBeLateMin_ = ps.getParameter<double>("crossCorrelationMinTimeToBeLateMin") / ecalcctiming::clockToNS;
    CCminTimeToBeLateMax_ = ps.getParameter<double>("crossCorrelationMinTimeToBeLateMax") / ecalcctiming::clockToNS;
    CCTimeShiftWrtRations_ = ps.getParameter<double>("crossCorrelationTimeShiftWrtRations");
    computeCC_ = std::make_unique<EcalUncalibRecHitTimingCCAlgo>(startTime, stopTime);
  } else if (timeAlgoName != "None") {
    edm::LogInfo("EcalUncalibRecHit") << "No time estimation algorithm requested";
  } else {
    edm::LogWarning("EcalUncalibRecHit") << "Unknown time estimation algorithm '" << timeAlgoName << "'";
  }
}

void EcalUncalibRecHitWorkerMultiFitPh2::set(const edm::EventSetup& es) {
  // common setup
  gains_ = es.getHandle(gainsToken_);
  peds_ = es.getHandle(pedsToken_);

  // for the multifit method
  if (!ampErrorCalculation_)
    multiFitMethod_.disableErrorCalculation();
  const auto& noisecovariances = es.getData(noiseConvariancesToken_);
  pulseshapes_ = es.getHandle(pulseShapesToken_);
  pulsecovariances_ = es.getHandle(pulseConvariancesToken_);

  // parameters for the time reconstruction
  if (timealgo_ == ratioMethod) {
    itime_ = es.getHandle(itimeToken_);
    offtime_ = es.getHandle(offtimeToken_);
    // which of the samples need be used
    sampleMaskHand_ = es.getHandle(sampleMaskToken_);
    // for the time correction methods
    timeCorrBias_ = es.getHandle(timeCorrBiasToken_);
  } else if (timealgo_ == weightsMethod) {
    grps_ = es.getHandle(grpsToken_);
    wgts_ = es.getHandle(wgtsToken_);
  }

  const int nnoise = SampleVector::RowsAtCompileTime;
  auto& noisecorg10 = noisecors_[ecalPh2::gainId10];
  auto& noisecorg1 = noisecors_[ecalPh2::gainId1];

  for (int i = 0; i < nnoise; ++i) {
    for (int j = 0; j < nnoise; ++j) {
      int vidx = std::abs(j - i);
      noisecorg10(i, j) = noisecovariances.g10SamplesCorrelation[vidx];
      noisecorg1(i, j) = noisecovariances.g1SamplesCorrelation[vidx];
    }
  }
}

void EcalUncalibRecHitWorkerMultiFitPh2::set(const edm::Event& evt) {
  unsigned int bunchspacing = 450;

  if (useLumiInfoRunHeader_) {
    bunchspacing = evt.get(bunchSpacing_);
  } else {
    bunchspacing = bunchSpacingManual_;
  }

  if (useLumiInfoRunHeader_ || bunchSpacingManual_ > 0) {
    if (bunchspacing == 25) {
      activeBX_.resize(ecalPh2::kNActiveLHCPeriods);
      activeBX_ << -2, -1, 0, 1, 2;
    } else {
      //50ns configuration otherwise (also for no pileup)
      activeBX_.resize(static_cast<unsigned int>(std::round(ecalPh2::kNActiveLHCPeriods / 2.)));
      activeBX_ << -2, 0, 2;
    }
  }
}

/**
 * Amplitude-dependent time corrections; EB have separate corrections:
 * EXtimeCorrAmplitudes (ADC) and EXtimeCorrShifts (ns) need to have the same number of elements
 * Bins must be ordered in amplitude. First-last bins take care of under-overflows.
 *
 * @return Jitter (in clock cycles) which will be added to UncalibRechit.setJitter(), 0 if no correction is applied.
 */
double EcalUncalibRecHitWorkerMultiFitPh2::timeCorrection(const float ampli,
                                                          const std::vector<float>& amplitudeBins,
                                                          const std::vector<float>& shiftBins) const {
  // computed initially in ns. Than turned in the BX's, as
  // EcalUncalibratedRecHit need be.
  double theCorrection = 0;

  // sanity check for arrays
  if (amplitudeBins.empty()) {
    edm::LogWarning("EcalUncalibRecHitProducer") << "timeCorrAmplitudeBins is empty, forcing no time bias corrections.";
    return 0;
  }

  if (amplitudeBins.size() != shiftBins.size()) {
    edm::LogWarning("EcalUncalibRecHitProducer")
        << "Size of timeCorrAmplitudeBins different from timeCorrShiftBins. Forcing no time bias corrections. ";
    return 0;
  }

  int myBin = -1;
  for (int bin = 0; bin < (int)amplitudeBins.size(); ++bin) {
    if (ampli > amplitudeBins[bin]) {
      myBin = bin;
    } else {
      break;
    }
  }

  if (myBin == -1) {
    theCorrection = shiftBins[0];
  } else if (myBin == ((int)(amplitudeBins.size() - 1))) {
    theCorrection = shiftBins[myBin];
  } else {
    // interpolate linearly between two assingned points
    theCorrection = (shiftBins[myBin + 1] - shiftBins[myBin]);
    theCorrection *= (((double)ampli) - amplitudeBins[myBin]) / (amplitudeBins[myBin + 1] - amplitudeBins[myBin]);
    theCorrection += shiftBins[myBin];
  }

  // convert ns into clocks
  constexpr double invPeriod = 1. / ecalPh2::Samp_Period;
  return theCorrection * invPeriod;
}

void EcalUncalibRecHitWorkerMultiFitPh2::run(const edm::Event& evt,
                                             const edm::DataFrameContainer& digis,
                                             EcalUncalibratedRecHitCollection& result) {
  if (digis.empty())
    return;

  const DetId detid(digis.begin()->id());
  const bool barrel = detid.subdetId() == EcalBarrel;
  if (!barrel) {
    return;
  }

  multiFitMethod_.setSimplifiedNoiseModelForGainSwitch(simplifiedNoiseModelForGainSwitch_);
  multiFitMethod_.setDoPrefit(doPrefit_);
  multiFitMethod_.setPrefitMaxChiSq(prefitMaxChiSq_);
  multiFitMethod_.setDynamicPedestals(dynamicPedestals_);
  multiFitMethod_.setMitigateBadSamples(mitigateBadSamples_);
  multiFitMethod_.setGainSwitchUseMaxSample(gainSwitchUseMaxSample_);
  multiFitMethod_.setSelectiveBadSampleCriteria(selectiveBadSampleCriteria_);
  multiFitMethod_.setAddPedestalUncertainty(addPedestalUncertainty_);

  FullSampleVector fullpulse(FullSampleVector::Zero());
  FullSampleMatrix fullpulsecov(FullSampleMatrix::Zero());

  assert(ecalPh2::kFullSampleVectorSize >= ecalPh2::kPulseShapeTemplateSampleSize);
  const auto indexOffset = ecalPh2::kNSamplesPerLHCPeriod * (ecalPh2::kNActiveLHCPeriods / 2);
  result.reserve(result.size() + digis.size());
  for (const auto digi : digis) {
    const auto& df = static_cast<EcalDataFrame_Ph2>(digi);
    const DetId detid(df.id());

    const unsigned int hashedIndex = EBDetId(detid).hashedIndex();
    const auto* aped = &peds_->barrel(hashedIndex);
    const auto* aGain = &gains_->barrel(hashedIndex);
    const auto* aPulse = &pulseshapes_->barrel(hashedIndex);
    const auto* aPulseCov = &pulsecovariances_->barrel(hashedIndex);

    double pedVec[ecalPh2::NGAINS];
    double pedRMSVec[ecalPh2::NGAINS];
    const double gainRatios[ecalPh2::NGAINS] = {1., *aGain};
    for (unsigned int i = 0; i < ecalPh2::NGAINS; ++i) {
      pedVec[i] = aped->mean(i);
      pedRMSVec[i] = aped->rms(i);
    }

    for (int i = 0; i < EcalPh2PulseShape::TEMPLATESAMPLES; ++i)
      fullpulse(i + indexOffset) = aPulse->val(i);

    for (int i = 0; i < EcalPh2PulseShape::TEMPLATESAMPLES; ++i)
      for (int j = 0; j < EcalPh2PulseShape::TEMPLATESAMPLES; ++j)
        fullpulsecov(i + indexOffset, j + indexOffset) = aPulseCov->val(i, j);

    int lastSampleBeforeSaturation = -2;
    for (unsigned int iSample = 0; iSample < EcalDataFrame_Ph2::MAXSAMPLES; ++iSample) {
      if (df.sample(iSample).gainId() == ecalPh2::gainId1 && df.sample(iSample).adc() == ecalPh2::MAXADC) {
        lastSampleBeforeSaturation = iSample - 1;
        break;
      }
    }

    // === amplitude computation ===

    if (lastSampleBeforeSaturation == ecalPh2::kMaxSampleIdx - 1) {  // saturation on the expected max sample
      result.emplace_back(df.id(), ecalPh2::MAXADC * ecalPh2::gains[ecalPh2::gainId10], 0, 0, 0);
      auto& uncalibRecHit = result.back();
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kSaturated);
      // do not propagate the default chi2 = -1 value to the calib rechit (mapped to 64), set it to 0 when saturation
      uncalibRecHit.setChi2(0);
    } else if (lastSampleBeforeSaturation >= -1) {  // saturation on other samples: cannot extrapolate from fourth one
      const auto gainId = df.sample(ecalPh2::kMaxSampleIdx).gainId();
      const auto pedestal = pedVec[gainId];
      const auto gainratio = gainRatios[gainId];
      const double amplitude = (static_cast<double>(df.sample(ecalPh2::kMaxSampleIdx).adc()) - pedestal) * gainratio;
      result.emplace_back(df.id(), amplitude, 0, 0, 0);
      auto& uncalibRecHit = result.back();
      uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kSaturated);
      // do not propagate the default chi2 = -1 value to the calib rechit (mapped to 64), set it to 0 when saturation
      uncalibRecHit.setChi2(0);
    } else {
      // multifit
      result.push_back(multiFitMethod_.makeRecHit(df, aped, aGain, noisecors_, fullpulse, fullpulsecov, activeBX_));
      auto& uncalibRecHit = result.back();

      // === time computation ===
      if (timealgo_ == ratioMethod) {
        // ratio method
        constexpr float clockToNsConstant = ecalPh2::Samp_Period;
        const EcalPh2SampleMask* sampleMask = sampleMaskHand_.product();
        ratioMethodAlgo_->init(df, *sampleMask, pedVec, pedRMSVec, gainRatios);
        ratioMethodAlgo_->fixMGPAslew(df);
        ratioMethodAlgo_->computeTime(timeFitParameters_, timeFitLimits_, amplitudeFitParameters_);
        ratioMethodAlgo_->computeAmplitude(amplitudeFitParameters_);
        const EcalUncalibRecHitRatioMethodAlgo<EcalDataFrame_Ph2, EcalPh2SampleMask>::CalculatedRecHit crh =
            ratioMethodAlgo_->getCalculatedRecHit();

        const double theTimeCorrection = timeCorrection(
            uncalibRecHit.amplitude(), timeCorrBias_->EBTimeCorrAmplitudeBins, timeCorrBias_->EBTimeCorrShiftBins);

        uncalibRecHit.setJitter(crh.timeMax - 5 + theTimeCorrection);
        uncalibRecHit.setJitterError(std::hypot(crh.timeError, timeConstantTerm_ / clockToNsConstant));

        // consider flagging as kOutOfTime only if above noise
        if (uncalibRecHit.amplitude() > pedRMSVec[0] * amplitudeThresh_) {
          float outOfTimeThreshP = outOfTimeThreshG10p_;
          float outOfTimeThreshM = outOfTimeThreshG10m_;
          // determine if gain has switched away from gainId==0 (x10 gain)
          // and determine cuts (number of 'sigmas') to ose for kOutOfTime
          // >3k ADC is necessasry condition for gain switch to occur
          if (uncalibRecHit.amplitude() > 3000.) {
            for (int iSample = 0; iSample < EcalDataFrame_Ph2::MAXSAMPLES; ++iSample) {
              int GainId = df.sample(iSample).gainId();
              if (GainId != ecalPh2::gainId1) {
                outOfTimeThreshP = outOfTimeThreshG1p_;
                outOfTimeThreshM = outOfTimeThreshG1m_;
                break;
              }
            }
          }

          // compute the right bin of the pulse shape using time calibration constants
          EcalTimeCalibConstantMap::const_iterator it = itime_->find(detid);
          EcalTimeCalibConstant itimeconst = 0;
          if (it != itime_->end()) {
            itimeconst = (*it);
          } else {
            edm::LogError("EcalUncalibRecHit") << "No time intercalib const found for xtal " << detid.rawId()
                                               << "! something wrong with EcalTimeCalibConstants in your DB? ";
          }
          // intelligence for recHit computation
          const float offsetTime = offtime_->getEBValue();

          float correctedTime = (crh.timeMax - 5) * clockToNsConstant + itimeconst + offsetTime;
          float cterm = timeConstantTerm_;
          float sigmaped = pedRMSVec[0];  // approx for lower gains
          float nterm = timeNconst_ * sigmaped / uncalibRecHit.amplitude();
          float sigmat = std::sqrt(nterm * nterm + cterm * cterm);
          if ((correctedTime > sigmat * outOfTimeThreshP) || (correctedTime < -sigmat * outOfTimeThreshM)) {
            uncalibRecHit.setFlagBit(EcalUncalibratedRecHit::kOutOfTime);
          }
        }
      } else if (timealgo_ == weightsMethod) {
        //  weights method on the PU subtracted pulse shape
        std::vector<double> amplitudes;
        for (unsigned int ibx = 0; ibx < activeBX_.size(); ++ibx)
          amplitudes.push_back(uncalibRecHit.outOfTimeAmplitude(ibx));

        const auto& gid = grps_->barrel(hashedIndex);
        EcalPh2TBWeights::EcalTDCId tdcid(1);
        EcalPh2TBWeights::EcalTBWeightMap const& wgtsMap = wgts_->getMap();
        EcalPh2TBWeights::EcalTBWeightMap::const_iterator wit;
        wit = wgtsMap.find(std::make_pair(gid, tdcid));
        if (wit == wgtsMap.end()) {
          edm::LogWarning("EcalUncalibRecHit")
              << "No weights found for EcalGroupId: " << gid.id() << " and  EcalTDCId: " << tdcid
              << "\n  skipping digi with id: " << detid.rawId();
          result.pop_back();
          continue;
        }
        const auto& wset = wit->second;  // this is the EcalWeightSet

        const auto& mat1 = wset.getWeightsBeforeGainSwitch();
        const auto& mat2 = wset.getWeightsAfterGainSwitch();

        weights_[0] = &mat1;
        weights_[1] = &mat2;

        const double timerh = weightsMethodAlgo_->time(df, amplitudes, aped, aGain, fullpulse, weights_);
        uncalibRecHit.setJitter(timerh);
        uncalibRecHit.setJitterError(0.);  // not computed with weights

      } else if (timealgo_ == crossCorrelationMethod) {
        std::vector<double> amplitudes(activeBX_.size());
        for (unsigned int ibx = 0; ibx < activeBX_.size(); ++ibx)
          amplitudes[ibx] = uncalibRecHit.outOfTimeAmplitude(ibx);

        float jitterError = 0.;
        const float jitter = 0.;  // cross correlation timing for Phase 2 not yet implemented
        const float noCorrectedJitter = 0.;
        //float jitter =
        //    computeCC_->computeTimeCC(
        //        *itdg, amplitudes, aped, aGain, fullpulse, uncalibRecHit, jitterError, CCtargetTimePrecision_, true) +
        //    CCTimeShiftWrtRations_ / ecalcctiming::clockToNS;
        //float noCorrectedJitter = computeCC_->computeTimeCC(*itdg,
        //                                                    amplitudes,
        //                                                    aped,
        //                                                    aGain,
        //                                                    fullpulse,
        //                                                    uncalibRecHit,
        //                                                    jitterError,
        //                                                    CCtargetTimePrecisionForDelayedPulses_,
        //                                                    false) +
        //                          CCTimeShiftWrtRations_ / ecalcctiming::clockToNS;

        uncalibRecHit.setJitter(jitter);
        uncalibRecHit.setNonCorrectedTime(jitter, noCorrectedJitter);
        uncalibRecHit.setJitterError(jitterError);

      } else {  // no time method;
        uncalibRecHit.setJitter(0.);
        uncalibRecHit.setJitterError(0.);
      }
    }
  }
}

edm::ParameterSetDescription EcalUncalibRecHitWorkerMultiFitPh2::getAlgoDescription() {
  edm::ParameterSetDescription psd;

  // parameters for the different time algos
  edm::ParameterSwitch<std::string> timeAlgoNode(
      edm::ParameterDescription<std::string>("timealgo", "RatioMethod", true),
      "RatioMethod" >>
              (edm::ParameterDescription<std::vector<double>>("timeFitParameters",
                                                              {-2.015452e+00,
                                                               3.130702e+00,
                                                               -1.234730e+01,
                                                               4.188921e+01,
                                                               -8.283944e+01,
                                                               9.101147e+01,
                                                               -5.035761e+01,
                                                               1.105621e+01},
                                                              true) and
               edm::ParameterDescription<std::vector<double>>("amplitudeFitParameters", {1.138, 1.652}, true) and
               edm::ParameterDescription<double>("timeFitLimits_Lower", 0.2, true) and
               edm::ParameterDescription<double>("timeFitLimits_Upper", 1.4, true) and
               edm::ParameterDescription<double>("timeConstantTerm", .6, true) and
               edm::ParameterDescription<double>("timeNconst", 28.5, true) and
               edm::ParameterDescription<double>("outOfTimeThresholdGain10p", 5, true) and
               edm::ParameterDescription<double>("outOfTimeThresholdGain10m", 5, true) and
               edm::ParameterDescription<double>("outOfTimeThresholdGain1p", 5, true) and
               edm::ParameterDescription<double>("outOfTimeThresholdGain1m", 5, true) and
               edm::ParameterDescription<double>("amplitudeThreshold", 10, true)) or
          "WeightsMethod" >> edm::EmptyGroupDescription() or
          "crossCorrelationMethod" >>
              (edm::ParameterDescription<double>("crossCorrelationStartTime", -15.0, true) and
               edm::ParameterDescription<double>("crossCorrelationStopTime", 25.0, true) and
               edm::ParameterDescription<double>("crossCorrelationTargetTimePrecision", 0.01, true) and
               edm::ParameterDescription<double>("crossCorrelationTargetTimePrecisionForDelayedPulses", 0.05, true) and
               edm::ParameterDescription<double>("crossCorrelationTimeShiftWrtRations", 1., true) and
               edm::ParameterDescription<double>("crossCorrelationMinTimeToBeLateMin", 2., true) and
               edm::ParameterDescription<double>("crossCorrelationMinTimeToBeLateMax", 5., true)));

  psd.addNode(edm::ParameterDescription<std::vector<int>>("activeBXs", {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4}, true) and
              edm::ParameterDescription<bool>("ampErrorCalculation", true, true) and
              edm::ParameterDescription<bool>("useLumiInfoRunHeader", true, true) and
              edm::ParameterDescription<int>("bunchSpacing", 0, true) and
              edm::ParameterDescription<bool>("doPrefit", false, true) and
              edm::ParameterDescription<double>("prefitMaxChiSq", 25., true) and
              edm::ParameterDescription<bool>("dynamicPedestals", false, true) and
              edm::ParameterDescription<bool>("mitigateBadSamples", false, true) and
              edm::ParameterDescription<bool>("gainSwitchUseMaxSample", false, true) and
              edm::ParameterDescription<bool>("selectiveBadSampleCriteria", false, true) and
              edm::ParameterDescription<double>("addPedestalUncertainty", 0., true) and
              edm::ParameterDescription<bool>("simplifiedNoiseModelForGainSwitch", true, true) and timeAlgoNode);

  return psd;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitWorkerFactory,
                  EcalUncalibRecHitWorkerMultiFitPh2,
                  "EcalUncalibRecHitWorkerMultiFitPh2");
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitFillDescriptionWorkerFactory.h"
DEFINE_EDM_PLUGIN(EcalUncalibRecHitFillDescriptionWorkerFactory,
                  EcalUncalibRecHitWorkerMultiFitPh2,
                  "EcalUncalibRecHitWorkerMultiFitPh2");
