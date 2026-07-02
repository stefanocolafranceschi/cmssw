#ifndef RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationCommonKernelsPhase2_h
#define RecoLocalCalo_EcalRecProducers_plugins_alpaka_AmplitudeComputationCommonKernelsPhase2_h

#include <cstdlib>
#include <limits>
#include <alpaka/alpaka.hpp>

#include "CondFormats/EcalObjects/interface/alpaka/EcalMultifitConditionsPhase2Device.h"
#include "DataFormats/EcalDigi/interface/alpaka/EcalDigiPhase2DeviceCollection.h"
#include "DataFormats/EcalRecHit/interface/alpaka/EcalUncalibratedRecHitDeviceCollection.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"
#include "DataFormats/EcalDigi/interface/EcalLiteDTUSample.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EigenMatrixTypes_gpu.h"

#include "DeclsForKernels.h"
#include "KernelHelpers.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit {

  // Phase-2 preparation kernels.
  //
  // CATIA 2-gain LiteDTU decode; no gain-switch noise machinery, mirroring the
  // CPU EcalUncalibRecHitMultiFitAlgoPh2 where hasGainSwitch is always false.
  // Saturation handling mirrors EcalUncalibRecHitWorkerMultiFitPh2:
  //  - saturation (gain 1 and ADC at max) on the expected max sample
  //    -> amplitude = MAXADC * gain10, flagged kSaturated, no fit;
  //  - saturation on any other sample -> max-sample amplitude, kSaturated, no fit.

  ALPAKA_FN_ACC ALPAKA_FN_INLINE ::ecal::multifit::Ph2::SampleVector::Scalar decodeSample(
      EcalMultifitConditionsPhase2Device::ConstView const& conditionsDev,
      uint16_t const digiSample,
      uint32_t const hashedId,
      ::ecal::multifit::Ph2::SampleVector::Scalar& pedestal) {
    using DataType = ::ecal::multifit::Ph2::SampleVector::Scalar;
    auto const adc = ecalLiteDTU::adc(digiSample);
    auto const gainId = ecalLiteDTU::gainId(digiSample);

    DataType gainratio;
    if (gainId == static_cast<int>(ecalPh2::gainId1)) {
      pedestal = conditionsDev.pedestals_mean_g1()[hashedId];
      gainratio = conditionsDev.gain10Over1()[hashedId];
    } else {
      pedestal = conditionsDev.pedestals_mean_g10()[hashedId];
      gainratio = 1.;
    }

    // amplitudes on the gain-10 scale
    // NOTE: the CPU EcalUncalibRecHitMultiFitAlgoPh2 multiplies ALL samples by the
    // gain ratio (also gain-10 ones); here the ratio is applied only to gain-1
    // samples, consistently with the worker's saturation path and with the Phase-2
    // weights reconstruction. To be confirmed with the reco leads.
    return (static_cast<DataType>(adc) - pedestal) * gainratio;
  }

  ///
  /// assume kernel launch configuration is
  /// (MAXSAMPLES * nchannels, blocks)
  ///
  class Kernel_prep_1d_and_initialize {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  EcalDigiPhase2DeviceCollection::ConstView digisDevEB,
                                  EcalUncalibratedRecHitDeviceCollection::View uncalibRecHitsEB,
                                  EcalMultifitConditionsPhase2Device::ConstView conditionsDev,
                                  ::ecal::multifit::Ph2::SampleVector* amplitudesForFit,
                                  char* acState) const {
      using DataType = ::ecal::multifit::Ph2::SampleVector::Scalar;
      constexpr auto nsamples = static_cast<int>(EcalDataFrame_Ph2::MAXSAMPLES);
      constexpr int sampleMax = static_cast<int>(ecalPh2::kMaxSampleIdx);

      auto const nchannels = digisDevEB.size();
      auto const totalElements = nchannels * nsamples;

      auto const* digis_in = digisDevEB.data().data()->data();
      auto const dids = digisDevEB.id();

      // shared memory: per-sample saturation flags
      bool* shr_isSaturatedSample = alpaka::getDynSharedMem<bool>(acc);

      for (auto block : cms::alpakatools::uniform_groups(acc, totalElements)) {
        // phase A: per-sample decode and initialization
        for (auto idx : cms::alpakatools::uniform_group_elements(acc, block, totalElements)) {
          // set the output collection size scalar
          if (idx.global == 0) {
            uncalibRecHitsEB.size() = nchannels;
          }

          auto const ch = idx.global / nsamples;
          auto const sample = static_cast<int>(idx.global % nsamples);
          int const inputCh = ch;
          int const inputTx = idx.global;

          auto const did = DetId{dids[inputCh]};
          auto const hashedId = reconstruction::hashedIndexEB(did.rawId());

          DataType pedestal;
          auto const amplitude = decodeSample(conditionsDev, digis_in[inputTx], hashedId, pedestal);
          amplitudesForFit[ch](sample) = amplitude;

          // per-sample saturation criterion (CPU worker: gain 1 and ADC at max range)
          auto const gainId = ecalLiteDTU::gainId(digis_in[inputTx]);
          auto const adc = ecalLiteDTU::adc(digis_in[inputTx]);
          shr_isSaturatedSample[idx.local] = (gainId == static_cast<int>(ecalPh2::gainId1)) &&
                                             (adc == static_cast<int>(ecalPh2::MAXADC));

          // zero the stored out-of-time amplitudes
          // FIXME (Phase 2): the data-format array holds ecalPh1::sampleSize entries
          if (sample < static_cast<int>(ecalPh1::sampleSize))
            uncalibRecHitsEB.outOfTimeAmplitudes().data()[inputCh][sample] = 0.;
        }

        alpaka::syncBlockThreads(acc);

        // phase B: one thread per channel finalizes the initialization
        for (auto idx : cms::alpakatools::uniform_group_elements(acc, block, totalElements)) {
          auto const ch = idx.global / nsamples;
          auto const sample = static_cast<int>(idx.global % nsamples);

          if (sample != sampleMax)
            continue;

          int const inputCh = ch;
          int const inputTx = idx.global;
          auto const chStart = idx.local - sampleMax;

          auto const did = DetId{dids[inputCh]};
          auto const hashedId = reconstruction::hashedIndexEB(did.rawId());

          // re-decode this thread's (max) sample for pedestal and amplitude
          DataType pedestal;
          auto const maxSampleAmplitude = decodeSample(conditionsDev, digis_in[inputTx], hashedId, pedestal);

          // mirror the CPU worker's lastSampleBeforeSaturation scan
          int lastSampleBeforeSaturation = -2;
          for (int iSample = 0; iSample < nsamples; ++iSample) {
            if (shr_isSaturatedSample[chStart + iSample]) {
              lastSampleBeforeSaturation = iSample - 1;
              break;
            }
          }

          auto energies = uncalibRecHitsEB.amplitude();
          auto amplitudeErrors = uncalibRecHitsEB.amplitudeError();
          auto chi2s = uncalibRecHitsEB.chi2();
          auto ootChi2s = uncalibRecHitsEB.OOTchi2();
          auto g_pedestal = uncalibRecHitsEB.pedestal();
          auto jitters = uncalibRecHitsEB.jitter();
          auto jitterErrors = uncalibRecHitsEB.jitterError();
          auto dids_out = uncalibRecHitsEB.id();
          auto flags = uncalibRecHitsEB.flags();
          auto auxs = uncalibRecHitsEB.aux();

          dids_out[inputCh] = did.rawId();
          g_pedestal[inputCh] = pedestal;  // pedestal of the max sample (CPU pedval)
          amplitudeErrors[inputCh] = 0.;
          chi2s[inputCh] = 0.;
          ootChi2s[inputCh] = 0.;
          jitters[inputCh] = 0.;  // timing computation deferred for Phase 2
          jitterErrors[inputCh] = 0.;
          auxs[inputCh] = 0;

          uint32_t flag = 0;
          if (lastSampleBeforeSaturation == sampleMax - 1) {
            // saturation on the expected max sample: amplitude at the ADC maximum
            energies[inputCh] = static_cast<float>(ecalPh2::MAXADC) * ecalPh2::gains[ecalPh2::gainId10];
            flag |= 0x1 << EcalUncalibratedRecHit::kSaturated;
            acState[ch] = static_cast<char>(MinimizationState::Precomputed);
          } else if (lastSampleBeforeSaturation >= -1) {
            // saturation elsewhere: use the max-sample amplitude, no extrapolation
            energies[inputCh] = maxSampleAmplitude;
            flag |= 0x1 << EcalUncalibratedRecHit::kSaturated;
            acState[ch] = static_cast<char>(MinimizationState::Precomputed);
          } else {
            energies[inputCh] = 0.;
            acState[ch] = static_cast<char>(MinimizationState::NotFinished);
          }
          flags[inputCh] = flag;
        }
      }
    }
  };

  ///
  /// assume kernel launch configuration is
  /// ([MAXSAMPLES, MAXSAMPLES], nchannels)
  ///
  class Kernel_prep_2d {
  public:
    ALPAKA_FN_ACC void operator()(Acc2D const& acc,
                                  EcalDigiPhase2DeviceCollection::ConstView digisDevEB,
                                  EcalMultifitConditionsPhase2Device::ConstView conditionsDev,
                                  ::ecal::multifit::Ph2::SampleMatrix* noisecov,
                                  SamplePulseMatrixPhase2* pulse_matrix) const {
      constexpr auto nsamples = static_cast<int>(EcalDataFrame_Ph2::MAXSAMPLES);
      constexpr int tsamples = static_cast<int>(EcalPh2PulseShape::TEMPLATESAMPLES);

      auto const blockDimX = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[1u];
      auto const elemsPerBlockX = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[1u];
      auto const elemsPerBlockY = alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u];
      Vec2D const size_2d = {elemsPerBlockY, blockDimX * elemsPerBlockX};  // {y, x} coordinates

      for (auto ndindex : cms::alpakatools::uniform_elements_nd(acc, size_2d)) {
        auto const ch = ndindex[1] / nsamples;
        auto const tx = static_cast<int>(ndindex[1] % nsamples);
        auto const ty = static_cast<int>(ndindex[0]);

        auto const dids = digisDevEB.id();
        auto const did = DetId{dids[ch]};
        auto const hashedId = ecal::reconstruction::hashedIndexEB(did.rawId());

        // noise covariance: no gain switch handling (CPU algo has hasGainSwitch
        // always false) -> gain-10 pedestal rms and gain-10 sample correlation
        auto const g10SamplesCorrelation = conditionsDev.sampleCorrelation_g10().data();
        auto const vidx = std::abs(ty - tx);
        auto const rms_g10 = conditionsDev.pedestals_rms_g10()[hashedId];
        noisecov[ch](ty, tx) = rms_g10 * rms_g10 * g10SamplesCorrelation[vidx];

        // pulse matrix: nsamples x kNPulsesPhase2 (not square as in Phase 1);
        // mirrors PulseChiSqSNNLS<ecalPh2>: template index = sample - shift with
        // shift = kTemplateStartSamplePhase2 + kBxToSampleShiftPhase2 * bx
        if (tx < kNPulsesPhase2) {
          int const bx = tx - kInTimePulseIdxPhase2;
          int const tmpl = ty - kTemplateStartSamplePhase2 - kBxToSampleShiftPhase2 * bx;
          float const value = (tmpl >= 0 && tmpl < tsamples) ? conditionsDev.pulseShapes()[hashedId][tmpl] : 0.f;
          pulse_matrix[ch](ty, tx) = value;
        }
      }
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit

namespace alpaka::trait {
  using namespace ALPAKA_ACCELERATOR_NAMESPACE;
  using namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit;

  //! The trait for getting the size of the block shared dynamic memory for Kernel_prep_1d_and_initialize.
  template <>
  struct BlockSharedMemDynSizeBytes<Kernel_prep_1d_and_initialize, Acc1D> {
    //! \return The size of the shared memory allocated for a block.
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(Kernel_prep_1d_and_initialize const&,
                                                                 TVec const& threadsPerBlock,
                                                                 TVec const& elemsPerThread,
                                                                 TArgs const&...) -> std::size_t {
      // per-sample saturation flags
      std::size_t bytes = threadsPerBlock[0u] * elemsPerThread[0u] * sizeof(bool);
      return bytes;
    }
  };
}  // namespace alpaka::trait

#endif  // RecoLocalCalo_EcalRecProducers_plugins_AmplitudeComputationCommonKernelsPhase2_h
