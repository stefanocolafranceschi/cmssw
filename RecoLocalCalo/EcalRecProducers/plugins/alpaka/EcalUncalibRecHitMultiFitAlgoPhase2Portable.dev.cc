#include <alpaka/alpaka.hpp>

#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "AmplitudeComputationCommonKernelsPhase2.h"
#include "AmplitudeComputationKernelsPhase2.h"
#include "EcalUncalibRecHitMultiFitAlgoPhase2Portable.h"
#include "DeclsForKernels.h"

//#define DEBUG
//#define ECAL_RECO_ALPAKA_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit {

  using namespace cms::alpakatools;

  void launchKernels(Queue& queue,
                     InputProduct const& digisDevEB,
                     OutputProduct& uncalibRecHitsDevEB,
                     EcalMultifitConditionsPhase2Device const& conditionsDev,
                     EcalMultifitParametersPhase2 const* /*paramsDev*/,  // needed again once timing is implemented
                     ConfigurationParametersPhase2 const& configParams) {
    auto constexpr kMaxSamples = EcalDataFrame_Ph2::MAXSAMPLES;

    auto const nChannels = static_cast<uint32_t>(uncalibRecHitsDevEB.const_view().metadata().size());

    EventDataForScratchDevicePhase2 scratch(configParams, nChannels, queue);

    //
    // 1d preparation kernel
    //
    uint32_t constexpr nchannels_per_block = 32;
    auto constexpr threads_1d = kMaxSamples * nchannels_per_block;
    auto const blocks_1d = cms::alpakatools::divide_up_by(nChannels * kMaxSamples, threads_1d);
    auto workDivPrep1D = cms::alpakatools::make_workdiv<Acc1D>(blocks_1d, threads_1d);
    // Since the ::ecal::multifit::Ph2::X objects are non-dynamic Eigen::Matrix types the returned
    // pointers from the buffers and the ::ecal::multifit::Ph2::X* both point to the data.
    alpaka::exec<Acc1D>(queue,
                        workDivPrep1D,
                        Kernel_prep_1d_and_initialize{},
                        digisDevEB.const_view(),
                        uncalibRecHitsDevEB.view(),
                        conditionsDev.const_view(),
                        reinterpret_cast<::ecal::multifit::Ph2::SampleVector*>(scratch.samplesDevBuf.data()),
                        scratch.acStateDevBuf.data());

    //
    // 2d preparation kernel
    //
    Vec2D const blocks_2d{1u, nChannels};  // {y, x} coordinates
    Vec2D const threads_2d{kMaxSamples, kMaxSamples};
    auto workDivPrep2D = cms::alpakatools::make_workdiv<Acc2D>(blocks_2d, threads_2d);
    alpaka::exec<Acc2D>(queue,
                        workDivPrep2D,
                        Kernel_prep_2d{},
                        digisDevEB.const_view(),
                        conditionsDev.const_view(),
                        reinterpret_cast<::ecal::multifit::Ph2::SampleMatrix*>(scratch.noisecovDevBuf.data()),
                        reinterpret_cast<SamplePulseMatrixPhase2*>(scratch.pulse_matrixDevBuf.data()));

    // run minimization kernels
    minimization_procedure(queue, digisDevEB, uncalibRecHitsDevEB, scratch, conditionsDev, configParams, nChannels);

    // The Phase-2 timing computation is deferred (the timing algorithm may still
    // change). The Phase-1-style timing kernels in TimeComputationKernelsPhase2.h
    // are not compatible with the Phase-2 conditions and are not launched.
    if (configParams.shouldRunTimingComputation) {
      throw cms::Exception("NotImplemented")
          << "The Phase-2 ECAL multifit timing computation on GPU is not implemented yet; "
          << "set shouldRunTimingComputation to False.";
    }
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit
