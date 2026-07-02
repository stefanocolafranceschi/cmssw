#include <cmath>
#include <limits>
#include <alpaka/alpaka.hpp>

#include "CondFormats/EcalObjects/interface/EcalPulseCovarianceT.h"
#include "DataFormats/CaloRecHit/interface/MultifitComputations.h"
#include "FWCore/Utilities/interface/CMSUnrollLoop.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "AmplitudeComputationKernelsPhase2.h"
#include "KernelHelpers.h"
#include "EcalUncalibRecHitMultiFitAlgoPhase2Portable.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit {

  using namespace ::ecal::multifit::Ph2;

  // Add the pulse-shape covariance contribution of every active pulse to the
  // noise covariance. Mirrors PulseChiSqSNNLS<ecalPh2>::updateCov, with the
  // 32x32 zero-padded full-pulse covariance replaced by direct, bounds-guarded
  // access to the 16x16 EcalPh2PulseCovariance template:
  //   template index = sample - shift, shift = kTemplateStartSamplePhase2 +
  //   kBxToSampleShiftPhase2 * bx, bx = ipulse - kInTimePulseIdxPhase2.
  template <typename MatrixType, typename AmplitudeVectorType>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void update_covariance(EcalPh2PulseCovariance const& pulse_covariance,
                                                        MatrixType& inverse_cov,
                                                        AmplitudeVectorType const& amplitudes) {
    constexpr int nsamples = SampleVector::RowsAtCompileTime;
    constexpr int npulses = kNPulsesPhase2;
    constexpr int tsamples = EcalPh2PulseCovariance::TEMPLATESAMPLES;

    CMS_UNROLL_LOOP
    for (int ipulse = 0; ipulse < npulses; ++ipulse) {
      auto const amplitude = amplitudes.coeff(ipulse);
      if (amplitude == 0)
        continue;

      int const bx = ipulse - kInTimePulseIdxPhase2;
      int const shift = kTemplateStartSamplePhase2 + kBxToSampleShiftPhase2 * bx;
      int const first_sample_t = std::max(0, shift);

      auto const value_sq = amplitude * amplitude;

      for (int col = first_sample_t; col < nsamples; ++col) {
        int const tc = col - shift;
        if (tc < 0 || tc >= tsamples)
          continue;
        for (int row = col; row < nsamples; ++row) {
          int const tr = row - shift;
          if (tr >= tsamples)
            break;
          inverse_cov(row, col) += value_sq * pulse_covariance.covval[tr][tc];
        }
      }
    }
  }

  ///
  /// launch ctx parameters are (nchannels / block, blocks)
  ///
  /// Phase-2 minimization: 16 samples, kNPulsesPhase2 (=5) pulses with
  /// activeBXs = {-2,-1,0,1,2}; in-time amplitude at index kInTimePulseIdxPhase2.
  ///
  /// Conventions:
  ///   - amplitudes -> solution vector, what we are fitting for
  ///   - samples -> raw detector responses
  ///   - passive constraint - satisfied constraint
  ///   - active constraint - unsatisfied (yet) constraint
  ///
  class Kernel_minimize {
  public:
    ALPAKA_FN_ACC void operator()(Acc1D const& acc,
                                  InputProduct::ConstView const& digisDevEB,
                                  OutputProduct::View uncalibRecHitsEB,
                                  EcalMultifitConditionsPhase2Device::ConstView conditionsDev,
                                  ::ecal::multifit::Ph2::SampleMatrix const* noisecov,
                                  SamplePulseMatrixPhase2 const* pulse_matrix,
                                  ::ecal::multifit::Ph2::SampleVector const* samples,
                                  char* acState,
                                  int max_iterations) const {
      constexpr int NSAMPLES = SampleMatrix::RowsAtCompileTime;
      constexpr int NPULSES = kNPulsesPhase2;
      static_assert(NPULSES <= NSAMPLES);

      using DataType = SampleVector::Scalar;

      // shared memory layout per element: one NSAMPLES symmetric matrix (noise
      // covariance, later reused as the fnnls decomposition) + one NPULSES
      // symmetric matrix (AtA)
      constexpr auto covTotal = calo::multifit::MapSymM<DataType, NSAMPLES>::total;
      constexpr auto pulTotal = calo::multifit::MapSymM<DataType, NPULSES>::total;

      auto const elemsPerBlock(alpaka::getWorkDiv<alpaka::Block, alpaka::Elems>(acc)[0u]);

      auto const nchannels = digisDevEB.size();

      auto const* pulse_covariance =
          reinterpret_cast<const EcalPh2PulseCovariance*>(conditionsDev.pulseCovariance().data());

      // shared memory
      DataType* shrmem = alpaka::getDynSharedMem<DataType>(acc);

      // channel
      for (auto idx : cms::alpakatools::uniform_elements(acc, nchannels)) {
        if (static_cast<MinimizationState>(acState[idx]) == MinimizationState::Precomputed)
          continue;

        auto const elemIdx = idx % elemsPerBlock;

        // shared memory pointers
        DataType* shrCovStorage = shrmem + (covTotal + pulTotal) * elemIdx;
        DataType* shrAtAStorage = shrCovStorage + covTotal;
        // the fnnls L matrix reuses the covariance storage (no longer needed there)
        DataType* shrMatrixLForFnnlsStorage = shrCovStorage;

        auto* amplitudes = uncalibRecHitsEB.outOfTimeAmplitudes().data();
        auto energies = uncalibRecHitsEB.amplitude();
        auto chi2s = uncalibRecHitsEB.chi2();

        // get the hash
        auto const dids = digisDevEB.id();
        auto const did = DetId{dids[idx]};
        auto const hashedId = ecal::reconstruction::hashedIndexEB(did.rawId());

        // inits
        int npassive = 0;

        calo::multifit::ColumnVector<NPULSES, int> pulseOffsets;
        CMS_UNROLL_LOOP
        for (int i = 0; i < NPULSES; ++i)
          pulseOffsets(i) = i;

        calo::multifit::ColumnVector<NPULSES, DataType> resultAmplitudes;
        CMS_UNROLL_LOOP
        for (int counter = 0; counter < NPULSES; ++counter)
          resultAmplitudes(counter) = 0;

        float chi2 = 0, chi2_now = 0;

        // loop for up to max_iterations
        for (int iter = 0; iter < max_iterations; ++iter) {
          DataType* covMatrixStorage = shrCovStorage;
          calo::multifit::MapSymM<DataType, NSAMPLES> covMatrix{covMatrixStorage};
          int counter = 0;
          CMS_UNROLL_LOOP
          for (int col = 0; col < NSAMPLES; ++col) {
            CMS_UNROLL_LOOP
            for (int row = col; row < NSAMPLES; ++row) {
              covMatrixStorage[counter++] = noisecov[idx].coeffRef(row, col);
            }
          }
          update_covariance(pulse_covariance[hashedId], covMatrix, resultAmplitudes);

          // compute actual covariance decomposition
          DataType matrixLStorage[covTotal];
          calo::multifit::MapSymM<DataType, NSAMPLES> matrixL{matrixLStorage};
          calo::multifit::compute_decomposition_unrolled(matrixL, covMatrix);

          // L * A = P
          calo::multifit::ColMajorMatrix<NSAMPLES, NPULSES> A;
          calo::multifit::solve_forward_subst_matrix(A, pulse_matrix[idx], matrixL);

          // L b = s
          float reg_b[NSAMPLES];
          calo::multifit::solve_forward_subst_vector(reg_b, samples[idx], matrixL);

          calo::multifit::MapSymM<DataType, NPULSES> AtA{shrAtAStorage};
          calo::multifit::ColumnVector<NPULSES, DataType> Atb;
          CMS_UNROLL_LOOP
          for (int icol = 0; icol < NPULSES; ++icol) {
            float reg_ai[NSAMPLES];

            // load column icol
            CMS_UNROLL_LOOP
            for (int counter = 0; counter < NSAMPLES; ++counter)
              reg_ai[counter] = A(counter, icol);

            // compute diagonal
            float sum = 0.f;
            CMS_UNROLL_LOOP
            for (int counter = 0; counter < NSAMPLES; ++counter)
              sum += reg_ai[counter] * reg_ai[counter];

            // store
            AtA(icol, icol) = sum;

            // go thru the other columns
            CMS_UNROLL_LOOP
            for (int j = icol + 1; j < NPULSES; ++j) {
              // load column j
              float reg_aj[NSAMPLES];
              CMS_UNROLL_LOOP
              for (int counter = 0; counter < NSAMPLES; ++counter)
                reg_aj[counter] = A(counter, j);

              // accum
              float sum = 0.f;
              CMS_UNROLL_LOOP
              for (int counter = 0; counter < NSAMPLES; ++counter)
                sum += reg_aj[counter] * reg_ai[counter];

              // store
              AtA(j, icol) = sum;
            }

            // Atb accum
            float sum_atb = 0.f;
            CMS_UNROLL_LOOP
            for (int counter = 0; counter < NSAMPLES; ++counter)
              sum_atb += reg_ai[counter] * reg_b[counter];

            // store atb
            Atb(icol) = sum_atb;
          }

          calo::multifit::MapSymM<DataType, NPULSES> matrixLForFnnls{shrMatrixLForFnnlsStorage};

          calo::multifit::fnnls(AtA,
                                Atb,
                                resultAmplitudes,
                                npassive,
                                pulseOffsets,
                                matrixLForFnnls,
                                1e-11,
                                500,
                                16,
                                2);

          calo::multifit::calculateChiSq(matrixL, pulse_matrix[idx], resultAmplitudes, samples[idx], chi2_now);

          auto const deltachi2 = chi2_now - chi2;
          chi2 = chi2_now;

          if (std::abs(deltachi2) < 1e-3)
            break;
        }

        // store to global output values
        chi2s[idx] = chi2;
        energies[idx] = resultAmplitudes(kInTimePulseIdxPhase2);

        // out-of-time amplitudes are stored at slot bx + 5 as in the CPU version
        // (slots 3..7 of the ecalPh1::sampleSize-entry data-format array)
        CMS_UNROLL_LOOP
        for (int counter = 0; counter < NPULSES; ++counter) {
          int const slot = (counter - kInTimePulseIdxPhase2) + 5;
          amplitudes[idx][slot] = resultAmplitudes(counter);
        }
      }
    }
  };

  void minimization_procedure(Queue& queue,
                              InputProduct const& digisDevEB,
                              OutputProduct& uncalibRecHitsDevEB,
                              EventDataForScratchDevicePhase2& scratch,
                              EcalMultifitConditionsPhase2Device const& conditionsDev,
                              ConfigurationParametersPhase2 const& configParams,
                              uint32_t const totalChannels) {
    // TODO: configure from python
    auto threads_min = configParams.kernelMinimizeThreads[0];
    auto blocks_min = cms::alpakatools::divide_up_by(totalChannels, threads_min);

    auto workDivMinimize = cms::alpakatools::make_workdiv<Acc1D>(blocks_min, threads_min);
    alpaka::exec<Acc1D>(queue,
                        workDivMinimize,
                        Kernel_minimize{},
                        digisDevEB.const_view(),
                        uncalibRecHitsDevEB.view(),
                        conditionsDev.const_view(),
                        reinterpret_cast<::ecal::multifit::Ph2::SampleMatrix*>(scratch.noisecovDevBuf.data()),
                        reinterpret_cast<SamplePulseMatrixPhase2*>(scratch.pulse_matrixDevBuf.data()),
                        reinterpret_cast<::ecal::multifit::Ph2::SampleVector*>(scratch.samplesDevBuf.data()),
                        scratch.acStateDevBuf.data(),
                        50);  // maximum number of fit iterations
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit

namespace alpaka::trait {
  using namespace ALPAKA_ACCELERATOR_NAMESPACE;
  using namespace ALPAKA_ACCELERATOR_NAMESPACE::ecal::multifit;

  //! The trait for getting the size of the block shared dynamic memory for Kernel_minimize.
  template <>
  struct BlockSharedMemDynSizeBytes<Kernel_minimize, Acc1D> {
    //! \return The size of the shared memory allocated for a block.
    template <typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(Kernel_minimize const&,
                                                                 TVec const& threadsPerBlock,
                                                                 TVec const& elemsPerThread,
                                                                 TArgs const&...) -> std::size_t {
      using ScalarType = ::ecal::multifit::Ph2::SampleVector::Scalar;

      // one NSAMPLES symmetric matrix + one NPULSES symmetric matrix per element
      constexpr auto covTotal =
          calo::multifit::MapSymM<ScalarType, ::ecal::multifit::Ph2::SampleVector::RowsAtCompileTime>::total;
      constexpr auto pulTotal = calo::multifit::MapSymM<ScalarType, kNPulsesPhase2>::total;
      std::size_t bytes = threadsPerBlock[0u] * elemsPerThread[0u] * (covTotal + pulTotal) * sizeof(ScalarType);
      return bytes;
    }
  };
}  // namespace alpaka::trait
