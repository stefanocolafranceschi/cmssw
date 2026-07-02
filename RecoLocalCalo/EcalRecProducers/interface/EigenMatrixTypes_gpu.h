#ifndef RecoLocalCalo_EcalRecProducers_EigenMatrixTypes_gpu_h
#define RecoLocalCalo_EcalRecProducers_EigenMatrixTypes_gpu_h

#include <array>
#include <Eigen/Dense>

#include "DataFormats/EcalRecHit/interface/RecoTypes.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"

namespace ecal {
  namespace multifit {

    using data_type = ::ecal::reco::ComputationScalarType;

    // Eigen matrix types for the multifit, templated on the ECAL phase constants
    // class (ecalPh1 / ecalPh2) so that Phase 1 (10 samples) and Phase 2
    // (16 samples) code can share the same definitions.
    //
    // TNGains is the number of gains used for the sample correlation / noise
    // covariance handling:
    //   Phase 1: ecalPh1::NGAINS - 1 (the MGPA zero gain is not counted)
    //   Phase 2: ecalPh2::NGAINS (both CATIA gains are counted)
    template <typename TEcalConsts, int TNGains>
    struct EigenMatrixTypesT {
      static constexpr int SampleVectorSize = static_cast<int>(TEcalConsts::sampleSize);
      static constexpr int FullSampleVectorSize = static_cast<int>(TEcalConsts::kFullSampleVectorSize);
      static constexpr int PulseVectorSize = static_cast<int>(TEcalConsts::kPulseShapeTemplateSampleSize);
      static constexpr int NGains = TNGains;

      using PulseMatrixType = Eigen::Matrix<data_type, SampleVectorSize, SampleVectorSize>;
      using BXVectorType = Eigen::Matrix<char, SampleVectorSize, 1>;
      using SampleMatrixD = Eigen::Matrix<double, SampleVectorSize, SampleVectorSize>;

      using SampleVector = Eigen::Matrix<data_type, SampleVectorSize, 1>;
      using FullSampleVector = Eigen::Matrix<data_type, FullSampleVectorSize, 1>;
      using PulseVector = Eigen::Matrix<data_type, Eigen::Dynamic, 1, 0, PulseVectorSize, 1>;
      using BXVector = Eigen::Matrix<char, Eigen::Dynamic, 1, 0, PulseVectorSize, 1>;
      using SampleGainVector = Eigen::Matrix<char, SampleVectorSize, 1>;
      using SampleMatrix = Eigen::Matrix<data_type, SampleVectorSize, SampleVectorSize>;
      using FullSampleMatrix = Eigen::Matrix<data_type, FullSampleVectorSize, FullSampleVectorSize>;
      using PulseMatrix = Eigen::Matrix<data_type, Eigen::Dynamic, Eigen::Dynamic, 0, PulseVectorSize, PulseVectorSize>;
      using SamplePulseMatrix =
          Eigen::Matrix<data_type, SampleVectorSize, Eigen::Dynamic, 0, SampleVectorSize, PulseVectorSize>;
      using SampleDecompLLT = Eigen::LLT<SampleMatrix>;
      using SampleDecompLLTD = Eigen::LLT<SampleMatrixD>;
      using PulseDecompLLT = Eigen::LLT<PulseMatrix>;
      using PulseDecompLDLT = Eigen::LDLT<PulseMatrix>;

      using SingleMatrix = Eigen::Matrix<data_type, 1, 1>;
      using SingleVector = Eigen::Matrix<data_type, 1, 1>;

      using SampleMatrixGainArray = std::array<SampleMatrixD, NGains>;

      using PermutationMatrix = Eigen::PermutationMatrix<SampleVectorSize>;
    };

    using Ph1Types = EigenMatrixTypesT<ecalPh1, static_cast<int>(ecalPh1::NGAINS) - 1>;  // do not count 0 gain
    using Ph2Types = EigenMatrixTypesT<ecalPh2, static_cast<int>(ecalPh2::NGAINS)>;

    // Phase 1 aliases at ecal::multifit namespace scope, kept for backward
    // compatibility with the existing Run 3 reconstruction code.
    constexpr int SampleVectorSize = Ph1Types::SampleVectorSize;
    constexpr int FullSampleVectorSize = Ph1Types::FullSampleVectorSize;
    constexpr int PulseVectorSize = Ph1Types::PulseVectorSize;
    constexpr int NGains = Ph1Types::NGains;

    using PulseMatrixType = Ph1Types::PulseMatrixType;
    using BXVectorType = Ph1Types::BXVectorType;
    using SampleMatrixD = Ph1Types::SampleMatrixD;

    using SampleVector = Ph1Types::SampleVector;
    using FullSampleVector = Ph1Types::FullSampleVector;
    using PulseVector = Ph1Types::PulseVector;
    using BXVector = Ph1Types::BXVector;
    using SampleGainVector = Ph1Types::SampleGainVector;
    using SampleMatrix = Ph1Types::SampleMatrix;
    using FullSampleMatrix = Ph1Types::FullSampleMatrix;
    using PulseMatrix = Ph1Types::PulseMatrix;
    using SamplePulseMatrix = Ph1Types::SamplePulseMatrix;
    using SampleDecompLLT = Ph1Types::SampleDecompLLT;
    using SampleDecompLLTD = Ph1Types::SampleDecompLLTD;
    using PulseDecompLLT = Ph1Types::PulseDecompLLT;
    using PulseDecompLDLT = Ph1Types::PulseDecompLDLT;

    using SingleMatrix = Ph1Types::SingleMatrix;
    using SingleVector = Ph1Types::SingleVector;

    using SampleMatrixGainArray = Ph1Types::SampleMatrixGainArray;

    using PermutationMatrix = Ph1Types::PermutationMatrix;

    // Phase 2 types (16 samples, 2 CATIA gains) for the Phase 2 kernels.
    namespace Ph2 {
      constexpr int SampleVectorSize = Ph2Types::SampleVectorSize;
      constexpr int FullSampleVectorSize = Ph2Types::FullSampleVectorSize;
      constexpr int PulseVectorSize = Ph2Types::PulseVectorSize;
      constexpr int NGains = Ph2Types::NGains;

      using PulseMatrixType = Ph2Types::PulseMatrixType;
      using BXVectorType = Ph2Types::BXVectorType;
      using SampleMatrixD = Ph2Types::SampleMatrixD;

      using SampleVector = Ph2Types::SampleVector;
      using FullSampleVector = Ph2Types::FullSampleVector;
      using PulseVector = Ph2Types::PulseVector;
      using BXVector = Ph2Types::BXVector;
      using SampleGainVector = Ph2Types::SampleGainVector;
      using SampleMatrix = Ph2Types::SampleMatrix;
      using FullSampleMatrix = Ph2Types::FullSampleMatrix;
      using PulseMatrix = Ph2Types::PulseMatrix;
      using SamplePulseMatrix = Ph2Types::SamplePulseMatrix;
      using SampleDecompLLT = Ph2Types::SampleDecompLLT;
      using SampleDecompLLTD = Ph2Types::SampleDecompLLTD;
      using PulseDecompLLT = Ph2Types::PulseDecompLLT;
      using PulseDecompLDLT = Ph2Types::PulseDecompLDLT;

      using SingleMatrix = Ph2Types::SingleMatrix;
      using SingleVector = Ph2Types::SingleVector;

      using SampleMatrixGainArray = Ph2Types::SampleMatrixGainArray;

      using PermutationMatrix = Ph2Types::PermutationMatrix;
    }  // namespace Ph2

  }  // namespace multifit
}  // namespace ecal

#endif  // RecoLocalCalo_EcalRecProducers_EigenMatrixTypes_gpu_h
