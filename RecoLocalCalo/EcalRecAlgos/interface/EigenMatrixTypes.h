#ifndef RecoLocalCalo_EcalRecAlgos_EigenMatrixTypes_h
#define RecoLocalCalo_EcalRecAlgos_EigenMatrixTypes_h

#include <Eigen/Dense>
#include <array>

template <class P>
struct EigenMatrixTypes {
  static constexpr int SampleVectorSize = static_cast<int>(P::sampleSize);
  static constexpr int FullSampleVectorSize = static_cast<int>(P::kFullSampleVectorSize);
  static constexpr int PulseVectorSize = static_cast<int>(P::kPulseShapeTemplateSampleSize);
  static constexpr int NGains =
      static_cast<int>((P::NGAINS > 2) ? P::NGAINS - 1 : P::NGAINS);  // Do not count gain 0 for Phase 1

  using SampleVector = Eigen::Matrix<double, SampleVectorSize, 1>;
  using FullSampleVector = Eigen::Matrix<double, FullSampleVectorSize, 1>;
  using PulseVector = Eigen::Matrix<double, Eigen::Dynamic, 1, 0, PulseVectorSize, 1>;
  using BXVector = Eigen::Matrix<char, Eigen::Dynamic, 1, 0, PulseVectorSize, 1>;
  using SampleGainVector = Eigen::Matrix<char, SampleVectorSize, 1>;
  using SampleMatrix = Eigen::Matrix<double, SampleVectorSize, SampleVectorSize>;
  using FullSampleMatrix = Eigen::Matrix<double, FullSampleVectorSize, FullSampleVectorSize>;
  using PulseMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, PulseVectorSize, PulseVectorSize>;
  using SamplePulseMatrix =
      Eigen::Matrix<double, SampleVectorSize, Eigen::Dynamic, 0, SampleVectorSize, PulseVectorSize>;
  using SampleDecompLLT = Eigen::LLT<SampleMatrix>;
  using PulseDecompLLT = Eigen::LLT<PulseMatrix>;
  using PulseDecompLDLT = Eigen::LDLT<PulseMatrix>;

  using SingleMatrix = Eigen::Matrix<double, 1, 1>;
  using SingleVector = Eigen::Matrix<double, 1, 1>;

  using SampleMatrixGainArray = std::array<SampleMatrix, NGains>;
};

#endif
