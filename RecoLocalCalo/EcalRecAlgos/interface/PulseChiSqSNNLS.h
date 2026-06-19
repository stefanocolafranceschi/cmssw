#ifndef RecoLocalCalo_EcalRecAlgos_PulseChiSqSNNLS_h
#define RecoLocalCalo_EcalRecAlgos_PulseChiSqSNNLS_h

#define EIGEN_NO_DEBUG  // kill throws in eigen code
#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include <set>
#include <array>

template <class P>
class PulseChiSqSNNLS {
public:
  using SampleVector = typename EigenMatrixTypes<P>::SampleVector;
  using FullSampleVector = typename EigenMatrixTypes<P>::FullSampleVector;
  using PulseVector = typename EigenMatrixTypes<P>::PulseVector;
  using BXVector = typename EigenMatrixTypes<P>::BXVector;
  using SampleGainVector = typename EigenMatrixTypes<P>::SampleGainVector;
  using SampleMatrix = typename EigenMatrixTypes<P>::SampleMatrix;
  using FullSampleMatrix = typename EigenMatrixTypes<P>::FullSampleMatrix;
  using PulseMatrix = typename EigenMatrixTypes<P>::PulseMatrix;
  using SamplePulseMatrix = typename EigenMatrixTypes<P>::SamplePulseMatrix;
  using SampleDecompLLT = typename EigenMatrixTypes<P>::SampleDecompLLT;
  using PulseDecompLDLT = typename EigenMatrixTypes<P>::PulseDecompLDLT;
  using SingleMatrix = typename EigenMatrixTypes<P>::SingleMatrix;
  using SingleVector = typename EigenMatrixTypes<P>::SingleVector;

  using Index = typename BXVector::Index;

  PulseChiSqSNNLS();
  ~PulseChiSqSNNLS();

  bool DoFit(const SampleVector &samples,
             const SampleMatrix &samplecov,
             const BXVector &bxs,
             const FullSampleVector &fullpulse,
             const FullSampleMatrix &fullpulsecov,
             const SampleGainVector &gains = -1 * SampleGainVector::Ones(),
             const SampleGainVector &badSamples = SampleGainVector::Zero());

  const SamplePulseMatrix &pulsemat() const { return _pulsemat; }
  const SampleMatrix &invcov() const { return _invcov; }

  const PulseVector &X() const { return _ampvecmin; }
  const PulseVector &Errors() const { return _errvec; }
  const BXVector &BXs() const { return _bxsmin; }

  double ChiSq() const { return _chisq; }
  void disableErrorCalculation() { _computeErrors = false; }
  void setMaxIters(int n) { _maxiters = n; }
  void setMaxIterWarnings(bool b) { _maxiterwarnings = b; }

protected:
  bool Minimize(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov);
  bool NNLS();
  void NNLSUnconstrainParameter(Index idxp);
  void NNLSConstrainParameter(Index minratioidx);
  bool OnePulseMinimize();
  bool updateCov(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov);
  double ComputeChiSq();
  double ComputeApproxUncertainty(unsigned int ipulse);

  SampleVector _sampvec;
  SampleMatrix _invcov;
  SamplePulseMatrix _pulsemat;
  PulseVector _ampvec;
  PulseVector _errvec;
  PulseVector _ampvecmin;

  SampleDecompLLT _covdecomp;
  SampleMatrix _covdecompLinv;
  PulseMatrix _topleft_work;
  PulseDecompLDLT _pulsedecomp;

  BXVector _bxs;
  BXVector _bxsmin;
  unsigned int _npulsetot;
  unsigned int _nP;

  SamplePulseMatrix invcovp;
  PulseMatrix aTamat;
  PulseVector aTbvec;
  PulseVector updatework;

  PulseVector ampvecpermtest;

  double _chisq;
  bool _computeErrors;
  int _maxiters;
  bool _maxiterwarnings;
};

#endif
