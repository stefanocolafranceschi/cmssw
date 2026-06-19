#include "RecoLocalCalo/EcalRecAlgos/interface/PulseChiSqSNNLS.h"
#include <cmath>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

template <class C>
void eigen_solve_submatrix(const typename C::PulseMatrix &mat,
                           const typename C::PulseVector &invec,
                           typename C::PulseVector &outvec,
                           const unsigned NP) {
  using namespace Eigen;
  switch (NP) {  // pulse matrix is always square.
    case 10: {
      // the template keywords are needed because the matrix types are dependent types of the template parameter C
      const auto temp = mat.template topLeftCorner<10, 10>();
      outvec.template head<10>() = temp.ldlt().solve(invec.template head<10>());
    } break;
    case 9: {
      const auto temp = mat.template topLeftCorner<9, 9>();
      outvec.template head<9>() = temp.ldlt().solve(invec.template head<9>());
    } break;
    case 8: {
      const auto temp = mat.template topLeftCorner<8, 8>();
      outvec.template head<8>() = temp.ldlt().solve(invec.template head<8>());
    } break;
    case 7: {
      const auto temp = mat.template topLeftCorner<7, 7>();
      outvec.template head<7>() = temp.ldlt().solve(invec.template head<7>());
    } break;
    case 6: {
      const auto temp = mat.template topLeftCorner<6, 6>();
      outvec.template head<6>() = temp.ldlt().solve(invec.template head<6>());
    } break;
    case 5: {
      const auto temp = mat.template topLeftCorner<5, 5>();
      outvec.template head<5>() = temp.ldlt().solve(invec.template head<5>());
    } break;
    case 4: {
      const auto temp = mat.template topLeftCorner<4, 4>();
      outvec.template head<4>() = temp.ldlt().solve(invec.template head<4>());
    } break;
    case 3: {
      const auto temp = mat.template topLeftCorner<3, 3>();
      outvec.template head<3>() = temp.ldlt().solve(invec.template head<3>());
    } break;
    case 2: {
      const auto temp = mat.template topLeftCorner<2, 2>();
      outvec.template head<2>() = temp.ldlt().solve(invec.template head<2>());
    } break;
    case 1: {
      const auto temp = mat.template topLeftCorner<1, 1>();
      outvec.template head<1>() = temp.ldlt().solve(invec.template head<1>());
    } break;
    default:
      throw cms::Exception("MultFitWeirdState")
          << "Weird number of pulses encountered in multifit, module is configured incorrectly!";
  }
}

template <class P>
PulseChiSqSNNLS<P>::PulseChiSqSNNLS() : _chisq(0.), _computeErrors(true), _maxiters(50), _maxiterwarnings(true) {}

template <class P>
PulseChiSqSNNLS<P>::~PulseChiSqSNNLS() {}

template <class P>
bool PulseChiSqSNNLS<P>::DoFit(const SampleVector &samples,
                               const SampleMatrix &samplecov,
                               const BXVector &bxs,
                               const FullSampleVector &fullpulse,
                               const FullSampleMatrix &fullpulsecov,
                               const SampleGainVector &gains,
                               const SampleGainVector &badSamples) {
  int npulse = bxs.rows();

  _sampvec = samples;
  _bxs = bxs;
  _pulsemat.resize(Eigen::NoChange, npulse);

  //construct dynamic pedestals if applicable
  int ngains = gains.maxCoeff() + 1;
  int nPedestals = 0;
  for (int gainidx = 0; gainidx < ngains; ++gainidx) {
    SampleGainVector mask = gainidx * SampleGainVector::Ones();
    SampleVector pedestal = (gains.array() == mask.array()).template cast<typename SampleVector::value_type>();
    if (pedestal.maxCoeff() > 0.) {
      ++nPedestals;
      _bxs.resize(npulse + nPedestals);
      _bxs[npulse + nPedestals - 1] = 100 + gainidx;  //bx values >=100 indicate dynamic pedestals
      _pulsemat.resize(Eigen::NoChange, npulse + nPedestals);
      _pulsemat.col(npulse + nPedestals - 1) = pedestal;
    }
  }

  //construct negative step functions for saturated or potentially slew-rate-limited samples
  for (int isample = 0; isample < SampleVector::RowsAtCompileTime; ++isample) {
    if (badSamples.coeff(isample) > 0) {
      SampleVector step = SampleVector::Zero();
      //step correction has negative sign for saturated or slew-limited samples which have been forced to zero
      step[isample] = -1.;

      ++nPedestals;
      _bxs.resize(npulse + nPedestals);
      _bxs[npulse + nPedestals - 1] =
          -100 - isample;  //bx values <=-100 indicate step corrections for saturated or slew-limited samples
      _pulsemat.resize(Eigen::NoChange, npulse + nPedestals);
      _pulsemat.col(npulse + nPedestals - 1) = step;
    }
  }

  _npulsetot = npulse + nPedestals;

  _ampvec = PulseVector::Zero(_npulsetot);
  _errvec = PulseVector::Zero(_npulsetot);
  _nP = 0;
  _chisq = 0.;

  if (_bxs.rows() == 1 && std::abs(_bxs.coeff(0)) < 100) {
    _ampvec.coeffRef(0) = _sampvec.coeff(_bxs.coeff(0) + 5);
  }

  aTamat.resize(_npulsetot, _npulsetot);

  //initialize pulse template matrix
  for (int ipulse = 0; ipulse < npulse; ++ipulse) {
    int bx = _bxs.coeff(ipulse);
    int offset = 7 - 3 - bx;
    _pulsemat.col(ipulse) = fullpulse.template segment<SampleVector::RowsAtCompileTime>(offset);
  }

  //unconstrain pedestals already for first iteration since they should always be non-zero
  if (nPedestals > 0) {
    for (int i = 0; i < _bxs.rows(); ++i) {
      int bx = _bxs.coeff(i);
      if (bx >= 100) {
        NNLSUnconstrainParameter(i);
      }
    }
  }

  //do the actual fit
  bool status = Minimize(samplecov, fullpulsecov);
  _ampvecmin = _ampvec;
  _bxsmin = _bxs;

  if (!status)
    return status;

  if (!_computeErrors)
    return status;

  //compute MINOS-like uncertainties for in-time amplitude
  bool foundintime = false;
  unsigned int ipulseintime = 0;
  for (unsigned int ipulse = 0; ipulse < _npulsetot; ++ipulse) {
    if (_bxs.coeff(ipulse) == 0) {
      ipulseintime = ipulse;
      foundintime = true;
      break;
    }
  }
  if (!foundintime)
    return status;

  const unsigned int ipulseintimemin = ipulseintime;

  double approxerr = ComputeApproxUncertainty(ipulseintime);
  double chisq0 = _chisq;
  double x0 = _ampvecmin[ipulseintime];

  //move in time pulse first to active set if necessary
  if (ipulseintime < _nP) {
    _pulsemat.col(_nP - 1).swap(_pulsemat.col(ipulseintime));
    std::swap(_ampvec.coeffRef(_nP - 1), _ampvec.coeffRef(ipulseintime));
    std::swap(_bxs.coeffRef(_nP - 1), _bxs.coeffRef(ipulseintime));
    ipulseintime = _nP - 1;
    --_nP;
  }

  SampleVector pulseintime = _pulsemat.col(ipulseintime);
  _pulsemat.col(ipulseintime).setZero();

  //two point interpolation for upper uncertainty when amplitude is away from boundary
  double xplus100 = x0 + approxerr;
  _ampvec.coeffRef(ipulseintime) = xplus100;
  _sampvec = samples - _ampvec.coeff(ipulseintime) * pulseintime;
  status &= Minimize(samplecov, fullpulsecov);
  if (!status)
    return status;
  double chisqplus100 = ComputeChiSq();

  double sigmaplus = std::abs(xplus100 - x0) / sqrt(chisqplus100 - chisq0);

  //if amplitude is sufficiently far from the boundary, compute also the lower uncertainty and average them
  if ((x0 / sigmaplus) > 0.5) {
    for (unsigned int ipulse = 0; ipulse < _npulsetot; ++ipulse) {
      if (_bxs.coeff(ipulse) == 0) {
        ipulseintime = ipulse;
        break;
      }
    }
    double xminus100 = std::max(0., x0 - approxerr);
    _ampvec.coeffRef(ipulseintime) = xminus100;
    _sampvec = samples - _ampvec.coeff(ipulseintime) * pulseintime;
    status &= Minimize(samplecov, fullpulsecov);
    if (!status)
      return status;
    double chisqminus100 = ComputeChiSq();

    double sigmaminus = std::abs(xminus100 - x0) / sqrt(chisqminus100 - chisq0);
    _errvec[ipulseintimemin] = 0.5 * (sigmaplus + sigmaminus);

  } else {
    _errvec[ipulseintimemin] = sigmaplus;
  }

  _chisq = chisq0;

  return status;
}

template <class P>
bool PulseChiSqSNNLS<P>::Minimize(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov) {
  const unsigned int npulse = _bxs.rows();

  int iter = 0;
  bool status = false;
  while (true) {
    if (iter >= _maxiters) {
      if (_maxiterwarnings) {
        LogDebug("PulseChiSqSNNLS::Minimize") << "Max Iterations reached at iter " << iter;
      }
      break;
    }

    status = updateCov(samplecov, fullpulsecov);
    if (!status)
      break;
    if (npulse > 1) {
      status = NNLS();
    } else {
      //special case for one pulse fit (performance optimized)
      status = OnePulseMinimize();
    }
    if (!status)
      break;

    double chisqnow = ComputeChiSq();
    double deltachisq = chisqnow - _chisq;

    _chisq = chisqnow;
    if (std::abs(deltachisq) < 1e-3) {
      break;
    }
    ++iter;
  }

  return status;
}

template <class P>
bool PulseChiSqSNNLS<P>::updateCov(const SampleMatrix &samplecov, const FullSampleMatrix &fullpulsecov) {
  const unsigned int nsample = SampleVector::RowsAtCompileTime;
  const unsigned int npulse = _bxs.rows();

  _invcov = samplecov;  //

  for (unsigned int ipulse = 0; ipulse < npulse; ++ipulse) {
    if (_ampvec.coeff(ipulse) == 0.)
      continue;
    int bx = _bxs.coeff(ipulse);
    if (std::abs(bx) >= 100)
      continue;  //no contribution to covariance from pedestal or saturation/slew step correction

    int firstsamplet = std::max(0, bx + 3);
    int offset = 7 - 3 - bx;

    const double ampveccoef = _ampvec.coeff(ipulse);
    const double ampsq = ampveccoef * ampveccoef;

    const unsigned int nsamplepulse = nsample - firstsamplet;
    _invcov.block(firstsamplet, firstsamplet, nsamplepulse, nsamplepulse) +=
        ampsq * fullpulsecov.block(firstsamplet + offset, firstsamplet + offset, nsamplepulse, nsamplepulse);
  }

  _covdecomp.compute(_invcov);

  bool status = true;
  return status;
}

template <class P>
double PulseChiSqSNNLS<P>::ComputeChiSq() {
  //   SampleVector resvec = _pulsemat*_ampvec - _sampvec;
  //   return resvec.transpose()*_covdecomp.solve(resvec);

  return _covdecomp.matrixL().solve(_pulsemat * _ampvec - _sampvec).squaredNorm();
}

template <class P>
double PulseChiSqSNNLS<P>::ComputeApproxUncertainty(unsigned int ipulse) {
  //compute approximate uncertainties
  //(using 1/second derivative since full Hessian is not meaningful in
  //presence of positive amplitude boundaries.)

  return 1. / _covdecomp.matrixL().solve(_pulsemat.col(ipulse)).norm();
}

template <class P>
bool PulseChiSqSNNLS<P>::NNLS() {
  //Fast NNLS (fnnls) algorithm as per http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.157.9203&rep=rep1&type=pdf

  const unsigned int npulse = _bxs.rows();
  constexpr unsigned int nsamples = SampleVector::RowsAtCompileTime;

  invcovp = _covdecomp.matrixL().solve(_pulsemat);
  aTamat.noalias() = invcovp.transpose().lazyProduct(invcovp);
  aTbvec.noalias() = invcovp.transpose().lazyProduct(_covdecomp.matrixL().solve(_sampvec));

  int iter = 0;
  Index idxwmax = 0;
  double wmax = 0.0;
  double threshold = 1e-11;
  while (true) {
    //can only perform this step if solution is guaranteed viable
    if (iter > 0 || _nP == 0) {
      if (_nP == std::min(npulse, nsamples))
        break;

      const unsigned int nActive = npulse - _nP;

      updatework = aTbvec - aTamat * _ampvec;
      Index idxwmaxprev = idxwmax;
      double wmaxprev = wmax;
      wmax = updatework.tail(nActive).maxCoeff(&idxwmax);

      //convergence
      if (wmax < threshold || (idxwmax == idxwmaxprev && wmax == wmaxprev))
        break;

      //worst case protection
      if (iter >= 500) {
        LogDebug("PulseChiSqSNNLS::NNLS()") << "Max Iterations reached at iter " << iter;
        break;
      }

      //unconstrain parameter
      Index idxp = _nP + idxwmax;
      NNLSUnconstrainParameter(idxp);
    }

    while (true) {
      //printf("iter in, idxsP = %i\n",int(_idxsP.size()));

      if (_nP == 0)
        break;

      ampvecpermtest = _ampvec;

      //solve for unconstrained parameters
      //need to have specialized function to call optimized versions
      // of matrix solver... this is truly amazing...
      eigen_solve_submatrix<EigenMatrixTypes<P>>(aTamat, aTbvec, ampvecpermtest, _nP);

      //check solution
      bool positive = true;
      for (unsigned int i = 0; i < _nP; ++i)
        positive &= (ampvecpermtest(i) > 0);
      if (positive) {
        _ampvec.head(_nP) = ampvecpermtest.head(_nP);
        break;
      }

      //update parameter vector
      Index minratioidx = 0;

      // no realizable optimization here (because it autovectorizes!)
      double minratio = std::numeric_limits<double>::max();
      for (unsigned int ipulse = 0; ipulse < _nP; ++ipulse) {
        if (ampvecpermtest.coeff(ipulse) <= 0.) {
          const double c_ampvec = _ampvec.coeff(ipulse);
          const double ratio = c_ampvec / (c_ampvec - ampvecpermtest.coeff(ipulse));
          if (ratio < minratio) {
            minratio = ratio;
            minratioidx = ipulse;
          }
        }
      }

      _ampvec.head(_nP) += minratio * (ampvecpermtest.head(_nP) - _ampvec.head(_nP));

      //avoid numerical problems with later ==0. check
      _ampvec.coeffRef(minratioidx) = 0.;

      //printf("removing index %i, orig idx %i\n",int(minratioidx),int(_bxs.coeff(minratioidx)));
      NNLSConstrainParameter(minratioidx);
    }
    ++iter;

    //adaptive convergence threshold to avoid infinite loops but still
    //ensure best value is used
    if (iter % 16 == 0) {
      threshold *= 2;
    }
  }

  return true;
}

template <class P>
void PulseChiSqSNNLS<P>::NNLSUnconstrainParameter(Index idxp) {
  aTamat.col(_nP).swap(aTamat.col(idxp));
  aTamat.row(_nP).swap(aTamat.row(idxp));
  _pulsemat.col(_nP).swap(_pulsemat.col(idxp));
  std::swap(aTbvec.coeffRef(_nP), aTbvec.coeffRef(idxp));
  std::swap(_ampvec.coeffRef(_nP), _ampvec.coeffRef(idxp));
  std::swap(_bxs.coeffRef(_nP), _bxs.coeffRef(idxp));
  ++_nP;
}

template <class P>
void PulseChiSqSNNLS<P>::NNLSConstrainParameter(Index minratioidx) {
  aTamat.col(_nP - 1).swap(aTamat.col(minratioidx));
  aTamat.row(_nP - 1).swap(aTamat.row(minratioidx));
  _pulsemat.col(_nP - 1).swap(_pulsemat.col(minratioidx));
  std::swap(aTbvec.coeffRef(_nP - 1), aTbvec.coeffRef(minratioidx));
  std::swap(_ampvec.coeffRef(_nP - 1), _ampvec.coeffRef(minratioidx));
  std::swap(_bxs.coeffRef(_nP - 1), _bxs.coeffRef(minratioidx));
  --_nP;
}

template <class P>
bool PulseChiSqSNNLS<P>::OnePulseMinimize() {
  //Fast NNLS (fnnls) algorithm as per http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.157.9203&rep=rep1&type=pdf

  //   const unsigned int npulse = 1;

  invcovp = _covdecomp.matrixL().solve(_pulsemat);
  //   aTamat = invcovp.transpose()*invcovp;
  //   aTbvec = invcovp.transpose()*_covdecomp.matrixL().solve(_sampvec);

  SingleMatrix aTamatval = invcovp.transpose() * invcovp;
  SingleVector aTbvecval = invcovp.transpose() * _covdecomp.matrixL().solve(_sampvec);
  _ampvec.coeffRef(0) = std::max(0., aTbvecval.coeff(0) / aTamatval.coeff(0));

  return true;
}

#include "DataFormats/EcalDigi/interface/EcalConstants.h"
template class PulseChiSqSNNLS<ecalPh1>;
template class PulseChiSqSNNLS<ecalPh2>;
