#ifndef CondFormats_EcalObjects_EcalSampleMask_H
#define CondFormats_EcalObjects_EcalSampleMask_H
/**
 * Author: Giovanni Franzoni, UMN
 * Created: 09 Apr 2012
 * $Id: EcalSampleMask.h,v 1.1 2012/05/10 08:22:10 argiro Exp $
 **/

#include "CondFormats/Serialization/interface/Serializable.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EcalConstants.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDataFrame_Ph2.h"

#include <iostream>
#include <vector>

class EcalSampleMask {
public:
  EcalSampleMask();

  // construct from pre-organized binary words
  EcalSampleMask(const unsigned int ebmask, const unsigned int eemask);
  // constructor from an ordered set of switches, one per sample
  EcalSampleMask(const std::vector<unsigned int> &ebmask, const std::vector<unsigned int> &eemask);

  ~EcalSampleMask();

  void setEcalSampleMaskRecordEB(const unsigned int mask) { sampleMaskEB_ = mask; }
  void setEcalSampleMaskRecordEE(const unsigned int mask) { sampleMaskEE_ = mask; }
  void setEcalSampleMaskRecordEB(const std::vector<unsigned int> &ebmask);
  void setEcalSampleMaskRecordEE(const std::vector<unsigned int> &eemask);

  unsigned int getEcalSampleMaskRecordEB() const { return sampleMaskEB_; }
  unsigned int getEcalSampleMaskRecordEE() const { return sampleMaskEE_; }
  void print(std::ostream &s) const { s << "EcalSampleMask: EB " << sampleMaskEB_ << "; EE " << sampleMaskEE_; }

  bool useSampleEB(const int sampleId) const;
  bool useSampleEE(const int sampleId) const;
  bool useSample(const int sampleId, DetId &theCrystal) const;

  static constexpr int MAXSAMPLES = EcalDataFrame::MAXSAMPLES;

private:
  unsigned int sampleMaskEB_;
  unsigned int sampleMaskEE_;

  COND_SERIALIZABLE;
};

// The Phase 2 version just has a different number of samples.
// The EE variable will not be used.
class EcalPh2SampleMask : public EcalSampleMask {
public:
  EcalPh2SampleMask();
  // construct from pre-organized binary words
  EcalPh2SampleMask(const unsigned int mask);

  ~EcalPh2SampleMask() = default;

  bool useSampleEB(const int sampleId) const;
  bool useSample(const int sampleId, DetId &theCrystal) const;

  static constexpr int MAXSAMPLES = EcalDataFrame_Ph2::MAXSAMPLES;

  COND_SERIALIZABLE;
};

#endif
