#ifndef DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
#define DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h

//#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometrysSoA.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrysSoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidatesSoA.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidatesSoACollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include <alpaka/alpaka.hpp>

constexpr int maxSubClusters = 10;    //max number of resulting clusters after the split (per cluster)
constexpr int maxPixels = 400;        //virtual number of pixel "created" during the split
constexpr int pixelsPerCluster = 50;  //max number of pixel per cluster (in the original data)

using namespace reco;

namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting {

  template <typename TrackerTraits>
  void runKernels(TrackingRecHitSoAView<TrackerTraits>& hits,
                  SiPixelDigisSoAView& digis,
                  SiPixelClustersSoAView& clusters,
                  CandidatesSoAView& candidates,
                  ClusterGeometrysSoAView& geoclusters,
                  double ptMin_,
                  double deltaR_,
                  double chargeFracMin_,
                  float expSizeXAtLorentzAngleIncidence_,
                  float expSizeXDeltaPerTanAlpha_,
                  float expSizeYAtNormalIncidence_,
                  double centralMIPCharge_,
                  double chargePerUnit_,
                  double fractionalWidth_,
                  SiPixelDigisSoAView& outputDigis,
                  SiPixelClustersSoAView& outputClusters,
                  //clusterProperties* clusterPropertiesDevice,
                  uint32_t* clusterCounterDevice,
                  uint32_t* pixelCounterDevice,                  
                  double forceXError_,
                  double forceYError_,
                  float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                  bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOf,
                  Queue& queue);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::Splitting

#endif  // DataFormats_SiPixelClusterSoA_test_alpaka_Hits_test_h
