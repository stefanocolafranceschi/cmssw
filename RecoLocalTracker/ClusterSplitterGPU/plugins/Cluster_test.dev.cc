#include <type_traits>

#include <alpaka/alpaka.hpp>

#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsSoA.h"
#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"

#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisDevice.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisSoA.h"
#include "DataFormats/SiPixelDigiSoA/interface/alpaka/SiPixelDigisSoACollection.h"

#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersDevice.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersSoA.h"
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/ZVertexDevice.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/GeometrySurface/interface/SOARotation.h"

#include "DataFormats/Math/interface/SSEVec.h"
#include "DataFormats/Math/interface/ExtVec.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometrysSoA.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrysSoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidatesSoA.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidatesSoACollection.h"

#include "Cluster_test.h"
#include "KernelShared.h"
#include "KernelRegister.h"

using namespace alpaka;
using namespace reco;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  namespace Splitting {



    template <typename TrackerTraits>
    void runKernels(//TrackingRecHitSoAView<TrackerTraits>& hitView,
                    SiPixelDigisSoAView& digiView,
                    //SiPixelClustersSoAView& clusterView,
                    CandidatesSoAView& candidateView,
                    ClusterGeometrysSoAView& geoclusterView,
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
                    //SiPixelClustersSoAView& outputClusters,
                    //clusterProperties* clusterPropertiesDevice,
                    uint32_t* clusterCounterDevice,
                    uint32_t* pixelCounterDevice,                    
                    //double forceXError_,
                    //double forceYError_,
                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                    bool verbose_,
                    bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                    uint16_t* workOnMe, 
                    uint16_t numClustersToRun, 
                    uint16_t maxPixels, 
                    Queue& queue) {

    // debugging printout
    //std::cout << "Running kernel on " << numClustersToRun << " clusters (maxPixels = " << maxPixels << ")\n";
    //std::cout << "Cluster indices to process: ";
    //for (uint16_t i = 0; i < numClustersToRun; ++i) {
    //    std::cout << workOnMe[i] << " ";
    //}
    //std::cout << std::endl;

    auto workOnMeHost = cms::alpakatools::make_host_buffer<uint16_t[]>(numClustersToRun);
    std::copy(workOnMe, workOnMe + numClustersToRun, alpaka::getPtrNative(workOnMeHost));

    auto workOnMeDevice = cms::alpakatools::make_device_buffer<uint16_t[]>(queue, numClustersToRun);
    alpaka::memcpy(queue, workOnMeDevice, workOnMeHost, numClustersToRun);
    alpaka::wait(queue); // Ensure memory copy is complete


    uint32_t threadsPerBlock;
    uint32_t numBlocks;

    if (maxPixels <= 16) {
        // Run one cluster in one thread, this launch will prioritize register memory
        // -------------------------------
        threadsPerBlock = 512;
        // Calculate how many needed blocks needed to cover all clusters
        numBlocks = (numClustersToRun + threadsPerBlock - 1) / threadsPerBlock;
    } else {
        // Run one cluster in one block, where threads work cooperatively this launch will prioritize shared memory
        // -------------------------------
        threadsPerBlock = 32;
        numBlocks = numClustersToRun;
    }

    const auto MyworkDiv = make_workdiv<Acc1D>(numBlocks, threadsPerBlock);
    //const auto MyworkDiv = make_workdiv<Acc1D>(1, 1);
    //const auto MyworkDiv = debugMode ? make_workdiv<Acc1D>(1, 1) : make_workdiv<Acc1D>(numBlocks, threadsPerBlock);

    ///if (verbose_) std::cout << "\nGot candidateView.metadata().size()=" << candidateView.metadata().size(); 
    ///if (verbose_) std::cout << "\nGot geoclusterView.metadata().size()=" << geoclusterView.metadata().size()
    ///                      << "\nExecuting with " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(MyworkDiv)[0u] << " blocks and " 
    ///                      << threadsPerBlock << " threads per block " 
    ///                      << " and " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(MyworkDiv)[0u] 
    ///                      << " threads in total" << std::endl;

    //std::cout << "Launching kernel with " << numBlocks << " blocks and " << threadsPerBlock << " threads per block." << std::endl;
    if (maxPixels < 4) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplitRegister<TrackerTraits, 4, 5, 100>{},
                                    //hitView, 
                                    digiView, 
                                    //clusterView, 
                                    candidateView, 
                                    geoclusterView,
                                    ptMin_,
                                    deltaR_,
                                    chargeFracMin_,
                                    expSizeXAtLorentzAngleIncidence_,
                                    expSizeXDeltaPerTanAlpha_,
                                    expSizeYAtNormalIncidence_,
                                    centralMIPCharge_,
                                    chargePerUnit_,
                                    fractionalWidth_,
                                    outputDigis,
                                    //outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset,
                                    alpaka::getPtrNative(workOnMeDevice), 
                                    numClustersToRun);                                    
            }
    else if (maxPixels < 8) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplitRegister<TrackerTraits, 8, 10, 200>{},
                                    //hitView, 
                                    digiView, 
                                    //clusterView, 
                                    candidateView, 
                                    geoclusterView,
                                    ptMin_,
                                    deltaR_,
                                    chargeFracMin_,
                                    expSizeXAtLorentzAngleIncidence_,
                                    expSizeXDeltaPerTanAlpha_,
                                    expSizeYAtNormalIncidence_,
                                    centralMIPCharge_,
                                    chargePerUnit_,
                                    fractionalWidth_,
                                    outputDigis,
                                    //outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset,
                                    alpaka::getPtrNative(workOnMeDevice), 
                                    numClustersToRun);                                    
            }
    else if (maxPixels < 16) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplitRegister<TrackerTraits, 16, 10, 2000>{},
                                    //hitView, 
                                    digiView, 
                                    //clusterView, 
                                    candidateView, 
                                    geoclusterView,
                                    ptMin_,
                                    deltaR_,
                                    chargeFracMin_,
                                    expSizeXAtLorentzAngleIncidence_,
                                    expSizeXDeltaPerTanAlpha_,
                                    expSizeYAtNormalIncidence_,
                                    centralMIPCharge_,
                                    chargePerUnit_,
                                    fractionalWidth_,
                                    outputDigis,
                                    //outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset,
                                    alpaka::getPtrNative(workOnMeDevice), 
                                    numClustersToRun);                                    
            }
    else if (maxPixels<32) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 32, 50, 5000>{},
                                    //hitView, 
                                    digiView, 
                                    //clusterView, 
                                    candidateView, 
                                    geoclusterView,
                                    ptMin_,
                                    deltaR_,
                                    chargeFracMin_,
                                    expSizeXAtLorentzAngleIncidence_,
                                    expSizeXDeltaPerTanAlpha_,
                                    expSizeYAtNormalIncidence_,
                                    centralMIPCharge_,
                                    chargePerUnit_,
                                    fractionalWidth_,
                                    outputDigis,
                                    //outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset,
                                    alpaka::getPtrNative(workOnMeDevice), 
                                    numClustersToRun);                                    
            }
    else if (maxPixels<256) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 256, 50, 5800>{},
                                    //hitView, 
                                    digiView, 
                                    //clusterView, 
                                    candidateView, 
                                    geoclusterView,
                                    ptMin_,
                                    deltaR_,
                                    chargeFracMin_,
                                    expSizeXAtLorentzAngleIncidence_,
                                    expSizeXDeltaPerTanAlpha_,
                                    expSizeYAtNormalIncidence_,
                                    centralMIPCharge_,
                                    chargePerUnit_,
                                    fractionalWidth_,
                                    outputDigis,
                                    //outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset,
                                    alpaka::getPtrNative(workOnMeDevice), 
                                    numClustersToRun);                                    
            }
    else {
            std::cout << "No kernel available for the given amount of pixels: " << std::endl;
        }

    }
    // Explicit template instantiation for Phase 1
    template void runKernels<pixelTopology::Phase1>(//TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    //SiPixelClustersSoAView& clusterView,
                                                    CandidatesSoAView& candidateView,
                                                    ClusterGeometrysSoAView& geoclusterView,
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
                                                    //SiPixelClustersSoAView& outputClusters,
                                                    //clusterProperties* clusterPropertiesDevice,
                                                    uint32_t* clusterCounterDevice,
                                                    uint32_t* pixelCounterDevice,                                                    
                                                    //double forceXError_,
                                                    //double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    uint16_t* workOnMe, uint16_t numClustersToRun, uint16_t maxPixels, 
                                                    Queue& queue);

    // Explicit template instantiation for Phase 2
    template void runKernels<pixelTopology::Phase2>(//TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    //SiPixelClustersSoAView& clusterView,
                                                    CandidatesSoAView& candidateView,
                                                    ClusterGeometrysSoAView& geoclusterView,
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
                                                    //SiPixelClustersSoAView& outputClusters,
                                                    //clusterProperties* clusterPropertiesDevice,
                                                    uint32_t* clusterCounterDevice, 
                                                    uint32_t* pixelCounterDevice,                                                     
                                                    //double forceXError_,
                                                    //double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,                                                    
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    uint16_t* workOnMe, uint16_t numClustersToRun, uint16_t maxPixels, 
                                                    Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
