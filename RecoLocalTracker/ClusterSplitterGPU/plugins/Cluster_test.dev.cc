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

using namespace alpaka;
using namespace reco;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  namespace Splitting {

    template <typename TrackerTraits>
    struct JetSplit {

        // Main operator function
        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                      TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                      SiPixelDigisSoAView digiView,
                                      SiPixelClustersSoAConstView clusterView,
                                      CandidatesSoAView candidateView,
                                      ClusterGeometrysSoAView geoclusterView,
                                      double ptMin_,
                                      double deltaR_,
                                      double chargeFracMin_,
                                      float expSizeXAtLorentzAngleIncidence_,
                                      float expSizeXDeltaPerTanAlpha_,
                                      float expSizeYAtNormalIncidence_,
                                      double centralMIPCharge_,
                                      double chargePerUnit_,
                                      double fractionalWidth_,
                                      SiPixelDigisSoAView outputDigis,
                                      SiPixelClustersSoAView outputClusters,
                                      //clusterProperties* clusterPropertiesDevice,
                                      uint32_t* clusterCounterDevice,
                                      uint32_t* pixelCounterDevice,                                      
                                      double forceXError_, double forceYError_,
                                      float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                      bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset) const {

            // Get thread and grid indices
            auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]; // Thread index within the block
            auto blockIdx  = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];   // Block index
            auto blockDim  = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]; // Threads per block

            // Compute the global thread ID
            uint32_t globalThreadId = blockIdx * blockDim + threadIdx;
            uint16_t moduleId;
            uint32_t clusterOffset;

/*
            /////////////////////////////////////////////////////
            if (globalThreadId == 0) {

                // Printout the entire DigiSoA              
                for (uint32_t pixel = 0; pixel < static_cast<uint32_t>(digiView.metadata().size()); pixel++) {
                    printf("Pixel %u | clus: %d | moduleID: %u | rawIdArr: %u | adc: %u | pdigi: %u | xx: %u | yy: %u\n",
                               pixel,
                               digiView.clus(pixel),
                               digiView.moduleId(pixel),
                               digiView.rawIdArr(pixel),
                               digiView.adc(pixel),                               
                               digiView.pdigi(pixel),
                               digiView.xx(pixel),
                               digiView.yy(pixel));

                }

                // Printout the entire ClusterSoA              
                for (uint32_t pixel = 0; pixel < static_cast<uint32_t>(clusterView.metadata().size()); pixel++) {
                    printf("Cluster %u | moduleStart: %u | clusInModule: %u | moduleId: %u | clusModuleStart: %u\n",
                               pixel,
                               clusterView.moduleStart(pixel),
                               clusterView.clusInModule(pixel),
                               clusterView.moduleId(pixel),
                               clusterView.clusModuleStart(pixel) );                              
                }

                // Printout the entire HitView              
                for (uint32_t jj = 0; jj < static_cast<uint32_t>(hitView.metadata().size()); jj++) {
                    printf("hit %u | charge: %u \n",
                               jj,
                                hitView.chargeAndStatus(jj).charge);                              
                }

            }
            /////////////////////////////////////////////////////
*/

            // Get total Clusters and Candidates
            uint32_t numClusters = static_cast<uint32_t>(geoclusterView.metadata().size());
            uint32_t numCandidates = static_cast<uint32_t>(candidateView.metadata().size());

            // Ensure only valid threads process clusters
            if (globalThreadId < numClusters-2) {

                if ( globalThreadId == 0 ) {
                    *clusterCounterDevice = 0;
                    *pixelCounterDevice = 0;
                }

                uint32_t clusterIdx = globalThreadId;      // Each thread handles exactly one cluster
                moduleId = geoclusterView.moduleId(clusterIdx);
                clusterOffset = geoclusterView.clusterOffset(clusterIdx);

                if (debugMode) {
                    //uint32_t clusterIdx = 327;
                    for (uint32_t j = 0; j < static_cast<uint32_t>(digiView.metadata().size()); j++) {
                        if ( static_cast<uint32_t>(digiView.rawIdArr(j)) == static_cast<uint32_t>(targetDetId)) {
                            moduleId = digiView.moduleId(j);
                        }
                    }
                    clusterOffset = targetClusterOffset;

                    for (uint32_t j = 0; j < static_cast<uint32_t>(geoclusterView.metadata().size()); j++) {
                        if ( static_cast<uint32_t>(geoclusterView.moduleId(j)) == static_cast<uint32_t>(moduleId)) {
                            if ( static_cast<uint32_t>(geoclusterView.clusterOffset(j)) == clusterOffset) {
                                clusterIdx = j;
                            }
                        }
                    }
                }

                uint32_t begin = geoclusterView.pixelStart(clusterIdx);
                uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);                
                uint32_t pixelCounter = end - begin;
                uint32_t ClusterCharge = geoclusterView.ClusterCharge(clusterIdx);

                if (static_cast<int>(begin) < 0 || static_cast<int>(end) < 0 || static_cast<int>(end) > static_cast<int>(digiView.metadata().size())) {
                    // Avoid crash if the end is kinda wrong/overflown
                    return;
                }


/*
                // Print all about this cluster under study.........
                for (uint32_t pixel = 0; pixel < static_cast<uint32_t>(digiView.metadata().size()); pixel++) {
                    if ( static_cast<uint32_t>(digiView.moduleId(pixel)) == moduleId) {
                        if ( static_cast<uint32_t>(digiView.clus(pixel)) == clusterOffset) {
                            printf("--  clus: %d | moduleID: %u | rawIdArr: %u | adc: %u | pdigi: %u | xx: %u | yy: %u CLX: %f CLY: %f CLZ: %f\n",
                                   digiView.clus(pixel),
                                   digiView.moduleId(pixel),
                                   digiView.rawIdArr(pixel),
                                   digiView.adc(pixel),                               
                                   digiView.pdigi(pixel),
                                   digiView.xx(pixel),
                                   digiView.yy(pixel),
                                   geoclusterView.x(clusterIdx),
                                   geoclusterView.y(clusterIdx), 
                                   geoclusterView.z(clusterIdx));                            
                        }
                    }
                }
*/

                //printf("I am in thread %u, analyzing cluster %u from module %u offset %u\n", 
                //       globalThreadId, clusterOffset, moduleId, clusterOffset);


                // Access fine-tuned Global position (previously saved into the GeoCluster SoA)
                float x = geoclusterView.x(clusterIdx);
                float y = geoclusterView.y(clusterIdx);
                float z = geoclusterView.z(clusterIdx);

/*
                int ClusterCharge = 0;
                for (uint32_t i = begin; i < end; ++i) {
                    ClusterCharge += digiView.adc(i);
                }
*/
                bool doSplit = false;

                for (uint16_t candIdx = 0; candIdx < numCandidates; ++candIdx) {
                    //printf("Processing Cluster: %u, Candidate: %u/%u Block index: %u, Threads per block: %u, Total threads: %u\n",
                    //    clusterIdx, candIdx, numCandidates-1, blockIdx, blockDim, blockDim * alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

                    // Debugging Candidate to be compared to the one originated in the other producer
                    //double testme = static_cast<double>(candidateView[candIdx].px());
                    //printf("Candidate %u px= %f \n", candIdx, testme);   

                    // Skip low-pt jets
                    if (candidateView.pt(candIdx) < ptMin_) {
                        //printf("SKIP: Candidates has low pt %f \n", candidateView.pt(candIdx));
                        return;
                    }

                    // Subtract the primary vertex position to obtain the relative position
                    float relX = x - vertexX;
                    float relY = y - vertexY;
                    float relZ = z - vertexZ;
                    if (verbose_) printf("Cluster direction (cPos - vertex):");
                    if (verbose_) printf(" dx = %.3f dy = %.3f dz = %.3f\n", relX, relY, relZ);

                    // Extract jet momentum components from candidateView
                    float jetPx = candidateView.px(candIdx);
                    float jetPy = candidateView.py(candIdx);
                    float jetPz = candidateView.pz(candIdx);


                    // Compute jet transverse momentum, eta, and phi
                    float jetPt = sqrt(jetPx * jetPx + jetPy * jetPy);
                    float jetP  = sqrt(jetPx * jetPx + jetPy * jetPy + jetPz * jetPz);
                    float jetEta = 0.5 * log((jetP + jetPz) / (jetP - jetPz));
                    float jetPhi = atan2(jetPy, jetPx);

                    // Print the jet information 
                    if (verbose_) printf("In globalThreadId=%u, Jet Information:\n", globalThreadId);
                    if (verbose_) printf("  jetPx = %.3f, jetPy = %.3f, jetPz = %.3f\n", jetPx, jetPy, jetPz);
                    if (verbose_) printf("  jetPt = %.3f, jetEta = %.3f, jetPhi = %.3f\n\n", jetPt, jetEta, jetPhi);


                    // Compute the cluster's relative eta and phi
                    float r = sqrt(relX * relX + relY * relY + relZ * relZ);
                    float clusterEta = 0.5 * log((r + relZ) / (r - relZ));  // Pseudorapidity formula
                    float clusterPhi = atan2(relY, relX);  // Azimuthal angle

                    // Compute differences and deltaR (assuming 'jetEta' and 'jetPhi' are known)
                    float deltaEta = clusterEta - jetEta;
                    float deltaPhi = atan2(sin(clusterPhi - jetPhi), cos(clusterPhi - jetPhi));  // Adjust for periodicity
                    float deltaR = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
                    //printf("  deltaEta = %f, deltaPhi = %f, deltaR = %f\n", deltaEta, deltaPhi, deltaR);

                    doSplit = deltaR < deltaR_;

                    // Check deltaR condition and split clusters if applicable
                    if (doSplit) {
                        if (verbose_) printf("This clusterOffset: %u has deltaR < deltaR_ and it might be split\n",clusterOffset);

/*
                        uint32_t begin = geoclusterView.pixelStart(clusterIdx);
                        uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);
                        uint32_t pixelCounter = end - begin;

                        for (uint32_t i = begin; i < end; ++i) {
                            ClusterCharge += digiView.adc(i);
                        }
*/

/*
                        if (static_cast<int>(begin) < 0 || static_cast<int>(end) < 0 || static_cast<int>(end) > static_cast<int>(digiView.metadata().size())) {
                            if (verbose_) printf("ERROR with digi index");
                            return;
                        }
*/

                        //printf("TEST, begin=%d end=%d size=%d\n", begin, end, digiView.metadata().size());

                        if (verbose_) {
                            printf("Working on Detector Module %u clusterOffset %u with these pixels: %u\n", moduleId, clusterOffset, pixelCounter);
                            for (uint32_t i = begin; i < end; ++i) {
                                printf(" pixel adc %d x=%d y=%d\n", digiView.adc(i), digiView.xx(i), digiView.yy(i));
                            }
                        }

                        //if (pixelCounter>10000) {
                        //    printf("Pixel Exceeded 10000, closing this cluster...");
                        //    return;
                        //}

                        splitCluster(acc,
                                     hitView,
                                     digiView,
                                     clusterView,
                                     clusterIdx,
                                     moduleId,
                                     clusterOffset,
                                     jetPx, jetPy, jetPz, 
                                     geoclusterView,                                     
                                     chargeFracMin_,
                                     expSizeXAtLorentzAngleIncidence_,
                                     expSizeXDeltaPerTanAlpha_,
                                     expSizeYAtNormalIncidence_,
                                     centralMIPCharge_,
                                     chargePerUnit_,
                                     fractionalWidth_,
                                     outputDigis,
                                     outputClusters,
                                     //clusterPropertiesDevice,
                                     clusterCounterDevice,
                                     pixelCounterDevice,                                     
                                     forceXError_,
                                     forceYError_,
                                     verbose_,
                                     ClusterCharge,
                                     pixelCounter, 
                                     begin, 
                                     end);
                    }
                }
                if (!doSplit) {                       
                    storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                }
            
            }
            //else {
            //    return;
            //}
        }


        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>        
        ALPAKA_FN_ACC void printDebug(
            TAcc const& acc,
            const SiPixelDigisSoAConstView digiView,
            const SiPixelClustersSoAConstView clusterView,
            const char* info) const {

            // Debugging printout ---------------------------------
            //int CalculatedClusters = 0;

            for (int n = 0; n < static_cast<int>(clusterView.metadata().size()); n++) {
                for (uint32_t foundCluster = 0; foundCluster < clusterView.clusInModule(n); foundCluster++) {
                    for (uint32_t pixel = 0; pixel < static_cast<uint32_t>(digiView.metadata().size()); pixel++) {
                        if ( clusterView.moduleId(n) == digiView.moduleId(pixel) ) {
                            if (static_cast<int>(foundCluster) == digiView.clus(pixel)) {
                                printf("%s Module = %u ", info, clusterView.moduleId(n) );
                                printf("%s Cluster = %u ",info, foundCluster);
                                printf("%s Pixel = %u \n", info, pixel);
                            }
                        }
                    }
                }
            }
            //if (verbose_) printf("CalculatedClusters = %d ", CalculatedClusters);            
            // Debugging printout ---------------------------------
        }


        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>        
        ALPAKA_FN_ACC void storeOutputDigis(
            TAcc const& acc,
            const SiPixelDigisSoAConstView digiView,
            SiPixelDigisSoAView outputDigis,
            uint32_t begin,
            uint32_t end,            
            uint32_t* clusterCounterDevice,
            uint32_t* pixelCounterDevice) const {

            uint32_t pixelCount = end - begin;

            //if (pixelCount == 0) {
            //    printf("[storeOutputDigis] WARNING: pixelCount is zero! begin=%u end=%u\n", begin, end);
            //    return;
            //}

            // Reserve a new cluster index (used for all pixels of this cluster)
            uint32_t clusterIndex = alpaka::atomicAdd(acc, clusterCounterDevice, 1u);

            // Reserve space for the pixels
            uint32_t pixelWriteOffset = alpaka::atomicAdd(acc, pixelCounterDevice, pixelCount);

            //uint32_t maxPixels = outputDigis.metadata().size();
            //if (pixelWriteOffset + pixelCount > maxPixels) {
            //    printf("[storeOutputDigis] ERROR: Attempt to write out of bounds! "
            //           "pixelWriteOffset=%u + pixelCount=%u > maxPixels=%u\n",
            //           pixelWriteOffset, pixelCount, maxPixels);
            //    return;
            //}

            // Copy all pixels of this cluster
            for (uint32_t i = 0; i < pixelCount; ++i) {
                uint32_t srcIdx = begin + i;
                uint32_t dstIdx = pixelWriteOffset + i;

                //if (srcIdx >= static_cast<uint32_t>(digiView.metadata().size())) {
                //    printf("[storeOutputDigis] ERROR: srcIdx out of bounds! srcIdx=%u >= %d\n",
                //           srcIdx, digiView.metadata().size());
                //    continue;
                //}

                //if (dstIdx >= static_cast<uint32_t>(maxPixels)) {
                //    printf("[storeOutputDigis] ERROR: dstIdx out of bounds! dstIdx=%u >= %d\n",
                //           dstIdx, maxPixels);
                //    continue;
                //}

                outputDigis.clus(dstIdx)      = clusterIndex;
                outputDigis.xx(dstIdx)        = digiView.xx(srcIdx);
                outputDigis.yy(dstIdx)        = digiView.yy(srcIdx);
                outputDigis.adc(dstIdx)       = digiView.adc(srcIdx);
                outputDigis.rawIdArr(dstIdx)  = digiView.rawIdArr(srcIdx);
                outputDigis.moduleId(dstIdx)  = digiView.moduleId(srcIdx);

                //printf("[storeOutputDigis] Wrote pixel %u dstIdx=%u: (x=%u y=%u adc=%u rawId=%u mod=%u)\n",
                //       i, dstIdx,
                //       digiView.xx(srcIdx), digiView.yy(srcIdx), digiView.adc(srcIdx),
                //       digiView.rawIdArr(srcIdx), digiView.moduleId(srcIdx));
            }
        }


        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>        
        ALPAKA_FN_ACC void splitCluster(TAcc const& acc,
                                        TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                        SiPixelDigisSoAConstView digiView,
                                        SiPixelClustersSoAConstView clusterView,
                                        uint32_t clusterIdx,
                                        uint16_t moduleId,
                                        uint32_t clusterOffset,  
                                        float jetPx, float jetPy, float jetPz,
                                        ClusterGeometrysSoAView geoclusterView,
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
                                        bool verbose_,
                                        uint32_t ClusterCharge,
                                        uint32_t pixelCounter,
                                        uint32_t begin,
                                        uint32_t end) const {

            //if (verbose_) printf("This cluster: %u now processed in SplitCluster routine\n",clusterIdx);

            // Local variables needed for the algo
            uint16_t pixels[maxPixels];                                    
            uint16_t pixel_X[maxPixels];
            uint16_t pixel_Y[maxPixels];
            uint32_t pixel_ADC[maxPixels];
            uint32_t rawIdArr;

            float clx[maxSubClusters];
            float cly[maxSubClusters];
            float cls[maxSubClusters];
            float oldclx[maxSubClusters];
            float oldcly[maxSubClusters];

            uint16_t scoresIndices[pixelsPerCluster];
            float scoresValues[pixelsPerCluster];

            uint16_t clusterForPixel[maxPixels];

            //float weightOfPixel[maxPixels];
            bool split = false;

            float expectedADC;
            float expSizeX, expSizeY;
            {
                //float pitchX = geoclusterView.pitchX(clusterIdx);
                float pitchY = geoclusterView.pitchY(clusterIdx);

                float thickness = geoclusterView.thickness(clusterIdx);
                float tanLorentzAngles = geoclusterView.tanLorentzAngles(clusterIdx);

                // Apply precomputed transformation matrix
                float jetDirLocalX = geoclusterView.transformXX(clusterIdx) * jetPx + geoclusterView.transformYX(clusterIdx) * jetPy + geoclusterView.transformZX(clusterIdx) * jetPz;
                float jetDirLocalY = geoclusterView.transformXY(clusterIdx) * jetPx + geoclusterView.transformYY(clusterIdx) * jetPy + geoclusterView.transformZY(clusterIdx) * jetPz;
                float jetDirLocalZ = geoclusterView.transformXZ(clusterIdx) * jetPx + geoclusterView.transformYZ(clusterIdx) * jetPy + geoclusterView.transformZZ(clusterIdx) * jetPz;

                // Now, proceed with your calculations
                float jetTanAlpha = jetDirLocalX / jetDirLocalZ;
                float jetTanBeta = jetDirLocalY / jetDirLocalZ;
                float jetZOverRho = std::sqrt(jetTanAlpha * jetTanAlpha + jetTanBeta * jetTanBeta);

                expSizeX = expSizeXAtLorentzAngleIncidence_ +
                                 std::abs(expSizeXDeltaPerTanAlpha_ * (jetTanAlpha - tanLorentzAngles));
                expSizeY = std::sqrt((expSizeYAtNormalIncidence_ * expSizeYAtNormalIncidence_) +
                                           thickness * thickness / (pitchY * pitchY) * jetTanBeta * jetTanBeta);
                
                if (expSizeX < 1.f) expSizeX = 1.f;
                if (expSizeY < 1.f) expSizeY = 1.f;

                expectedADC = std::sqrt(1.08f + jetZOverRho * jetZOverRho) * centralMIPCharge_;
 

                if (verbose_) printf("Trying to split: charge=%d expSizeX=%f expSizeY=%f\n",
                        static_cast<int>(ClusterCharge), expSizeX, expSizeY);

                if ( ClusterCharge > expectedADC * chargeFracMin_ &&
                       ( ClusterCharge > expSizeX + 1 || ClusterCharge > expSizeY + 1)) {
                    split = true;
                }
            }

            if (split) {

                // Aligning to the original "fittingSplit" variables..
                uint16_t sizeX = expSizeX;
                uint16_t sizeY = expSizeY;
                uint16_t meanExp = std::floor( ClusterCharge / expectedADC + 0.5f);


                if (meanExp <= 1) {
                    if (verbose_) printf("meanExp <= 1 writing cluster");
                    storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                }
                else {
                    // Splitting the pixels and writing them for the current clusterIdx
                    if (verbose_) printf("cluster has meanExp=%d\n", meanExp);

                    uint16_t pixelsSize = 0;
                    //uint32_t begin = geoclusterView.pixelStart(clusterIdx);
                    //uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);
                    uint32_t firstOccurrence = begin;
                    rawIdArr = digiView.rawIdArr(begin);

                    int j = -1;

                    for (uint32_t jj = begin; jj < end; ++jj) {
                        j++;

                        uint16_t sub = static_cast<int>(digiView.adc(jj)) / chargePerUnit_ * expectedADC / centralMIPCharge_;
                        if (sub < 1) sub = 1;

                        uint16_t perDiv = digiView.adc(jj) / sub;

                        if (verbose_) {
                            printf("Splitting %d in [ %d , %d ], expected numb of clusters: %u original pixel (x,y) %d %d sub %d\n",
                                   j, pixelsSize, pixelsSize + sub, meanExp, digiView.xx(jj), digiView.yy(jj), sub);
                        }

                        for (uint8_t k = 0; k < sub; ++k) {
                            if (k == sub - 1) {
                                perDiv = digiView.adc(jj) - perDiv * k;
                            }

                            if (pixelsSize >= maxPixels - 1) return;

                            pixels[pixelsSize]      = j; 
                            pixel_X[pixelsSize]     = digiView.xx(jj);
                            pixel_Y[pixelsSize]     = digiView.yy(jj);
                            pixel_ADC[pixelsSize]   = perDiv;
                            //rawIdArr    = digiView.rawIdArr(jj); //moved up for perfomance
                            ++pixelsSize;
                        }
                    }

                    // Compute the initial values, set all distances and centers to -999
                    if (verbose_) printf("Computing initial values, set all distances");
                    for (uint8_t j = 0; j < meanExp; j++) {
                        oldclx[j] = -999;
                        oldcly[j] = -999;
                        clx[j] = digiView.xx(firstOccurrence) + j;
                        cly[j] = digiView.yy(firstOccurrence) + j;
                        cls[j] = 0;
                    }
                    bool stop = false;
                    uint16_t remainingSteps = 100;


                    // Refactored kernel with corrected distance scoring logic
                    while (!stop && remainingSteps > 0) {
                    //for (int step = 0; step < 100; ++step) {
                        if (verbose_) printf("---------------\n");
                        if (verbose_) printf("REMAINING STEPS : %d\n", remainingSteps);
                        remainingSteps--;

                        for (uint32_t pixelIdx = 0; pixelIdx < pixelCounter; pixelIdx++) {
                            if (pixelIdx < maxPixels) {
                                float minDist = std::numeric_limits<float>::max();
                                float secondMinDist = std::numeric_limits<float>::max();

                                uint32_t j = 0;
                                int temp_originalpixels_x = -1;
                                int temp_originalpixels_y = -1;

                                //uint32_t begin = geoclusterView.pixelStart(clusterIdx);
                                //uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);

                                for (uint32_t jj = begin; jj < end; ++jj) {
                                    if (j == pixelIdx) {
                                        temp_originalpixels_x = digiView.xx(jj);
                                        temp_originalpixels_y = digiView.yy(jj);
                                        break;
                                    }
                                    ++j;
                                }

                                // If not found, skip this pixelIdx
                                if (temp_originalpixels_x == -1 || temp_originalpixels_y == -1) {
                                    continue;
                                }

                                for (uint8_t subClusterIdx = 0; subClusterIdx < meanExp; subClusterIdx++) {
                                    float distanceX = static_cast<float>(temp_originalpixels_x) - clx[subClusterIdx];
                                    float distanceY = static_cast<float>(temp_originalpixels_y) - cly[subClusterIdx];

                                    float distX = 0.f;
                                    if (std::abs(distanceX) > sizeX / 2.f) {
                                        distX = std::pow(std::abs(distanceX) - sizeX / 2.f + 1.f, 2.f);
                                    } else {
                                        distX = std::pow(2.f * distanceX / sizeX, 2.f);
                                    }

                                    float distY = 0.f;
                                    if (std::abs(distanceY) > sizeY / 2.f) {
                                        distY = std::pow(std::abs(distanceY) - sizeY / 2.f + 1.f, 2.f);
                                    } else {
                                        distY = std::pow(2.f * distanceY / sizeY, 2.f);
                                    }

                                    float dist = std::sqrt(distX + distY);

                                    if (dist < minDist) {
                                        secondMinDist = minDist;
                                        minDist = dist;
                                    } else if (dist < secondMinDist) {
                                        secondMinDist = dist;
                                    }
                                }

                                // Score is based on -secondMinDist
                                scoresIndices[pixelIdx] = pixelIdx;
                                scoresValues[pixelIdx] = -secondMinDist;
                            }
                        }

                        // SORT SCORES ------------------------------------
                        for (uint16_t i = 0; i < pixelCounter - 1; i++) {
                            for (uint16_t j = 0; j < pixelCounter - i - 1; j++) {
                                if (scoresValues[j] > scoresValues[j + 1]) { // ascending order
                                    std::swap(scoresValues[j], scoresValues[j + 1]);
                                    std::swap(scoresIndices[j], scoresIndices[j + 1]);
                                }
                            }
                        }

                        if (verbose_) {
                            printf("Scores:\n");
                            for (uint16_t k = 0; k < pixelCounter; k++) {
                                printf("Score = %.5f, Index = %d\n", scoresValues[k], scoresIndices[k]);
                            }
                        }

                        // Caching coordinates for faster access in the following Scores Indices/Values iterations
                        float pixelX_cache[maxPixels];
                        float pixelY_cache[maxPixels];

                        // Precompute coordinates
                        for (uint16_t jj = begin, matchIdx = 0; jj < end && matchIdx < maxPixels; ++jj, ++matchIdx) {
                            pixelX_cache[matchIdx] = digiView.xx(jj);
                            pixelY_cache[matchIdx] = digiView.yy(jj);

                        }

                        // Iterating over Scores Indices and Values
                        for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
                            uint16_t pixel_index = scoresIndices[i];
                            int subpixel_counter = 0;

                            for (uint16_t subpixel = 0; subpixel < pixelsSize && subpixel < maxPixels; ++subpixel, ++subpixel_counter) {

                                if (pixels[subpixel] > static_cast<uint32_t>(pixel_index)) break;
                                if (pixels[subpixel] != static_cast<uint32_t>(pixel_index)) continue;

                                float maxEst = 0.f;
                                int cl = -1;

                                // Cached coordinates - to reenable this later
                                float temp_originalpixels_x = pixelX_cache[pixel_index];
                                float temp_originalpixels_y = pixelY_cache[pixel_index];

                                for (uint8_t subcluster_index = 0; subcluster_index < meanExp && subcluster_index < maxSubClusters; ++subcluster_index) {
                                    // Cache clx, cly and cls in registers
                                    float cx = clx[subcluster_index];
                                    float cy = cly[subcluster_index];
                                    float clusterSignal = cls[subcluster_index];

                                    float dx = temp_originalpixels_x - cx;
                                    float dy = temp_originalpixels_y - cy;

                                    float absX = std::abs(dx);
                                    float absY = std::abs(dy);

                                    float dist = 0.f;

                                    if (absX > sizeX / 2.f) {
                                        float delta = absX - sizeX / 2.f + 1.f;
                                        dist += delta * delta;
                                    } else {
                                        float norm = 2.f * dx / sizeX;
                                        dist += norm * norm;
                                    }

                                    if (absY > sizeY / 2.f) {
                                        float delta = absY - sizeY / 2.f + 1.f;
                                        dist += delta * delta;
                                    } else {
                                        float norm = 2.f * dy / sizeY;
                                        dist += norm * norm;
                                    }

                                    float distance = std::sqrt(dist);

                                    float nsig = (clusterSignal - expectedADC) / (expectedADC * fractionalWidth_);
                                    float clQest = 1.f / (1.f + std::exp(nsig)) + 1e-6f;
                                    float clDest = 1.f / (distance + 0.05f);
                                    float est = clQest * clDest;

                                    if (est > maxEst) {
                                        cl = subcluster_index;
                                        maxEst = est;
                                    }
                                }

                                // Update best-fit cluster assignment
                                cls[cl] += pixel_ADC[subpixel];
                                clusterForPixel[subpixel_counter] = cl;
                                //weightOfPixel[subpixel_counter] = maxEst;
                                //if (verbose_) printf("Pixel weight weightOfPixel[%d]=%.4f  cl=%d\n",
                                //                     subpixel_counter, weightOfPixel[subpixel_counter], cl);

                            }
                        }


                        // Recompute cluster centers
                        if (verbose_) printf("Recomputing cluster centers.........\n");

                        stop = true;
                        for (uint8_t subcluster_index = 0; subcluster_index < meanExp; subcluster_index++) {
                            //if (subcluster_index < maxSubClusters-1) {
                                if (std::abs(clx[subcluster_index] - oldclx[subcluster_index]) > 0.01f)
                                    stop = false; // still moving
                                if (std::abs(cly[subcluster_index] - oldcly[subcluster_index]) > 0.01f)
                                    stop = false;
                                oldclx[subcluster_index] = clx[subcluster_index];
                                oldcly[subcluster_index] = cly[subcluster_index];
                                clx[subcluster_index] = 0;
                                cly[subcluster_index] = 0;
                                cls[subcluster_index] = 1e-38f;//1e-99;
                            //}
                        }

                        for (uint16_t pixel_index = 0; pixel_index < pixelsSize; pixel_index++) {
                            //if (pixel_index < maxPixels-1) {
                            //    if (clusterForPixel[pixel_index] < 0)
                            //        continue;

                                clx[ clusterForPixel[pixel_index] ] += pixel_X[pixel_index] * pixel_ADC[pixel_index];
                                cly[ clusterForPixel[pixel_index] ] += pixel_Y[pixel_index] * pixel_ADC[pixel_index];
                                cls[ clusterForPixel[pixel_index] ] += pixel_ADC[pixel_index];
                            //}
                        }

                        for (uint8_t subcluster_index = 0; subcluster_index < meanExp; subcluster_index++) {
                            //if (subcluster_index < maxSubClusters-1) {                            
                                if (cls[subcluster_index] != 0) {
                                    clx[subcluster_index] /= cls[subcluster_index];
                                    cly[subcluster_index] /= cls[subcluster_index];
                                }
                                if (verbose_) printf("Center for cluster, clx[%u]=%f cly[%u]=%f\n",subcluster_index, clx[subcluster_index], subcluster_index, cly[subcluster_index]);

                                cls[subcluster_index] = 0;
                            //}
                        }
                    }
                    

                    //storeOutputDigis
                    for (uint8_t cl = 0; cl < static_cast<int>(meanExp); ++cl) {
                        // Reserve a new cluster index (only once per subcluster)
                        uint32_t clusterIndex = alpaka::atomicAdd(acc, clusterCounterDevice, 1u);

                        // Count pixels to be written for this subcluster
                        uint16_t validPixelCount = 0;
                        for (uint16_t j = 0; j < static_cast<int>(pixelsSize); ++j) {
                            if (j < maxPixels && clusterForPixel[j] == cl && pixel_ADC[j] != 0) {
                                ++validPixelCount;
                            }
                        }

                        // Reserve output space for all pixels in this subcluster
                        uint32_t pixelWriteOffset = alpaka::atomicAdd(acc, pixelCounterDevice, static_cast<uint32_t>(validPixelCount));

                        // Write pixels for this subcluster
                        int pixelOffset = 0;
                        for (uint16_t j = 0; j < static_cast<int>(pixelsSize); ++j) {
                            if (j < maxPixels && clusterForPixel[j] == cl && pixel_ADC[j] != 0) {
                                // Merge duplicate pixels (same x/y within same subcluster)
                                for (uint16_t k = j + 1; k < pixelsSize; ++k) {
                                    if (k < maxPixels &&
                                        pixel_ADC[k] != 0 &&
                                        pixel_X[k] == pixel_X[j] &&
                                        pixel_Y[k] == pixel_Y[j] &&
                                        clusterForPixel[k] == cl) {
                                            pixel_ADC[j] += pixel_ADC[k];
                                            pixel_ADC[k] = 0;
                                    }
                                }

                                // Write to output using reserved space
                                uint32_t outIdx = pixelWriteOffset + pixelOffset;

                                outputDigis.clus(outIdx)      = clusterIndex;
                                outputDigis.xx(outIdx)        = pixel_X[j];
                                outputDigis.yy(outIdx)        = pixel_Y[j];
                                outputDigis.adc(outIdx)       = pixel_ADC[j];
                                outputDigis.rawIdArr(outIdx)  = rawIdArr;
                                outputDigis.moduleId(outIdx)  = moduleId;

                                printf("NSplit cl=%d pixel_X[%d]=%d pixel_Y[%d]=%d ADC=%d \n",
                                           cl, j, pixel_X[j], j, pixel_Y[j], pixel_ADC[j]);

                                if (verbose_) {
                                    printf("Split cl=%d pixel_X[%d]=%d pixel_Y[%d]=%d ADC=%d\n",
                                           cl, j, pixel_X[j], j, pixel_Y[j], pixel_ADC[j]);
                                }
                                ++pixelOffset;
                            }
                        }
                    }
                }
            }
        }
    };



    template <typename TrackerTraits>
    void runKernels(TrackingRecHitSoAView<TrackerTraits>& hitView,
                    SiPixelDigisSoAView& digiView,
                    SiPixelClustersSoAView& clusterView,
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
                    SiPixelClustersSoAView& outputClusters,
                    //clusterProperties* clusterPropertiesDevice,
                    uint32_t* clusterCounterDevice,
                    uint32_t* pixelCounterDevice,                    
                    double forceXError_,
                    double forceYError_,
                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                    bool verbose_,
                    bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                    Queue& queue) {

    // Get the number of items per block (threads per block)
    const uint32_t threadsPerBlock = 128;

    // Calculate how many groups (blocks) you need for each view
    const uint32_t numBlocks = (geoclusterView.metadata().size() + threadsPerBlock - 1) / threadsPerBlock;
  
    //const auto MyworkDiv = make_workdiv<Acc1D>(numBlocks, threadsPerBlock);
    //const auto MyworkDiv = make_workdiv<Acc1D>(1, 1);
    const auto MyworkDiv = debugMode ? make_workdiv<Acc1D>(1, 1) : make_workdiv<Acc1D>(numBlocks, threadsPerBlock);

    if (verbose_) std::cout << "\nGot candidateView.metadata().size()=" << candidateView.metadata().size(); 
    if (verbose_) std::cout << "\nGot geoclusterView.metadata().size()=" << geoclusterView.metadata().size()
                          << "\nExecuting with " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(MyworkDiv)[0u] << " blocks and " 
                          << threadsPerBlock << " threads per block " 
                          << " and " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(MyworkDiv)[0u] 
                          << " threads in total" << std::endl;


    if (verbose_) std::cout << "In the kernel... " << std::endl;

    // std::cout << "Launching kernel with " << groups << " blocks and " << items << " threads per block." << std::endl;

                // Kernel executions AccCpuSerial should be Acc1D
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits>{}, 
                                    hitView, 
                                    digiView, 
                                    clusterView, 
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
                                    outputClusters,
                                    //clusterPropertiesDevice,
                                    clusterCounterDevice,
                                    pixelCounterDevice,                                    
                                    forceXError_,
                                    forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
            }

    // Explicit template instantiation for Phase 1
    template void runKernels<pixelTopology::Phase1>(TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clusterView,
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
                                                    SiPixelClustersSoAView& outputClusters,
                                                    //clusterProperties* clusterPropertiesDevice,
                                                    uint32_t* clusterCounterDevice,
                                                    uint32_t* pixelCounterDevice,                                                    
                                                    double forceXError_,
                                                    double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    Queue& queue);

    // Explicit template instantiation for Phase 2
    template void runKernels<pixelTopology::Phase2>(TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
                                                    SiPixelDigisSoAView& digiView,
                                                    SiPixelClustersSoAView& clusterView,
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
                                                    SiPixelClustersSoAView& outputClusters,
                                                    //clusterProperties* clusterPropertiesDevice,
                                                    uint32_t* clusterCounterDevice, 
                                                    uint32_t* pixelCounterDevice,                                                     
                                                    double forceXError_,
                                                    double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,                                                    
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
