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

    template <typename TrackerTraits, uint32_t maxPixels>
    struct JetSplit {

        // Main operator function
        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                      //TrackingRecHitSoAConstView<TrackerTraits> hitView,
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
                                      //double forceXError_, double forceYError_,
                                      float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                      bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset) const {

            // Get thread and grid indices
            const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]; // Thread index within the block
            const auto blockIdx  = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];   // Block index
            const auto blockDim  = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]; // Threads per block
            //auto gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

            // Compute the global thread ID
            //uint32_t globalThreadId = blockIdx * blockDim + threadIdx;

            __attribute__((shared)) float clx[maxSubClusters];
            __attribute__((shared)) float cly[maxSubClusters];
            __attribute__((shared)) float cls[maxSubClusters];
            __attribute__((shared)) float oldclx[maxSubClusters];
            __attribute__((shared)) float oldcly[maxSubClusters];

            uint16_t pixelX_cache[maxPixels];
            uint16_t pixelY_cache[maxPixels];
            uint16_t pixelADC_cache[maxPixels];
            uint16_t pixel_info[maxPixels];

            __attribute__((shared)) uint8_t scoresIndices[maxPixels];
            __attribute__((shared)) float scoresValues[maxPixels];

            __attribute__((shared)) uint8_t clusterForPixel[extendedMaxPixels]; 

            __attribute__((shared)) uint8_t subpixelOffset[maxPixels];  //MAKE IT REGISTER
            __attribute__((shared)) bool blockShouldStop;


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
            const uint32_t numClusters = static_cast<uint32_t>(geoclusterView.metadata().size());
            const uint32_t numCandidates = static_cast<uint32_t>(candidateView.metadata().size());

            // Ensure only valid threads process clusters

                //if ( globalThreadId == 0 ) {
                //    *clusterCounterDevice = 0;
                //    *pixelCounterDevice = 0;
                //}

                //uint32_t clusterIdx = globalThreadId;      // Each thread handles exactly one cluster
                const uint32_t clusterIdx = blockIdx;            // Each block handles exactly one cluster
                //clusterIdx=651; // test 18 494 651 387;   //sample test

                //uint32_t clusterOffset = geoclusterView.clusterOffset(clusterIdx);
/*
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
*/
                const uint32_t begin = geoclusterView.pixelStart(clusterIdx);
                const uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);                
                const uint32_t pixelCounter = end - begin;
                uint32_t ClusterCharge = geoclusterView.ClusterCharge(clusterIdx);

                // NO ENOUGH SPACE ON THE TEMPORARY ARRAY
                if (pixelCounter > maxPixels) return;

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
                //       globalThreadId, clusterOffset, moduleId, clusterOffset);e

                // Access fine-tuned Global position (previously saved into the GeoCluster SoA)
                float x = geoclusterView.x(clusterIdx);
                float y = geoclusterView.y(clusterIdx);
                float z = geoclusterView.z(clusterIdx);

                // Subtract the primary vertex position to obtain the relative position
                float relX = x - vertexX;
                float relY = y - vertexY;
                float relZ = z - vertexZ;
                ///if (verbose_) printf("Cluster direction (cPos - vertex):");
                ///if (verbose_) printf(" dx = %.3f dy = %.3f dz = %.3f\n", relX, relY, relZ);

                bool doSplit = false;
                bool alreadySplit = false;

                for (uint32_t candIdx = 0; candIdx < numCandidates; ++candIdx) {


                    //printf("Processing Cluster: %u, Candidate: %u/%u Block index: %u, Threads per block: %u, Total threads: %u\n",
                    //    clusterIdx, candIdx, numCandidates-1, blockIdx, blockDim, blockDim * alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

                    // Debugging Candidate to be compared to the one originated in the other producer
                    //printf("Candidate %u px= %f \n", candIdx, static_cast<double>(candidateView[candIdx].px() );   

                    // Skip low-pt jets
                    if (candidateView.pt(candIdx) < ptMin_) {
                        //printf("SKIP: Candidates has low pt %f \n", candidateView.pt(candIdx));
                        continue;
                    }

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
                    ///if (verbose_) printf("In globalThreadId=%u, Jet Information:\n", globalThreadId);
                    ///if (verbose_) printf("  jetPx = %.3f, jetPy = %.3f, jetPz = %.3f\n", jetPx, jetPy, jetPz);
                    ///if (verbose_) printf("  jetPt = %.3f, jetEta = %.3f, jetPhi = %.3f\n\n", jetPt, jetEta, jetPhi);

                    // Compute the cluster's relative eta and phi
                    float r = sqrt(relX * relX + relY * relY + relZ * relZ);
                    float clusterEta = 0.5 * log((r + relZ) / (r - relZ));  // Pseudorapidity formula
                    float clusterPhi = atan2(relY, relX);  // Azimuthal angle

                    // Compute differences and deltaR (assuming 'jetEta' and 'jetPhi' are known)
                    float deltaEta = clusterEta - jetEta;
                    float deltaPhi = atan2(sin(clusterPhi - jetPhi), cos(clusterPhi - jetPhi));  // Adjust for periodicity
                    float deltaR = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

                    doSplit = deltaR < deltaR_;

                    // Check deltaR condition and split clusters if applicable
                    if (doSplit) {
                        ///if (verbose_) printf("This clusterOffset: %u has deltaR < deltaR_ and it might be split\n",clusterOffset);

                        ///if (verbose_) {
                        ///    printf("Working on Detector Module %u clusterOffset %u with these pixels: %u\n", moduleId, clusterOffset, pixelCounter);
                        ///    for (uint32_t i = begin; i < end; ++i) {
                        ///        printf(" pixel adc %d x=%d y=%d\n", digiView.adc(i), digiView.xx(i), digiView.yy(i));
                        ///    }
                        ///}

                        ///if (verbose_) printf("This cluster: %u now processed in SplitCluster routine\n",clusterIdx);

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
             
                            //printf("Trying to split: charge=%d expSizeX=%f expSizeY=%f\n",
                            //        static_cast<int>(ClusterCharge), expSizeX, expSizeY);

                            if ( ClusterCharge > expectedADC * chargeFracMin_ &&
                                   ( geoclusterView.sizeX(clusterIdx) > expSizeX + 1 || geoclusterView.sizeY(clusterIdx) > expSizeY + 1)) {
                                split = true;
                            }
                        }

                        if (split) {


                            // Filling local cache for faster access and sub (for avoiding repeated pixels large arrays)
                            for (uint16_t jj = begin, matchIdx = 0; jj < end && matchIdx < maxPixels; ++jj, ++matchIdx) {
                                pixelX_cache[matchIdx] = digiView.xx(jj);
                                pixelY_cache[matchIdx] = digiView.yy(jj);
                                pixelADC_cache[matchIdx] = digiView.adc(jj);

                                uint8_t sub = static_cast<int>(pixelADC_cache[matchIdx]) / chargePerUnit_ * expectedADC / centralMIPCharge_;
                                if (sub < 1) sub = 1;
                                pixel_info[matchIdx] = sub;
                                //printf("%u x=%u y=%u c=%u, sub=%u \n", matchIdx, pixelX_cache[matchIdx], pixelY_cache[matchIdx], pixelADC_cache[matchIdx], pixel_info[matchIdx] );
                            }

/*
                            // Filling local cache for faster access and sub (for avoiding repeated pixels large arrays)
                            for (uint16_t idx = threadIdx; idx < pixelCounter && idx < maxPixels; idx += blockDim) {
                                uint16_t jj = begin + idx;
                                pixelX_cache[idx] = digiView.xx(jj);
                                pixelY_cache[idx] = digiView.yy(jj);
                                uint16_t charge = digiView.adc(jj);
                                pixelADC_cache[idx] = charge;

                                uint16_t sub = static_cast<int>( charge) / chargePerUnit_ * expectedADC / centralMIPCharge_;
                                
                                if (sub < 1) sub = 1;
                                pixel_info[idx] = sub;
                                //printf("%u x=%u y=%u c=%u, sub=%u \n", idx, pixelX_cache[idx], pixelY_cache[idx], pixelADC_cache[idx], pixel_info[idx] );                                
                            }
                            alpaka::syncBlockThreads(acc);                            
                            //-----------------------------------
*/

                            // Aligning to the original "fittingSplit" variables..
                            uint16_t sizeX = expSizeX;
                            uint16_t sizeY = expSizeY;
                            uint8_t meanExp = std::floor( ClusterCharge / expectedADC + 0.5f);

                            if (meanExp <= 1) {
                                ///if (verbose_) printf("meanExp <= 1 writing cluster");
                                storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                            }
                            else {
                                // Splitting the pixels and writing them for the current clusterIdx
                                ///if (verbose_) printf("cluster has meanExp=%d\n", meanExp);

                                // Compute the initial values, set all distances and centers to -999
                                ///if (verbose_) printf("Computing initial values, set all distances");

/*
                                // Resetting temp variables
                                for (uint16_t j = 0; j < meanExp; j++) {
                                    oldclx[j] = -999;
                                    oldcly[j] = -999;
                                    clx[j] = pixelX_cache[0] + j;
                                    cly[j] = pixelY_cache[0] + j;
                                    cls[j] = 0;
                                }
*/

                                // zeroing
                                // Parallel reset of clx, cly, cls, oldclx, oldcly
                                for (uint16_t i = threadIdx; i < maxSubClusters; i += blockDim) {
                                    clx[i] = 0;
                                    cly[i] = 0;
                                    cls[i] = 0;
                                    oldclx[i] = 0;
                                    oldcly[i] = 0;
                                }

                                // Parallel reset of clusterForPixel
                                for (uint16_t i = threadIdx; i < extendedMaxPixels; i += blockDim) {
                                    clusterForPixel[i] = 0;
                                }

                                alpaka::syncBlockThreads(acc);  // Ensure reset is complete before continuing

                                // Resetting temp variables (initialization)
                                for (uint16_t j = threadIdx; j < meanExp; j += blockDim) {
                                    oldclx[j] = -999;
                                    oldcly[j] = -999;
                                    clx[j]    = pixelX_cache[0] + j;
                                    cly[j]    = pixelY_cache[0] + j;
                                    cls[j]    = 0;
                                }

                                alpaka::syncBlockThreads(acc);  // Again ensure sync after re-initialization


                                // Main recursive algorithm ----------------
                                blockShouldStop = false;
                                for (uint8_t remainingSteps = 100; remainingSteps > 0 && !blockShouldStop; --remainingSteps) {

                                    // Reseting temporary variables
                                    alpaka::syncBlockThreads(acc);
                                    for (uint16_t i = threadIdx; i < maxPixels; i += blockDim) {
                                        scoresIndices[i] = 0;
                                        scoresValues[i] = 0.f;
                                    }
                                    alpaka::syncBlockThreads(acc);

                                    //if (verbose_) printf("---------------\n");
                                    //if (verbose_) printf("REMAINING STEPS : %d\n", remainingSteps);

/*
                                    for (uint16_t pixelIdx = 0; pixelIdx < pixelCounter; pixelIdx++) {
                                        if (pixelIdx < maxPixels) {
                                            float minDist = std::numeric_limits<float>::max();
                                            float secondMinDist = std::numeric_limits<float>::max();

                                            uint16_t j = 0;
                                            float temp_originalpixels_x = -1;
                                            float temp_originalpixels_y = -1;

                                            for (uint16_t jj = 0; jj < pixelCounter; jj++) {
                                                if (j == pixelIdx) {
                                                    temp_originalpixels_x = pixelX_cache[jj];
                                                    temp_originalpixels_y = pixelY_cache[jj];
                                                    break;
                                                }
                                                ++j;
                                            }

                                            // If not found, skip this pixelIdx
                                            if (temp_originalpixels_x == -1 || temp_originalpixels_y == -1) continue;

                                            for (uint16_t subClusterIdx = 0; subClusterIdx < meanExp; subClusterIdx++) {
                                                float distanceX = static_cast<float>(temp_originalpixels_x) - clx[subClusterIdx];
                                                float distanceY = static_cast<float>(temp_originalpixels_y) - cly[subClusterIdx];

                                                //printf("Pixel %d: clx[%d]=%.4f cly[%d]=%.4f original_x=%f original_y=%f\n",
                                                //       pixelIdx, subClusterIdx, clx[subClusterIdx], subClusterIdx, cly[subClusterIdx],
                                                //       temp_originalpixels_x, temp_originalpixels_y);

                                                float distX = 0.f;
                                                if (std::abs(distanceX) > sizeX / 2.f) {
                                                    float diff = std::abs(distanceX) - sizeX / 2.f + 1.f;
                                                    distX = diff * diff;
                                                } else {
                                                    float scaled = 2.f * distanceX / sizeX;
                                                    distX = scaled * scaled;
                                                }

                                                float distY = 0.f;
                                                if (std::abs(distanceY) > sizeY / 2.f) {
                                                    float diff = std::abs(distanceY) - sizeY / 2.f + 1.f;
                                                    distY = diff * diff;
                                                } else {
                                                    float scaled = 2.f * distanceY / sizeY;
                                                    distY = scaled * scaled;
                                                }

                                                float dist = std::sqrt(distX + distY);
                                                //printf("subClusterIdx=%u distX=%f distanceY=%f dist=%f\n", subClusterIdx, distanceX, distanceY, dist );

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
*/

                                    for (uint16_t pixelIdx = threadIdx; pixelIdx < pixelCounter; pixelIdx += blockDim) {
                                        if (pixelIdx < maxPixels) {
                                            float minDist = std::numeric_limits<float>::max();
                                            float secondMinDist = std::numeric_limits<float>::max();

                                            // Find the (x, y) of this pixelIdx
/*
                                            float temp_originalpixels_x = -1.f;
                                            float temp_originalpixels_y = -1.f;
                                            
                                            for (uint16_t jj = 0; jj < pixelCounter; ++jj) {
                                                if (jj == pixelIdx) {
                                                    temp_originalpixels_x = pixelX_cache[jj];
                                                    temp_originalpixels_y = pixelY_cache[jj];
                                                    break;
                                                }
                                            }
*/
                                            float temp_originalpixels_x = pixelX_cache[pixelIdx];
                                            float temp_originalpixels_y = pixelY_cache[pixelIdx];


                                            // If not found, skip this pixelIdx
                                            if (temp_originalpixels_x == -1.f || temp_originalpixels_y == -1.f) continue;

                                            // Compute distances to all subclusters
                                            for (uint8_t subClusterIdx = 0; subClusterIdx < meanExp; ++subClusterIdx) {
                                                float distanceX = temp_originalpixels_x - clx[subClusterIdx];
                                                float distanceY = temp_originalpixels_y - cly[subClusterIdx];

                                                float distX = 0.f;
                                                if (std::abs(distanceX) > sizeX / 2.f) {
                                                    float diff = std::abs(distanceX) - sizeX / 2.f + 1.f;
                                                    distX = diff * diff;
                                                } else {
                                                    float scaled = 2.f * distanceX / sizeX;
                                                    distX = scaled * scaled;
                                                }

                                                float distY = 0.f;
                                                if (std::abs(distanceY) > sizeY / 2.f) {
                                                    float diff = std::abs(distanceY) - sizeY / 2.f + 1.f;
                                                    distY = diff * diff;
                                                } else {
                                                    float scaled = 2.f * distanceY / sizeY;
                                                    distY = scaled * scaled;
                                                }

                                                float dist = std::sqrt(distX + distY);

                                                if (dist < minDist) {
                                                    secondMinDist = minDist;
                                                    minDist = dist;
                                                } else if (dist < secondMinDist) {
                                                    secondMinDist = dist;
                                                }
                                            }

                                            // Save score: negative secondMinDist
                                            scoresIndices[pixelIdx] = pixelIdx;
                                            scoresValues[pixelIdx] = -secondMinDist;
                                        }
                                    }
                                    alpaka::syncBlockThreads(acc);                            
                                    //-----------------------------------

                                    // SORT SCORES ------------------------------------
                                    // (not much gain in parallel, so sequential)
                                    if (threadIdx==0) {
                                        for (uint16_t i = 0; i < pixelCounter - 1; i++) {
                                            for (uint16_t j = 0; j < pixelCounter - i - 1; j++) {
                                                if (scoresValues[j] > scoresValues[j + 1] ||
                                                    (scoresValues[j] == scoresValues[j + 1] && scoresIndices[j] > scoresIndices[j + 1])) {
                                                    std::swap(scoresValues[j], scoresValues[j + 1]);
                                                    std::swap(scoresIndices[j], scoresIndices[j + 1]);
                                                }
                                            }
                                        }

                                        if (verbose_) {
                                            printf("Cluster %u Scores:\n", clusterIdx);
                                            for (uint16_t k = 0; k < pixelCounter; k++) {
                                                printf("Cluster %u Score = %.5f, Index = %d\n", clusterIdx, scoresValues[k], scoresIndices[k]);
                                            }
                                        }
                                    }
                                    alpaka::syncBlockThreads(acc);                            
                                    //-----------------------------------



                                    /////////////////////////////////////////////////////////////////
                                    // Thread 0 computes subpixel offset for each pixel
                                    if (threadIdx == 0) {
                                        uint32_t offset = 0;
                                        for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
                                            subpixelOffset[i] = offset;
                                            offset += pixel_info[i];  // accumulate subpixels
                                        }
                                    }
                                    alpaka::syncBlockThreads(acc);  // Make sure subpixelOffset[] is visible to all


                                    // Iterating over Scores Indices and Values                                    
                                    // Each thread handles one score index
                                    for (uint16_t i = threadIdx; i < pixelCounter && i < maxPixels; i += blockDim) {
                                        uint16_t pixel_index = scoresIndices[i];

                                        uint8_t sub = pixel_info[pixel_index];
                                        uint16_t adc = pixelADC_cache[pixel_index];
                                        uint16_t perDiv = adc / sub;

                                        float temp_x = pixelX_cache[pixel_index];
                                        float temp_y = pixelY_cache[pixel_index];

                                        for (uint8_t k = 0; k < sub; ++k) {

                                            float maxEst = 0.f;
                                            int cl = -1;

                                            for (uint16_t subcluster_index = 0; subcluster_index < meanExp && subcluster_index < maxSubClusters; ++subcluster_index) {
                                                //printf("pixelCounter=%u subcluster_index=%u clx[%u]=%f cly[%u]=%f cls[%u]=%f sizeX=%u sizeY=%u\n ",pixel_index,subcluster_index,subcluster_index,clx[subcluster_index],subcluster_index,cly[subcluster_index],subcluster_index,cls[subcluster_index],sizeX,sizeY);

                                                float cx = clx[subcluster_index];
                                                float cy = cly[subcluster_index];
                                                float clusterSignal = cls[subcluster_index];

                                                float dx = temp_x - cx;
                                                float dy = temp_y - cy;

                                                float dist = 0.f;
                                                float absX = std::abs(dx), absY = std::abs(dy);

                                                if (absX > sizeX / 2.f)
                                                    dist += (absX - sizeX / 2.f + 1.f) * (absX - sizeX / 2.f + 1.f);
                                                else
                                                    dist += (2.f * dx / sizeX) * (2.f * dx / sizeX);

                                                if (absY > sizeY / 2.f)
                                                    dist += (absY - sizeY / 2.f + 1.f) * (absY - sizeY / 2.f + 1.f);
                                                else
                                                    dist += (2.f * dy / sizeY) * (2.f * dy / sizeY);

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

                                            uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;

                                            // Atomically add charge to cluster total
                                            alpaka::atomicAdd(acc, &cls[cl], static_cast<float>(charge));

                                            // Write cluster assignment
                                            clusterForPixel[subpixelOffset[pixel_index] + k] = cl;

                                            //printf("remainingSteps=%u GPU: pixel_index=%u k=%u -> subpixel=%u -> cl=%d "
                                            //        "charge=%u est=%.4f (maxEst=%.4f) pixel=(%.2f,%.2f) -> cl center=(%.2f,%.2f)\n",
                                            //        remainingSteps, pixel_index, k, subpixelOffset[pixel_index] + k, cl,
                                            //        charge, maxEst, maxEst, temp_x, temp_y, clx[cl], cly[cl]);
                                            //printf("remainingSteps=%u GPU: pixel_index=%u k=%u -> subpixel=%u -> cl=%d charge=%u\n",
                                            //        remainingSteps, pixel_index, k, subpixelOffset[pixel_index] + k, cl, charge);
                                        }
                                    }
                                    alpaka::syncBlockThreads(acc);


                                    // Recompute cluster centers
                                    if (verbose_) printf("Recomputing cluster centers.........\n");
                                    if (threadIdx == 0) blockShouldStop = true;  // Assume we will stop
                                    alpaka::syncBlockThreads(acc);


                                    if (threadIdx == 0) {
                                        //for (uint16_t oo = 0; oo < meanExp; oo++) {
                                        //    printf("remainingSteps=%u oldclx[%u]=%f oldcly[%u]=%f clx[%u]=%f cly[%u]=%f\n",remainingSteps, oo,oldclx[oo],oo,oldcly[oo],oo,clx[oo],oo,oldcly[oo]);
                                        //}
                                    }

                                    // Each thread checks part of the subcluster array
                                    for (uint16_t i = threadIdx; i < meanExp; i += blockDim) {
                                        if (std::abs(clx[i] - oldclx[i]) > 0.01f || std::abs(cly[i] - oldcly[i]) > 0.01f) {
                                            // Someone observed movement; mark as not ready to stop
                                            blockShouldStop = false;
                                        }

                                        // Update old cluster centers
                                        oldclx[i] = clx[i];
                                        oldcly[i] = cly[i];

                                        // Reset centers and signal
                                        clx[i] = 0.f;
                                        cly[i] = 0.f;
                                        cls[i] = 1e-38f;
                                    }

                                    alpaka::syncBlockThreads(acc);
                                    /////////////////////////////////////////////////////////////////





                                    /*
                                    /////////////////////////////////////////////////////////////////
                                    if (threadIdx==0) {

                                        // Iterating over Scores Indices and Values
                                        for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
                                            uint16_t pixel_index = scoresIndices[i];
                                            int subpixel_counter=-1;


                                            for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
          
                                                // Subpixel "simulation"
                                                uint16_t sub = pixel_info[i];
                                                uint16_t adc = pixelADC_cache[i];
                                                uint16_t perDiv = adc / sub;

                                                for (uint16_t k = 0; k < sub; ++k) {

                                                    subpixel_counter++;

                                                    if (i > static_cast<uint32_t>(pixel_index)) break;
                                                    if (i != static_cast<uint32_t>(pixel_index)) continue;

                                                    float maxEst = 0.f;
                                                    int cl = -1;

                                                    float temp_x = pixelX_cache[i];
                                                    float temp_y = pixelY_cache[i];

                                                    for (uint16_t subcluster_index = 0; subcluster_index < meanExp && subcluster_index < maxSubClusters; ++subcluster_index) {
                                                        printf("pixelCounter=%u subcluster_index=%u clx[%u]=%f cly[%u]=%f cls[%u]=%f sizeX=%u sizeY=%u\n ",i,subcluster_index,subcluster_index,clx[subcluster_index],subcluster_index,cly[subcluster_index],subcluster_index,cls[subcluster_index],sizeX,sizeY);

                                                        float cx = clx[subcluster_index];
                                                        float cy = cly[subcluster_index];
                                                        float clusterSignal = cls[subcluster_index];

                                                        float dx = temp_x - cx;
                                                        float dy = temp_y - cy;

                                                        float dist = 0.f;
                                                        float absX = std::abs(dx), absY = std::abs(dy);

                                                        if (absX > sizeX / 2.f)
                                                            dist += (absX - sizeX / 2.f + 1.f) * (absX - sizeX / 2.f + 1.f);
                                                        else
                                                            dist += (2.f * dx / sizeX) * (2.f * dx / sizeX);

                                                        if (absY > sizeY / 2.f)
                                                            dist += (absY - sizeY / 2.f + 1.f) * (absY - sizeY / 2.f + 1.f);
                                                        else
                                                            dist += (2.f * dy / sizeY) * (2.f * dy / sizeY);

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

                                                    uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;
                                                    cls[cl] += charge;
                                                    clusterForPixel[subpixel_counter] = cl;

                                                    //printf("remainingSteps=%u GPU: pixel_index=%u k=%u -> subpixel=%u -> cl=%d "
                                                    //        "charge=%u est=%.4f (maxEst=%.4f) pixel=(%.2f,%.2f) -> cl center=(%.2f,%.2f)\n",
                                                    //        remainingSteps, pixel_index, k, subpixel_counter, cl,
                                                    //        charge, maxEst, maxEst, temp_x, temp_y, clx[cl], cly[cl]);

                                                    //printf("remainingSteps=%u SER: pixel_index=%u k=%u -> subpixel=%u -> cl=%d charge=%u\n",
                                                    //    remainingSteps, pixel_index, k, subpixelOffset[pixel_index] + k, cl, charge);
                                                }
                                            }
                                        }

                                        for (uint16_t oo = 0; oo < meanExp; oo++) {
                                            printf("remainingSteps=%u oldclx[%u]=%f oldcly[%u]=%f clx[%u]=%f cly[%u]=%f\n",remainingSteps, oo,oldclx[oo],oo,oldcly[oo],oo,clx[oo],oo,oldcly[oo]);
                                        }


                                        // Recompute cluster centers
                                        if (verbose_) printf("Recomputing cluster centers.........\n");

                                        blockShouldStop = true;
                                        for (uint16_t subcluster_index = 0; subcluster_index < meanExp; subcluster_index++) {
                                            if (std::abs(clx[subcluster_index] - oldclx[subcluster_index]) > 0.01f)
                                                blockShouldStop = false; // still moving
                                            if (std::abs(cly[subcluster_index] - oldcly[subcluster_index]) > 0.01f)
                                                blockShouldStop = false;
                                            oldclx[subcluster_index] = clx[subcluster_index];
                                            oldcly[subcluster_index] = cly[subcluster_index];
                                            clx[subcluster_index] = 0;
                                            cly[subcluster_index] = 0;
                                            cls[subcluster_index] = 1e-38f;//1e-99;
                                        }
                                    }
                                    /////////////////////////////////////////////////////////////////
                                    */




if (threadIdx==0) {

                                    //for (uint16_t i = 0; i < extendedMaxPixels; ++i) {
                                    //    printf("i=%d clusterForPixel=%d\n",i,clusterForPixel[i]);
                                    //}

                                    int nnn = 0 ;
                                    for (uint16_t i = 0; i < pixelCounter && i < maxPixels; i++) {

                                        uint16_t x = pixelX_cache[i];
                                        uint16_t y = pixelY_cache[i];                                

                                        // Subpixel "simulation"
                                        uint8_t sub = pixel_info[i];
                                        uint16_t adc = pixelADC_cache[i];
                                        uint16_t perDiv = adc / sub;

                                        for (uint8_t k = 0; k < sub; k++) {
                                            uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;
                                            //printf("nnn=%d x*charge=%d  y*charge=%d  clx=%f  cly=%f  cls=%f\n", nnn, x * charge, y * charge, clx[ clusterForPixel[nnn] ],cly[ clusterForPixel[nnn] ],cls[ clusterForPixel[nnn] ]);

                                            clx[ clusterForPixel[nnn] ] += x * charge;
                                            cly[ clusterForPixel[nnn] ] += y * charge;
                                            cls[ clusterForPixel[nnn] ] += charge;
                                            nnn++;
                                        }
                                    }

                                    for (uint8_t subcluster_index = 0; subcluster_index < meanExp; subcluster_index++) {
                                        if (cls[subcluster_index] != 0) {
                                            clx[subcluster_index] /= cls[subcluster_index];
                                            cly[subcluster_index] /= cls[subcluster_index];
                                        }
                                        //printf("Center for cluster, clx[%u]=%f cly[%u]=%f\n",subcluster_index, clx[subcluster_index], subcluster_index, cly[subcluster_index]);
                                        cls[subcluster_index] = 0;
                                    }
}
                                    alpaka::syncBlockThreads(acc);  // Make sure all threads completed and wrote to blockShouldStop
                                }



                                if (threadIdx==0) {          
                                    //storeOutputDigis
                                    uint32_t kkk = 0;
                                    for (uint16_t cl = 0; cl < static_cast<int>(meanExp); ++cl) {

                                        uint32_t clusterIndex = alpaka::atomicAdd(acc, clusterCounterDevice, 1u);
                                        uint32_t pixelStartWritingAt = alpaka::atomicAdd(acc, pixelCounterDevice, 0u);
                                        uint32_t pixelOffset = 0;

                                        int nnn = 0;

                                        for (uint16_t i = 0; i < pixelCounter && i < maxPixels; i++) {

                                            uint16_t x = pixelX_cache[i];
                                            uint16_t y = pixelY_cache[i];                                

                                            uint8_t sub = pixel_info[i];
                                            uint16_t adc = pixelADC_cache[i];
                                            uint16_t perDiv = adc / sub;

                                            uint32_t writeCharge = 0;

                                            for (uint8_t k = 0; k < sub; k++, nnn++) {
                                                if (clusterForPixel[nnn] == cl) {
                                                    uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;
                                                    writeCharge += charge;
                                                }
                                            }

                                            // Only write if this pixel contributed subpixels to cluster cl
                                            if (writeCharge > 0) {
                                                uint32_t outIdx = pixelStartWritingAt + pixelOffset;
                                                uint32_t rawIdArr = digiView.rawIdArr(begin);

                /*
                                                if (kkk > (end - begin)) {
                                                    printf("Error, more clusters than original");
                                                    return;
                                                }

                                                outputDigis.clus(outIdx)      = clusterIndex;
                                                outputDigis.xx(outIdx)        = x;
                                                outputDigis.yy(outIdx)        = y;
                                                outputDigis.adc(outIdx)       = writeCharge;
                                                outputDigis.rawIdArr(outIdx)  = rawIdArr;
                                                outputDigis.moduleId(outIdx)  = moduleId;

                                                alpaka::atomicAdd(acc, pixelCounterDevice, 1u);
                */
                                                if (verbose_) {
                                                    uint16_t moduleId = geoclusterView.moduleId(clusterIdx);
                                                    printf("candIdx=%u/%u moduleId=%u NSplit cl=%d rawIdArr %d pixel_X[%d]=%u pixel_Y[%d]=%u ADC=%d \n",
                                                       candIdx,numCandidates, moduleId, cl, rawIdArr, i, x, i, y, writeCharge);
                                                }
                                                pixelOffset++;
                                                kkk++;
                                            }
                                        }
                                    }
                                    return;
                                }
                               



                            }

                        }
                        alreadySplit = true;

                    }
                }

                if (threadIdx == 0) {    
                    if (!doSplit && !alreadySplit) {                       
                        storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                    }
                }            
            //}
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

            uint32_t MaxPixels = outputDigis.metadata().size();
            uint32_t currentPixelCount = alpaka::atomicAdd(acc, pixelCounterDevice, 0u);  // peek current

            // NO ENOUGH SPACE ON THE OUTPUT SOA
            if (currentPixelCount + pixelCount > MaxPixels) {
                // Log and return safely
                //printf("[storeOutputDigis] ERROR: Attempt to write out of bounds! "
                //       "currentPixelCount=%u + pixelCount=%u > maxPixels=%u\n",
                //       currentPixelCount, pixelCount, maxPixels);
                return;
            }

            // Reserve a new cluster index (used for all pixels of this cluster)
            uint32_t clusterIndex = alpaka::atomicAdd(acc, clusterCounterDevice, 1u);

            // Reserve space for the pixels
            uint32_t pixelWriteOffset = alpaka::atomicAdd(acc, pixelCounterDevice, pixelCount);

            // Copy all pixels of this cluster
            for (uint32_t i = 0; i < pixelCount; i++) {
                uint32_t srcIdx = begin + i;
                uint32_t dstIdx = currentPixelCount + i;

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
    };



    template <typename TrackerTraits>
    void runKernels(//TrackingRecHitSoAView<TrackerTraits>& hitView,
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
                    //double forceXError_,
                    //double forceYError_,
                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                    bool verbose_,
                    bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                    uint16_t maxPixelsRetrieved,
                    Queue& queue) {

    // Get the number of items per block (threads per block)
    const uint32_t threadsPerBlock = 32;

    // Calculate how many groups (blocks) you need for each view
    //const uint32_t numBlocks = (geoclusterView.metadata().size() + threadsPerBlock - 1) / threadsPerBlock;
    const uint32_t numBlocks = geoclusterView.metadata().size();
  
    //const auto MyworkDiv = make_workdiv<Acc1D>(numBlocks, threadsPerBlock);
    //const auto MyworkDiv = make_workdiv<Acc1D>(1, 1);
    const auto MyworkDiv = debugMode ? make_workdiv<Acc1D>(1, 1) : make_workdiv<Acc1D>(numBlocks, threadsPerBlock);

    ///if (verbose_) std::cout << "\nGot candidateView.metadata().size()=" << candidateView.metadata().size(); 
    ///if (verbose_) std::cout << "\nGot geoclusterView.metadata().size()=" << geoclusterView.metadata().size()
    ///                      << "\nExecuting with " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(MyworkDiv)[0u] << " blocks and " 
    ///                      << threadsPerBlock << " threads per block " 
    ///                      << " and " << alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(MyworkDiv)[0u] 
    ///                      << " threads in total" << std::endl;


    ///if (verbose_) std::cout << "In the kernel... " << std::endl;
    //std::cout << "MaxPixelRetrieved " << maxPixelsRetrieved << std::endl;
    // std::cout << "Launching kernel with " << groups << " blocks and " << items << " threads per block." << std::endl;

/*
    if (maxPixelsRetrieved<32) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 32>{},
                                    //hitView, 
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
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
            }
    else if (maxPixelsRetrieved<64) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 64>{},
                                    //hitView, 
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
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
            }
    else if (maxPixelsRetrieved<128) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 128>{},
                                    //hitView, 
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
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
            }
    else if (maxPixelsRetrieved<256) {
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 256>{},
                                    //hitView, 
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
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
            }
    else if (maxPixelsRetrieved<512) {
*/
                alpaka::exec<Acc1D>(queue, 
                                    MyworkDiv, 
                                    JetSplit<TrackerTraits, 512>{},
                                    //hitView, 
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
                                    //forceXError_,
                                    //forceYError_,
                                    vertexX, vertexY, vertexZ, vertexEta, vertexPhi, 
                                    verbose_, debugMode, targetDetId, targetClusterOffset);
/*
            }
    else {
            std::cout << "No kernel available for the given amount of pixels: " << maxPixelsRetrieved << std::endl;
        }
*/

    }
    // Explicit template instantiation for Phase 1
    template void runKernels<pixelTopology::Phase1>(//TrackingRecHitSoAView<pixelTopology::Phase1>& hitView,
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
                                                    //double forceXError_,
                                                    //double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    uint16_t maxPixelsRetrieved,
                                                    Queue& queue);

    // Explicit template instantiation for Phase 2
    template void runKernels<pixelTopology::Phase2>(//TrackingRecHitSoAView<pixelTopology::Phase2>& hitView,
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
                                                    //double forceXError_,
                                                    //double forceYError_,
                                                    float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,                                                    
                                                    bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                                    uint16_t maxPixelsRetrieved,
                                                    Queue& queue);
  }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
