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

#include <cmath>
#include <alpaka/math/MathStdLib.hpp>

using namespace alpaka;
using namespace reco;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;
  namespace Splitting {

    template <typename TrackerTraits, uint32_t maxPixels, uint8_t maxSubClusters, uint16_t extendedMaxPixels>
    struct JetSplitShared {

        // Main operator function
        template <typename TAcc, typename = std::enable_if_t<isAccelerator<TAcc>>>
        ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                      //TrackingRecHitSoAConstView<TrackerTraits> hitView,
                                      SiPixelDigisSoAView digiView,
                                      //SiPixelClustersSoAConstView clusterView,
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
                                      //SiPixelClustersSoAView outputClusters,
                                      //clusterProperties* clusterPropertiesDevice,
                                      uint32_t* clusterCounterDevice,
                                      uint32_t* pixelCounterDevice,                                      
                                      //double forceXError_, double forceYError_,
                                      float vertexX, float vertexY, float vertexZ, float vertexEta, float vertexPhi,
                                      bool verbose_, bool debugMode, int targetDetId, uint16_t targetClusterOffset,
                                      uint16_t* workOnMe, uint16_t numClustersToRun) const {

            // Get thread and grid indices
            const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]; // Thread index within the block
            const auto blockIdx  = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];   // Block index
            const auto blockDim  = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]; // Threads per block
            //const auto gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

            //constexpr uint8_t maxSubClusters = maxSubClusters_large;
            //constexpr uint16_t extendedMaxPixels = extendedMaxPixels_large;

            // Run cooperatively one block with several threads
            // handles one cluster

            if (blockIdx >= numClustersToRun) return;
            const uint32_t clusterIdx = workOnMe[blockIdx];

//printf("KernelShared; maxPixels=%u Running on blockIdx=%u threadIdx=%u clusterIdx=%u\n", maxPixels, blockIdx, threadIdx, clusterIdx);

            __attribute__((shared)) uint16_t pixelX_cache[maxPixels];
            __attribute__((shared)) uint16_t pixelY_cache[maxPixels];
            __attribute__((shared)) uint16_t pixelADC_cache[maxPixels];
            __attribute__((shared)) uint16_t pixel_info[maxPixels];
            __attribute__((shared)) uint8_t subpixelOffset[maxPixels];

            __attribute__((shared)) uint8_t scoresIndices[maxPixels];
            __attribute__((shared)) float scoresValues[maxPixels];

            __attribute__((shared)) float clx[maxSubClusters];
            __attribute__((shared)) float cly[maxSubClusters];
            __attribute__((shared)) float cls[maxSubClusters];
            __attribute__((shared)) float oldclx[maxSubClusters];
            __attribute__((shared)) float oldcly[maxSubClusters];

            __attribute__((shared)) uint8_t clusterForPixel[extendedMaxPixels]; 

            __attribute__((shared)) bool blockShouldStop;

            // Temporary variables used for accumulators
            __attribute__((shared)) float temp_clx[maxPixels][maxSubClusters];
            __attribute__((shared)) float temp_cly[maxPixels][maxSubClusters];
            __attribute__((shared)) float temp_cls[maxPixels][maxSubClusters];

//            __attribute__((shared)) float shared_cls[maxSubClusters];
            


            const uint32_t begin = geoclusterView.pixelStart(clusterIdx);
            const uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);                
            const uint32_t pixelCounter = end - begin;
            uint32_t ClusterCharge = geoclusterView.ClusterCharge(clusterIdx);

            if (pixelCounter > maxPixels) return;

            // NO ENOUGH SPACE ON THE TEMPORARY ARRAY
            if (static_cast<int>(begin) < 0 || static_cast<int>(end) < 0 || static_cast<int>(end) > static_cast<int>(digiView.metadata().size())) {
                // Avoid crash if the end is kinda wrong/overflown
                return;
            }


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

            const uint32_t numCandidates = static_cast<uint32_t>(candidateView.metadata().size());
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
                float jetPt = sqrtf(jetPx * jetPx + jetPy * jetPy);
                float jetP  = sqrtf(jetPx * jetPx + jetPy * jetPy + jetPz * jetPz);
                float jetEta = 0.5 * logf((jetP + jetPz) / (jetP - jetPz));
                float jetPhi = atan2f(jetPy, jetPx);

                // Print the jet information 
                ///if (verbose_) printf("In globalThreadId=%u, Jet Information:\n", globalThreadId);
                ///if (verbose_) printf("  jetPx = %.3f, jetPy = %.3f, jetPz = %.3f\n", jetPx, jetPy, jetPz);
                ///if (verbose_) printf("  jetPt = %.3f, jetEta = %.3f, jetPhi = %.3f\n\n", jetPt, jetEta, jetPhi);

                // Compute the cluster's relative eta and phi
                float r = sqrtf(relX * relX + relY * relY + relZ * relZ);
                float clusterEta = 0.5 * logf((r + relZ) / (r - relZ));  // Pseudorapidity formula
                float clusterPhi = atan2f(relY, relX);  // Azimuthal angle

                // Compute differences and deltaR (assuming 'jetEta' and 'jetPhi' are known)
                float deltaEta = clusterEta - jetEta;
                float deltaPhi = atan2f(sinf(clusterPhi - jetPhi), cosf(clusterPhi - jetPhi));  // Adjust for periodicity
                float deltaR = sqrtf(deltaEta * deltaEta + deltaPhi * deltaPhi);

                doSplit = deltaR < deltaR_;

                // Check deltaR condition and split clusters if applicable
                if (doSplit) {
                    //printf("This clusterOffset: %u has deltaR < deltaR_ and it might be split\n",clusterIdx);

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
                        float jetZOverRho = sqrtf(jetTanAlpha * jetTanAlpha + jetTanBeta * jetTanBeta);

                        expSizeX = expSizeXAtLorentzAngleIncidence_ +
                                         fabsf(expSizeXDeltaPerTanAlpha_ * (jetTanAlpha - tanLorentzAngles));

                        expSizeY = sqrtf((expSizeYAtNormalIncidence_ * expSizeYAtNormalIncidence_) +
                                                   thickness * thickness / (pitchY * pitchY) * jetTanBeta * jetTanBeta);
                        
                        if (expSizeX < 1.f) expSizeX = 1.f;
                        if (expSizeY < 1.f) expSizeY = 1.f;

                        expectedADC = sqrtf(1.08f + jetZOverRho * jetZOverRho) * centralMIPCharge_;
         
                        //printf("Trying to split: charge=%d expSizeX=%f expSizeY=%f\n",
                        //        static_cast<int>(ClusterCharge), expSizeX, expSizeY);

                        if ( ClusterCharge > expectedADC * chargeFracMin_ &&
                               ( geoclusterView.sizeX(clusterIdx) > expSizeX + 1 || geoclusterView.sizeY(clusterIdx) > expSizeY + 1)) {
                            split = true;
                        }
                    }

                    if (split) {

                        // Aligning to the original "fittingSplit" variables..
                        uint16_t sizeX = expSizeX;
                        uint16_t sizeY = expSizeY;
                        uint8_t meanExp = std::floor( ClusterCharge / expectedADC + 0.5f);
//printf("CHECK meanExp=%u ClusterCharge=%u\n", meanExp,ClusterCharge);

                        if (meanExp <= 1) {
                            ///if (verbose_) printf("meanExp <= 1 writing cluster");
                            storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                        }
                        else {
                            // Loading the data from the SoA --------------------------

                            // Filling local cache for faster access and sub (for avoiding repeated pixels large arrays)
                            for (uint16_t idx = threadIdx; idx < pixelCounter && idx < maxPixels; idx += blockDim) {
                                uint16_t jj = begin + idx;
                                if (jj>end) break;
                                pixelX_cache[idx] = digiView.xx(jj);
                                pixelY_cache[idx] = digiView.yy(jj);
                                uint16_t charge = digiView.adc(jj);
                                pixelADC_cache[idx] = charge;

                                uint16_t sub = static_cast<int>( charge) / chargePerUnit_ * expectedADC / centralMIPCharge_;
                                
                                if (sub < 1) sub = 1;
                                pixel_info[idx] = sub;
                            }

                            if (threadIdx == 0) {
                                uint32_t offset = 0;
                                for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
                                    subpixelOffset[i] = offset;
                                    offset += pixel_info[i];  // accumulate subpixels
                                }
                            }
                            alpaka::syncBlockThreads(acc);                            
                            //-----------------------------------



                            // Splitting the pixels and writing them for the current clusterIdx
                            ///if (verbose_) printf("cluster has meanExp=%d\n", meanExp);

                            // Compute the initial values, set all distances and centers to -999
                            ///if (verbose_) printf("Computing initial values, set all distances");


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
                                //cls[j]    = 0;
                            }

                            alpaka::syncBlockThreads(acc);  // Again ensure sync after re-initialization


                            // Main recursive algorithm ----------------
                            blockShouldStop = false;
                            for (uint8_t remainingSteps = 100; remainingSteps > 0 && !blockShouldStop; --remainingSteps) {


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
                                            if ( fabsf(distanceX) > sizeX / 2.f) {
                                                float diff = fabsf(distanceX) - sizeX / 2.f + 1.f;
                                                distX = diff * diff;
                                            } else {
                                                float scaled = 2.f * distanceX / sizeX;
                                                distX = scaled * scaled;
                                            }

                                            float distY = 0.f;
                                            if ( fabsf(distanceY) > sizeY / 2.f) {
                                                float diff = fabsf(distanceY) - sizeY / 2.f + 1.f;
                                                distY = diff * diff;
                                            } else {
                                                float scaled = 2.f * distanceY / sizeY;
                                                distY = scaled * scaled;
                                            }

                                            float dist = sqrtf(distX + distY);
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


                                // Reseting temporary variables
                                for (uint16_t i = threadIdx; i < maxPixels; i += blockDim) {
                                    scoresIndices[i] = 0;
                                    scoresValues[i] = 0.f;
                                }
                                alpaka::syncBlockThreads(acc);

                                float sizeX_half = sizeX / 2.f;
                                float sizeY_half = sizeY / 2.f;
                                float invsizeX = 1.f / sizeX;
                                float invsizeY = 1.f / sizeY;

                                for (uint16_t pixelIdx = threadIdx; pixelIdx < pixelCounter; pixelIdx += blockDim) {
                                    if (pixelIdx < maxPixels) {
                                        float minDist = std::numeric_limits<float>::max();
                                        float secondMinDist = std::numeric_limits<float>::max();

                                        float temp_originalpixels_x = pixelX_cache[pixelIdx];
                                        float temp_originalpixels_y = pixelY_cache[pixelIdx];

                                        // If not found, skip this pixelIdx
                                        if (temp_originalpixels_x == -1.f || temp_originalpixels_y == -1.f) continue;

                                        // Compute distances to all subclusters
                                        for (uint8_t subClusterIdx = 0; subClusterIdx < meanExp; ++subClusterIdx) {
                                            float distanceX = temp_originalpixels_x - clx[subClusterIdx];
                                            float distanceY = temp_originalpixels_y - cly[subClusterIdx];

                                            float distX = 0.f;
                                            if ( fabsf(distanceX) > sizeX_half ) {
                                                float diff = fabsf(distanceX) - sizeX_half + 1.f;
                                                distX = diff * diff;
                                            } else {
                                                float scaled = 2.f * distanceX * invsizeX;
                                                distX = scaled * scaled;
                                            }

                                            float distY = 0.f;
                                            if ( fabsf(distanceY) > sizeY_half) {
                                                float diff = fabsf(distanceY) - sizeY_half + 1.f;
                                                distY = diff * diff;
                                            } else {
                                                float scaled = 2.f * distanceY * invsizeY;
                                                distY = scaled * scaled;
                                            }

                                            float dist = sqrtf(distX + distY);

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


/*
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
*/
// Optimized odd-even sort for small arrays
for (uint16_t phase = 0; phase < pixelCounter; phase++) {
    bool isEvenPhase = (phase % 2 == 0);
    
    // Calculate starting index for this thread
    uint16_t startIdx = isEvenPhase ? threadIdx * 2 : threadIdx * 2 + 1;
    
    if (startIdx < pixelCounter - 1) {
        if (scoresValues[startIdx] > scoresValues[startIdx + 1] ||
            (scoresValues[startIdx] == scoresValues[startIdx + 1] && 
             scoresIndices[startIdx] > scoresIndices[startIdx + 1])) {
            std::swap(scoresValues[startIdx], scoresValues[startIdx + 1]);
            std::swap(scoresIndices[startIdx], scoresIndices[startIdx + 1]);
        }
    }
    alpaka::syncBlockThreads(acc);
}


/*
                                // Odd-even sort simpler than bitonic, better than bubble sort
                                for (uint16_t phase = 0; phase < pixelCounter; phase++) {
                                    if (phase % 2 == 0) {
                                        // Even phase: compare (0,1), (2,3), (4,5), ...
                                        for (uint16_t idx = threadIdx * 2; idx < pixelCounter - 1; idx += blockDim * 2) {
                                            if (scoresValues[idx] > scoresValues[idx + 1] ||
                                                (scoresValues[idx] == scoresValues[idx + 1] && scoresIndices[idx] > scoresIndices[idx + 1])) {
                                                std::swap(scoresValues[idx], scoresValues[idx + 1]);
                                                std::swap(scoresIndices[idx], scoresIndices[idx + 1]);
                                            }
                                        }
                                    } else {
                                        // Odd phase: compare (1,2), (3,4), (5,6), ...
                                        for (uint16_t idx = threadIdx * 2 + 1; idx < pixelCounter - 1; idx += blockDim * 2) {
                                            if (scoresValues[idx] > scoresValues[idx + 1] ||
                                                (scoresValues[idx] == scoresValues[idx + 1] && scoresIndices[idx] > scoresIndices[idx + 1])) {
                                                std::swap(scoresValues[idx], scoresValues[idx + 1]);
                                                std::swap(scoresIndices[idx], scoresIndices[idx + 1]);
                                            }
                                        }
                                    }
                                    alpaka::syncBlockThreads(acc);
                                }
*/
                                //if (threadIdx == 0 && verbose_) {
                                //    printf("Cluster %u Scores:\n", clusterIdx);
                                //    for (uint16_t k = 0; k < pixelCounter; k++) {
                                //        printf("Cluster %u Score = %.5f, Index = %d\n", clusterIdx, scoresValues[k], scoresIndices[k]);
                                //    }
                                //}
/*
if (threadIdx==0) {

                                //float localCls[maxSubClusters] = {0.f};
                                // Iterating over Scores Indices and Values                                    
                                // Each thread handles one score index
                                //for (uint16_t i = 0; i < pixelCounter && i < maxPixels; i++) {
                                for (uint16_t i = threadIdx; i < pixelCounter && i < maxPixels; i += blockDim) {

                                    uint16_t pixel_index = scoresIndices[i];

                                    uint8_t sub = pixel_info[pixel_index];
                                    uint16_t adc = pixelADC_cache[pixel_index];
                                    uint16_t perDiv = adc / sub;

                                    float temp_x = pixelX_cache[pixel_index];
                                    float temp_y = pixelY_cache[pixel_index];

                                    float scaledFrac = expectedADC * fractionalWidth_;
                                    float invScaledFrac = 1.f / scaledFrac;
                                    float sizeX_half = sizeX / 2.f;
                                    float sizeY_half = sizeY / 2.f;

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
                                            //float absX = std::abs(dx), absY = std::abs(dy);
                                            float absX = fabsf(dx), absY = fabsf(dy);
                                            if ( absX > sizeX_half ) {
                                                //dist += (absX - sizeX / 2.f + 1.f) * (absX - sizeX / 2.f + 1.f);
                                                float norm_dx = absX - sizeX / 2.f + 1.f;
                                                dist = fmaf(norm_dx, norm_dx, dist);
                                            }

                                            else {
                                                //dist += (2.f * dx / sizeX) * (2.f * dx / sizeX);
                                                float norm_dx = 2.f * dx / sizeX;
                                                dist = fmaf(norm_dx, norm_dx, dist);
                                            }

                                            if (absY > sizeY_half) {
                                                //dist += (absY - sizeY / 2.f + 1.f) * (absY - sizeY / 2.f + 1.f);
                                                float norm_dy = absY - sizeY / 2.f + 1.f;
                                                dist = fmaf(norm_dy, norm_dy, dist);
                                            }

                                            else {
                                                //dist += (2.f * dy / sizeY) * (2.f * dy / sizeY);
                                                float norm_dy = 2.f * dy / sizeY;
                                                dist = fmaf(norm_dy, norm_dy, dist);
                                            }

                                            float distance = sqrtf(dist);
                                            float nsig = (clusterSignal - expectedADC) * invScaledFrac;
                                            float clQest = 1.f / (1.f + expf(nsig)) + 1e-6f;
                                            float clDest = 1.f / (distance + 0.05f);

                                            //float inv_sqrt_dist = alpaka::math::rsqrt(dist);
                                            //float clDest = inv_sqrt_dist / (1.f + 0.05f * inv_sqrt_dist);

                                            float est = clQest * clDest;

                                            if (est > maxEst) {
                                                cl = subcluster_index;
                                                maxEst = est;
                                            }
                                        }

                                        uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;

                                        // Atomically add charge to cluster total
                                        //alpaka::atomicAdd(acc, &cls[cl], static_cast<float>(charge));

                                        cls[cl] += static_cast<float>(charge);

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
}
*/




                                //Hierarchical Processing
                                 __attribute__((shared)) float shared_estimates[maxPixels];
                                 __attribute__((shared)) int shared_best_clusters[maxPixels];

                                // Process pixels sequentially to maintain exact dependency chain
                                for (uint16_t i = 0; i < pixelCounter && i < maxPixels; i++) {
                                    
                                    uint16_t pixel_index = scoresIndices[i];
                                    
                                    uint8_t sub = pixel_info[pixel_index];
                                    uint16_t adc = pixelADC_cache[pixel_index];
                                    uint16_t perDiv = adc / sub;
                                    
                                    float temp_x = pixelX_cache[pixel_index];
                                    float temp_y = pixelY_cache[pixel_index];
                                    
                                    float scaledFrac = expectedADC * fractionalWidth_;
                                    float invScaledFrac = 1.f / scaledFrac;
                                    float sizeX_half = sizeX / 2.f;
                                    float sizeY_half = sizeY / 2.f;
                                    
                                    // Process each subpixel k sequentially (to maintain exact semantics)
                                    for (uint8_t k = 0; k < sub; ++k) {
                                        
                                        float maxEst = 0.f;
                                        int cl = -1;
                                        
                                        // PARALLELIZE THIS INNER LOOP: Each thread handles different subclusters
                                        for (uint16_t sc_start = 0; sc_start < meanExp && sc_start < maxSubClusters; sc_start += blockDim) {
                                            
                                            uint16_t subcluster_index = sc_start + threadIdx;
                                            float est = 0.f;
                                            
                                            if (subcluster_index < meanExp && subcluster_index < maxSubClusters) {
                                                
                                                float cx = clx[subcluster_index];
                                                float cy = cly[subcluster_index];
                                                float clusterSignal = cls[subcluster_index]; // Use current cls values
                                                
                                                float dx = temp_x - cx;
                                                float dy = temp_y - cy;
                                                
                                                float dist = 0.f;
                                                float absX = fabsf(dx), absY = fabsf(dy);
                                                if (absX > sizeX_half) {
                                                    float norm_dx = absX - sizeX / 2.f + 1.f;
                                                    dist = fmaf(norm_dx, norm_dx, dist);
                                                }
                                                else {
                                                    float norm_dx = 2.f * dx / sizeX;
                                                    dist = fmaf(norm_dx, norm_dx, dist);
                                                }
                                                
                                                if (absY > sizeY_half) {
                                                    float norm_dy = absY - sizeY / 2.f + 1.f;
                                                    dist = fmaf(norm_dy, norm_dy, dist);
                                                }
                                                else {
                                                    float norm_dy = 2.f * dy / sizeY;
                                                    dist = fmaf(norm_dy, norm_dy, dist);
                                                }
                                                
                                                float distance = sqrtf(dist);
                                                float nsig = (clusterSignal - expectedADC) * invScaledFrac;
                                                float clQest = 1.f / (1.f + expf(nsig)) + 1e-6f;
                                                float clDest = 1.f / (distance + 0.05f);
                                                
                                                est = clQest * clDest;
                                            }
                                            
                                            // Store results in shared memory for reduction
                                            shared_estimates[threadIdx] = est;
                                            shared_best_clusters[threadIdx] = subcluster_index;
                                            alpaka::syncBlockThreads(acc);
                                            
                                            // Reduction to find maximum estimate within this batch
                                            for (uint16_t stride = blockDim / 2; stride > 0; stride /= 2) {
                                                if (threadIdx < stride) {
                                                    if (shared_estimates[threadIdx + stride] > shared_estimates[threadIdx]) {
                                                        shared_estimates[threadIdx] = shared_estimates[threadIdx + stride];
                                                        shared_best_clusters[threadIdx] = shared_best_clusters[threadIdx + stride];
                                                    }
                                                }
                                                alpaka::syncBlockThreads(acc);
                                            }
                                            
                                            // Thread 0 updates the global maximum
                                            if (threadIdx == 0 && shared_estimates[0] > maxEst) {
                                                maxEst = shared_estimates[0];
                                                cl = shared_best_clusters[0];
                                            }
                                            alpaka::syncBlockThreads(acc);
                                        }
                                        
                                        // Thread 0 assigns the charge and writes cluster assignment
                                        if (threadIdx == 0) {
                                            uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;
                                            cls[cl] += static_cast<float>(charge);
                                            clusterForPixel[subpixelOffset[pixel_index] + k] = cl;
                                        }
                                        alpaka::syncBlockThreads(acc);
                                    }
                                }

                                alpaka::syncBlockThreads(acc);





                                // Recompute cluster centers
                                if (verbose_) printf("Recomputing cluster centers.........\n");
                                if (threadIdx == 0) blockShouldStop = true;  // Assume we will stop
                                alpaka::syncBlockThreads(acc);


                                //if (threadIdx == 0) {
                                    //for (uint16_t oo = 0; oo < meanExp; oo++) {
                                    //    printf("remainingSteps=%u oldclx[%u]=%f oldcly[%u]=%f clx[%u]=%f cly[%u]=%f\n",remainingSteps, oo,oldclx[oo],oo,oldcly[oo],oo,clx[oo],oo,oldcly[oo]);
                                    //}
                                //}

                                // Each thread checks part of the subcluster array
                                for (uint16_t i = threadIdx; i < meanExp; i += blockDim) {
                                    if (fabsf(clx[i] - oldclx[i]) > 0.01f || fabsf(cly[i] - oldcly[i]) > 0.01f) {
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

                                //for (uint16_t i = 0; i < extendedMaxPixels; ++i) {
                                //    printf("i=%d clusterForPixel=%d\n",i,clusterForPixel[i]);
                                //}


                                // Initialize thread-local accumulators
                                for (uint16_t c = 0; c < maxSubClusters; c++) {
                                    temp_clx[threadIdx][c] = 0.0f;
                                    temp_cly[threadIdx][c] = 0.0f;
                                    temp_cls[threadIdx][c] = 0.0f;
                                }
                                alpaka::syncBlockThreads(acc);

                                for (uint16_t i = threadIdx; i < pixelCounter && i < maxPixels; i += blockDim) {
                                    uint16_t x = pixelX_cache[i];
                                    uint16_t y = pixelY_cache[i];
                                    uint8_t sub = pixel_info[i];
                                    uint16_t adc = pixelADC_cache[i];
                                    uint16_t perDiv = adc / sub;
                                    
                                    uint32_t baseIdx = subpixelOffset[i];
                                    
                                    for (uint8_t k = 0; k < sub; k++) {
                                        uint16_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;
                                        uint32_t subpixelIdx = baseIdx + k;
                                        
                                        uint16_t clusterId = clusterForPixel[subpixelIdx];
                                        temp_clx[threadIdx][clusterId] += x * charge;
                                        temp_cly[threadIdx][clusterId] += y * charge;
                                        temp_cls[threadIdx][clusterId] += charge;
                                    }
                                }
                                alpaka::syncBlockThreads(acc);

                                for (uint16_t c = threadIdx; c < maxSubClusters; c += blockDim) {
                                    float sum_clx = 0.0f, sum_cly = 0.0f, sum_cls = 0.0f;
                                    for (uint16_t t = 0; t < blockDim; t++) {
                                        sum_clx += temp_clx[t][c];
                                        sum_cly += temp_cly[t][c];
                                        sum_cls += temp_cls[t][c];
                                    }
                                    clx[c] += sum_clx;
                                    cly[c] += sum_cly;
                                    cls[c] += sum_cls;
                                }                                
                                alpaka::syncBlockThreads(acc);


                                // Only use the threads we need
                                if (threadIdx < meanExp) {
                                    uint8_t subcluster_index = threadIdx;
                                    if (cls[subcluster_index] != 0) {
                                        clx[subcluster_index] /= cls[subcluster_index];
                                        cly[subcluster_index] /= cls[subcluster_index];
                                    }
                                    //printf("Center for cluster, clx[%u]=%f cly[%u]=%f\n", subcluster_index, clx[subcluster_index], subcluster_index, cly[subcluster_index]);
                                    cls[subcluster_index] = 0;
                                }
                                // No need for sync here since each thread works on independent data

                            }
                            // end of recursive stuff




                            // Parallel over clusters (outer loop) each thread handles one cluster
                            for (uint16_t cl = threadIdx; cl < static_cast<int>(meanExp); cl += blockDim) {
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
                                        outputDigis.clus(outIdx)      = clusterIndex;
                                        outputDigis.xx(outIdx)        = x;
                                        outputDigis.yy(outIdx)        = y;
                                        outputDigis.adc(outIdx)       = writeCharge;
                                        outputDigis.rawIdArr(outIdx)  = rawIdArr;
                                        outputDigis.moduleId(outIdx)  = moduleId;
                            */
                                        alpaka::atomicAdd(acc, pixelCounterDevice, 1u);

                                        //if (verbose_) {
                                            uint16_t moduleId = geoclusterView.moduleId(clusterIdx);            
                                            printf("candIdx=%u/%u moduleId=%u NSplit cl=%d rawIdArr %d pixel_X[%d]=%u pixel_Y[%d]=%u ADC=%d \n",
                                                   candIdx, numCandidates, moduleId, cl, rawIdArr, i, x, i, y, writeCharge);
                                        //}
                                        pixelOffset++;
                                    }
                                }
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


 }  // namespace Splitting
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
