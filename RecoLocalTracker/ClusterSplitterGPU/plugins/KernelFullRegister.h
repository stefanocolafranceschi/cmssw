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

    template <typename TrackerTraits, uint32_t maxPixels, uint8_t maxSubClusters>
    struct JetSplitFullRegister {

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
                                      uint16_t* workOnMe, 
                                      uint32_t numClustersToRun) const {

            // Get thread and grid indices
            //const auto threadIdx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]; // Thread index within the block
            //const auto blockIdx  = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];   // Block index
            //const auto blockDim  = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]; // Threads per block
            //const auto gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

            const auto globalThreadIdx = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

            //constexpr uint8_t maxSubClusters = maxSubClusters_small;
            //constexpr uint16_t extendedMaxPixels = extendedMaxPixels_small;

            if (globalThreadIdx >= numClustersToRun) return;
            const uint32_t clusterIdx = workOnMe[globalThreadIdx];

            //printf("KernelRegister; Running on cluster=%u \n", clusterIdx );

            uint16_t pixelX_cache[maxPixels];
            uint16_t pixelY_cache[maxPixels];
            uint16_t pixelADC_cache[maxPixels];
            uint16_t pixel_info[maxPixels];

            uint8_t scoresIndices[maxPixels];
            float scoresValues[maxPixels];

            float clx[maxSubClusters];
            float cly[maxSubClusters];
            float cls[maxSubClusters];
            float oldclx[maxSubClusters];
            float oldcly[maxSubClusters];
            float clx_temp[maxSubClusters];
            float cly_temp[maxSubClusters];
            float cls_temp[maxSubClusters];
            //uint8_t clusterForPixel[extendedMaxPixels]; 

            bool blockShouldStop;


            const uint32_t begin = geoclusterView.pixelStart(clusterIdx);
            const uint32_t end = begin + geoclusterView.pixelCount(clusterIdx);                
            const uint32_t pixelCounter = end - begin;
            uint32_t ClusterCharge = geoclusterView.ClusterCharge(clusterIdx);

            //if (pixelCounter > maxPixels) return;

            // NO ENOUGH SPACE ON THE TEMPORARY ARRAY
            //if (static_cast<int>(begin) < 0 || static_cast<int>(end) < 0 || static_cast<int>(end) > static_cast<int>(digiView.metadata().size())) {
            //    // Avoid crash if the end is kinda wrong/overflown
            //    return;
            //}

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

                        // Aligning to the original "fittingSplit" variables..
                        uint16_t sizeX = expSizeX;
                        uint16_t sizeY = expSizeY;
                        uint8_t meanExp = std::floor( ClusterCharge / expectedADC + 0.5f);


                        if (meanExp <= 1) {
                            ///if (verbose_) printf("meanExp <= 1 writing cluster");
                            storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
                        }
                        else {

                            // Loading the data from the SoA --------------------------
                            uint32_t offset = 0;
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
                            //-----------------------------------



                            // Splitting the pixels and writing them for the current clusterIdx
                            ///if (verbose_) printf("cluster has meanExp=%d\n", meanExp);

                            // Compute the initial values, set all distances and centers to -999
                            ///if (verbose_) printf("Computing initial values, set all distances");


                            // Resetting temp variables
                            for (uint16_t j = 0; j < meanExp; j++) {
                                oldclx[j] = -999;
                                oldcly[j] = -999;
                                clx[j] = pixelX_cache[0] + j;
                                cly[j] = pixelY_cache[0] + j;
                                cls[j] = 0;
                            }


                            // Main recursive algorithm ----------------
                            blockShouldStop = false;
                            for (uint8_t remainingSteps = 100; remainingSteps > 0 && !blockShouldStop; --remainingSteps) {


                                //if (verbose_) printf("---------------\n");
                                //if (verbose_) printf("REMAINING STEPS : %d\n", remainingSteps);

                                for (uint16_t pixelIdx = 0; pixelIdx < maxPixels; pixelIdx++) {
                                        scoresIndices[pixelIdx] = 0;
                                        scoresValues[pixelIdx] = 0;
                                }

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
                            

                                // SORT SCORES ------------------------------------
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


                                // Save current cluster centers before resetting
                                for (uint16_t subcluster_index = 0; subcluster_index < maxSubClusters; subcluster_index++) {
                                    clx_temp[subcluster_index] = clx[subcluster_index];
                                    cly_temp[subcluster_index] = cly[subcluster_index];  
                                    cls_temp[subcluster_index] = cls[subcluster_index];
                                }

                                // Iterating over Scores Indices and Values
                                for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
                                    uint16_t pixel_index = scoresIndices[i];
                                    int subpixel_counter=-1;

                                    for (uint16_t i = 0; i < pixelCounter && i < maxPixels; ++i) {
  
                                        // Subpixel "simulation"
                                        uint8_t sub = pixel_info[i];
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
                                                //printf("pixelCounter=%u subcluster_index=%u clx[%u]=%f cly[%u]=%f cls[%u]=%f sizeX=%u sizeY=%u\n ",i,subcluster_index,subcluster_index,clx[subcluster_index],subcluster_index,cly[subcluster_index],subcluster_index,cls[subcluster_index],sizeX,sizeY);

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
                                            //clusterForPixel[subpixel_counter] = cl;

                                            //printf("DEBUG, i=%u subpixel_counter=%u cl=%u\n",i, subpixel_counter, cl);                                            
                                            //printf("remainingSteps=%u GPU: pixel_index=%u k=%u -> subpixel=%u -> cl=%d "
                                            //        "charge=%u est=%.4f (maxEst=%.4f) pixel=(%.2f,%.2f) -> cl center=(%.2f,%.2f)\n",
                                            //        remainingSteps, pixel_index, k, subpixel_counter, cl,
                                            //        charge, maxEst, maxEst, temp_x, temp_y, clx[cl], cly[cl]);

                                            //printf("remainingSteps=%u SER: pixel_index=%u k=%u -> subpixel=%u -> cl=%d charge=%u\n",
                                            //    remainingSteps, pixel_index, k, subpixelOffset[pixel_index] + k, cl, charge);
                                        }
                                    }
                                }



                                // Recompute cluster centers
                                //if (verbose_) printf("Recomputing cluster centers.........\n");


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
                                //for (uint16_t oo = 0; oo < meanExp; oo++) {
                                //    printf("remainingSteps=%u oldclx[%u]=%f oldcly[%u]=%f clx[%u]=%f cly[%u]=%f\n",remainingSteps, oo,oldclx[oo],oo,oldcly[oo],oo,clx[oo],oo,cly[oo]);
                                //}




                                // Updating clx, cly and clx
                                // but recalculating on the fly the clusterForPixel

                                for (uint16_t ii = 0; ii < pixelCounter && ii < maxPixels; ++ii) {
                                    uint16_t pixel_index = scoresIndices[ii];

                                    int subpixel_counter=-1;

                                    for (uint16_t i = 0; i < pixelCounter; ++i) {

                                        // Subpixel "simulation"
                                        uint8_t sub = pixel_info[i];
                                        uint16_t adc = pixelADC_cache[i];
                                        uint16_t perDiv = adc / sub;


                                        for (uint16_t k = 0; k < sub; ++k) {

                                            if (i > static_cast<uint32_t>(pixel_index)) break;
                                            if (i != static_cast<uint32_t>(pixel_index)) continue;

                                            float maxEst = 0.f;
                                            int cl = -1;

                                            float temp_x= pixelX_cache[i];
                                            float temp_y= pixelY_cache[i];

                                            uint16_t x = pixelX_cache[i];
                                            uint16_t y = pixelY_cache[i];

                                            for (uint16_t subcluster_index = 0; subcluster_index < meanExp && subcluster_index < maxSubClusters; ++subcluster_index) {
                                                //printf("pixelCounter=%u subcluster_index=%u clx[%u]=%f cly[%u]=%f cls[%u]=%f sizeX=%u sizeY=%u\n ",i,subcluster_index,subcluster_index,clx[subcluster_index],subcluster_index,cly[subcluster_index],subcluster_index,cls[subcluster_index],sizeX,sizeY);

                                                float cx = clx_temp[subcluster_index];
                                                float cy = cly_temp[subcluster_index];
                                                float clusterSignal = cls_temp[subcluster_index];

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
                                            clx[ cl ] += pixelX_cache[i] * charge;
                                            cly[ cl ] += pixelY_cache[i] * charge;
                                            cls_temp[ cl ] += charge;
                                            cls[ cl ] += charge;

                                            
                                 //           uint8_t subj = pixel_info[iii];
                                 //           uint16_t adcj = pixelADC_cache[iii];
                                //            if ( kkk == subj - 1 ) {
                                //                kkk=0;
                                 //               iii++;
                                 //           }
                                            
                             //               uint16_t perDivj = adcj / subj;
                               //             uint16_t chargej = (kkk == subj - 1) ? adcj - perDivj * kkk : perDivj;
                                           

                                //            if ( iii < pixelCounter) {
                               //             clx[ cl ] += pixelX_cache[iii] * chargej;
                              //              cly[ cl ] += pixelY_cache[iii] * chargej;
                              //              cls[ cl ] += chargej;
                            //                kkk++;

                                            //printf("DUMPnew clusterForPixel[pixel_index]=%u i=%u  x=%u  y=%u  charge=%u\n", cl, i, pixelX_cache[i], pixelY_cache[i], charge);
                                            //}


                                        }
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


/*
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

        
                                        //if (kkk > (end - begin)) {
                                        //    printf("Error, more clusters than original");
                                        //    return;
                                        //}

                                        //outputDigis.clus(outIdx)      = clusterIndex;
                                        //outputDigis.xx(outIdx)        = x;
                                        //outputDigis.yy(outIdx)        = y;
                                        //outputDigis.adc(outIdx)       = writeCharge;
                                        //outputDigis.rawIdArr(outIdx)  = rawIdArr;
                                        //outputDigis.moduleId(outIdx)  = moduleId;

                                        alpaka::atomicAdd(acc, pixelCounterDevice, 1u);
        
                                        //if (verbose_) {
                                            uint16_t moduleId = geoclusterView.moduleId(clusterIdx);
                                            printf("candIdx=%u/%u moduleId=%u NSplit cl=%d rawIdArr %d pixel_X[%d]=%u pixel_Y[%d]=%u ADC=%d \n",
                                               candIdx,numCandidates, moduleId, cl, rawIdArr, i, x, i, y, writeCharge);
                                        //}
                                        pixelOffset++;
                                        kkk++;
                                    }
                                }
                            }
*/


                            //storeOutputDigis
                            uint32_t kkk = 0;
                            for (uint16_t mycl = 0; mycl < static_cast<int>(meanExp); ++mycl) {

                                uint32_t clusterIndex = alpaka::atomicAdd(acc, clusterCounterDevice, 1u);
                                uint32_t pixelStartWritingAt = alpaka::atomicAdd(acc, pixelCounterDevice, 0u);
                                uint32_t pixelOffset = 0;



                                // Iterating over Scores Indices and Values
                                for (uint16_t ii = 0; ii < pixelCounter && ii < maxPixels; ++ii) {
                                    uint16_t pixel_index = scoresIndices[ii];

                                    for (uint16_t i = 0, j=0; i < pixelCounter && j < pixelCounter; ++i,j++) {

                                        // Subpixel "simulation"
                                        uint8_t sub = pixel_info[i];
                                        uint8_t subj = pixel_info[j];

                                        uint16_t adc = pixelADC_cache[i];
                                        uint16_t adcj = pixelADC_cache[j];

                                        uint16_t perDiv = adc / sub;
                                        uint16_t perDivj = adcj / subj;

                                        uint32_t writeCharge = 0;
                                        int cl = -1;

                                        for (uint16_t k = 0; k < sub; ++k) {

                                            uint32_t charge = (k == sub - 1) ? adc - perDiv * k : perDiv;

                                            if (i > static_cast<uint32_t>(pixel_index)) break;
                                            if (i != static_cast<uint32_t>(pixel_index)) continue;

                                            float maxEst = 0.f;

                                            float temp_x = pixelX_cache[i];
                                            float temp_y = pixelY_cache[i];

                                            for (uint16_t subcluster_index = 0; subcluster_index < meanExp && subcluster_index < maxSubClusters; ++subcluster_index) {
                                                //printf("pixelCounter=%u subcluster_index=%u clx[%u]=%f cly[%u]=%f cls[%u]=%f sizeX=%u sizeY=%u\n ",i,subcluster_index,subcluster_index,clx[subcluster_index],subcluster_index,cly[subcluster_index],subcluster_index,cls[subcluster_index],sizeX,sizeY);

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
                                        }
                                        for (uint16_t k = 0; k < subj; ++k) {
                                            uint32_t chargej = (k == subj - 1) ? adcj - perDivj * k : perDivj;
                                            if (cl==mycl) writeCharge += chargej;
                                        }

                                        if (writeCharge>0) {

                                            //uint32_t outIdx = pixelStartWritingAt + pixelOffset;
                                            //uint32_t rawIdArr = digiView.rawIdArr(begin);

                                            //if (kkk > (end - begin)) {
                                            //    printf("Error, more clusters than original");
                                            //    return;
                                            //}

                                            //outputDigis.clus(outIdx)      = clusterIndex;
                                            //outputDigis.xx(outIdx)        = x;
                                            //outputDigis.yy(outIdx)        = y;
                                            //outputDigis.adc(outIdx)       = writeCharge;
                                            //outputDigis.rawIdArr(outIdx)  = rawIdArr;
                                            //outputDigis.moduleId(outIdx)  = moduleId;

                                            alpaka::atomicAdd(acc, pixelCounterDevice, 1u);

                                            if (verbose_) {
                                                        uint16_t moduleId = geoclusterView.moduleId(clusterIdx);
                                                        printf("candIdx=%u/%u moduleId=%u NSplit cl=%d  pixel_X[%d]=%u pixel_Y[%d]=%u ADC=%d \n",
                                                            candIdx,numCandidates, moduleId, mycl, j, pixelX_cache[j], j, pixelY_cache[j], writeCharge);
                                            }
                                        }
                                    }
                                }
                            }
                            // done storing
                        }

                    }
                    alreadySplit = true;
                }
            }


            if (!doSplit && !alreadySplit) {                       
                storeOutputDigis(acc, digiView, outputDigis, begin, end, clusterCounterDevice, pixelCounterDevice);
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
