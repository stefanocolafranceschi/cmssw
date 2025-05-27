#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED  // Ensure this file is only used for CUDA

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/stringize.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TFile.h"
#include "TString.h"
#include <vector>
#include <tuple>
#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

//#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsDevice.h"
//#include "DataFormats/TrackingRecHitSoA/interface/TrackingRecHitsHost.h"
//#include "DataFormats/TrackingRecHitSoA/interface/alpaka/TrackingRecHitsSoACollection.h"
//#include "DataFormats/TrackingRecHitSoA/interface/SiPixelHitStatus.h"

#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisDevice.h"
#include "DataFormats/SiPixelDigiSoA/interface/SiPixelDigisHost.h"
#include "DataFormats/SiPixelDigiSoA/interface/alpaka/SiPixelDigisSoACollection.h"

#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersDevice.h"
#include "DataFormats/SiPixelClusterSoA/interface/SiPixelClustersHost.h"
#include "DataFormats/SiPixelClusterSoA/interface/alpaka/SiPixelClustersSoACollection.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/VertexSoA/interface/ZVertexSoA.h"
#include "DataFormats/VertexSoA/interface/ZVertexHost.h"
#include "DataFormats/VertexSoA/interface/ZVertexDevice.h"
#include "DataFormats/VertexSoA/interface/alpaka/ZVertexSoACollection.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "HeterogeneousCore/AlpakaInterface/interface/Backend.h"
#include "HeterogeneousCore/AlpakaInterface/interface/devices.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/SynchronizingEDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDMetadata.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDMetadataSentry.h"
#include "Cluster_test.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometrysSoA.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrysSoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidatesSoA.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidatesSoACollection.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;

class trial : public global::EDProducer<> {
public:
  explicit trial(edm::ParameterSet const& iConfig);
  ~trial() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& setup) const override;

  //uint32_t nHits_;
  //int32_t offset_;
  const double ptMin_;
  TFile* rootFile_;

  double deltaR_;
  double chargeFracMin_;
  float tanLorentzAngle_;
  float tanLorentzAngleBarrelLayer1_;  

  float expSizeXAtLorentzAngleIncidence_;
  float expSizeXDeltaPerTanAlpha_;
  float expSizeYAtNormalIncidence_;
  double centralMIPCharge_;
  double chargePerUnit_;
  double forceXError_;
  double forceYError_;  
  double fractionalWidth_;
  const edm::EDGetTokenT<SiPixelClustersHost> clusterToken_;
  const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelDigisSoACollection> digisToken_;
  //const edm::EDGetTokenT<SiPixelClustersHost> digisToken_;
  //const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelClustersSoACollection> clusterToken_;
  ////const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::TrackingRecHitsSoACollection<pixelTopology::Phase1>> recHitsToken_;
  const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::CandidatesSoACollection> candidateToken_;
  //const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::ZVertexSoACollection> zVertexToken_;
  const device::EDGetToken<ALPAKA_ACCELERATOR_NAMESPACE::ClusterGeometrysSoACollection> geometryToken_;
  bool verbose_;
  bool debugMode;
  int targetDetId;
  uint16_t targetClusterOffset;
  int targetEvent;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  //edm::EDGetTokenT<uint32_t> maxPixelsToken_;
  //edm::EDGetTokenT<std::vector<std::pair<uint32_t, uint32_t>>> clusterPixelCountsToken_;
  edm::EDGetTokenT<std::vector<uint32_t>> clusterPixelCountsToken_;
  const device::EDPutToken<ALPAKA_ACCELERATOR_NAMESPACE::SiPixelDigisSoACollection> outputdigisToken_;
  std::vector<Device> devices_;  
};

trial::trial(edm::ParameterSet const& iConfig)
    : EDProducer(iConfig),
      //nHits_(iConfig.getParameter<uint32_t>("nHits")),
      //offset_(iConfig.getParameter<int32_t>("offset")),
      ptMin_(iConfig.getParameter<double>("ptMin")),
      rootFile_(nullptr),
      deltaR_(iConfig.getParameter<double>("deltaR")),
      chargeFracMin_(iConfig.getParameter<double>("chargeFracMin")),
      tanLorentzAngle_(iConfig.getParameter<double>("tanLorentzAngle")),
      tanLorentzAngleBarrelLayer1_(iConfig.getParameter<double>("tanLorentzAngleBarrelLayer1")),
      expSizeXAtLorentzAngleIncidence_(iConfig.getParameter<double>("expSizeXAtLorentzAngleIncidence")),
      expSizeXDeltaPerTanAlpha_(iConfig.getParameter<double>("expSizeXDeltaPerTanAlpha")),
      expSizeYAtNormalIncidence_(iConfig.getParameter<double>("expSizeYAtNormalIncidence")),
      centralMIPCharge_(iConfig.getParameter<double>("centralMIPCharge")),
      chargePerUnit_(iConfig.getParameter<double>("chargePerUnit")),
      forceXError_(iConfig.getParameter<double>("forceXError")),
      forceYError_(iConfig.getParameter<double>("forceYError")),
      fractionalWidth_(iConfig.getParameter<double>("fractionalWidth")),
      clusterToken_(consumes(iConfig.getParameter<edm::InputTag>("siPixelClusters"))),
      digisToken_(consumes(iConfig.getParameter<edm::InputTag>("siPixelDigis"))),
      ////recHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("trackingRecHits"))),
      candidateToken_(consumes(iConfig.getParameter<edm::InputTag>("candidateInput"))),
      //zVertexToken_(consumes(iConfig.getParameter<edm::InputTag>("zVertex"))),
      geometryToken_(consumes(iConfig.getParameter<edm::InputTag>("geometryInput"))),
      verbose_(iConfig.getParameter<bool>("verbose")),
      debugMode(iConfig.getParameter<bool>("debugMode")),
      targetDetId(iConfig.getParameter<int>("targetDetId")),
      targetClusterOffset(iConfig.getParameter<int>("targetClusterOffset")),
      targetEvent(iConfig.getParameter<int>("targetEvent")),
      vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      //maxPixelsToken_(consumes<uint32_t>(iConfig.getParameter<edm::InputTag>("maxPixels"))),
      clusterPixelCountsToken_(consumes<std::vector<uint32_t>>(iConfig.getParameter<edm::InputTag>("clusterPixelCounts"))),

      outputdigisToken_{produces()}
      {
          devices_ = cms::alpakatools::devices<alpaka::PlatformCudaRt>();
          rootFile_ = new TFile("config_output.root", "RECREATE");
      }


trial::~trial() {
  if (rootFile_) {
    rootFile_->Close();
    delete rootFile_;
  }
}

void trial::produce(edm::StreamID sid, device::Event& deviceEvent, device::EventSetup const& iSetup) const {


    //if ((debugMode) && (deviceEvent.id().event() != static_cast<edm::EventNumber_t>(targetEvent))) {
    //    std::cout << "Skipping this event" << std::endl;
    //    return;
    //}


    // Retrieve the value of maxPixels
    auto const& clusterPixelCounts = deviceEvent.get(clusterPixelCountsToken_);

    //for (size_t i = 0; i < clusterPixelCounts.size(); ++i) {
    //  printf("Pixel count at index %zu: %u\n", i, clusterPixelCounts[i]);
    //}


    //if (deviceEvent.id().event() !=25) return;

    //if (verbose_) std::cout << "Entering in produce method.. testing" << std::endl;

    // Ensure we're selecting the first available GPU device PlatformCudaRt vs PlatformCpu
    auto const& deviceList = cms::alpakatools::devices<alpaka::PlatformCudaRt>();
    if (deviceList.empty()) {
        throw cms::Exception("Configuration") << "No available Alpaka GPU devices found!";
    }

    // Select the first GPU device
    auto const& device = deviceList[0];
    ///if (verbose_) std::cout << "Using GPU device: " << alpaka::getName(device) << std::endl;

    // ---------------------------------------------------------------
    // RETRIEVE THE SOA COLLECTIONS TO BE USED IN THE KERNEL DEVICE
    // (THE FOLLOWING DATA ARE ALREADY ALPAKA-FRIENDLY)

    auto const& clusters = deviceEvent.get(clusterToken_);
    auto const& digis = deviceEvent.get(digisToken_);

    ////auto const& recHits = deviceEvent.get(recHitsToken_);
    //auto const& zVertices = deviceEvent.get(zVertexToken_);
    auto const& candidates = deviceEvent.get(candidateToken_);
    auto const& clustergeometry = deviceEvent.get(geometryToken_);
    auto const& vertices = deviceEvent.get(vertices_);

    const reco::Vertex& pv = vertices[0];

    float vertexX = pv.position().x();
    float vertexY = pv.position().y();
    float vertexZ = pv.position().z();
    GlobalPoint ppv(vertexX, vertexY, vertexZ);
    float vertexEta = ppv.eta();
    float vertexPhi = ppv.phi();


    ///if (verbose_) std::cout << "All Things retrievied..." << std::endl;

    // Use event ID as the offset
    int32_t eventOffset = deviceEvent.id().event();
    ///if (verbose_) std::cout << "Event offset: " << eventOffset << std::endl;
    for (const auto& device : devices_) {
        Queue queue(device);

        // ------------- CREATE DEVICE BUFFERS -------------------------------

        ////size_t nHits = recHits.nHits();
        ////TrackingRecHitsSoACollection<pixelTopology::Phase1> tkHit(queue, nHits, eventOffset, moduleStartD.data());
        ////if (verbose_) std::cout << "TrackingRecHitsSoACollection done " << nHits << std::endl;
        //- - - - - - - - - - - - - - - - - - -

        /* Digis 
        the SiPixelDigisSoACollection is an alias for: SiPixelDigisDevice (gpu) or 
                                                          SiPixelDigisHost (cpu)
        but it's not templated so <pixelTopology> won't work
        I could also: SiPixelDigisDevice<Device> digisDevice(nDigis, queue); */

        size_t nDigis = digis.view().metadata().size()-1;
        SiPixelDigisSoACollection tkDigi(nDigis, queue);
        ///if (verbose_) std::cout << "SiPixelDigisSoACollection done " << nDigis << std::endl;

        //- - - - - - - - - - - - - - - - - - -

        /* Clusters
           the SiPixelClustersSoACollection is an alias for: SiPixelClustersDevice (gpu) 
                                                             SiPixelClustersHost (cpu)  */
        //size_t nClusters = clusters.view().metadata().size()-1;
        //SiPixelClustersSoACollection tkClusters(nClusters, queue); // It seems the above class has no topology and no Modules.. not sure why
        ///if (verbose_) std::cout << "SiPixelClustersSoACollection done " << nClusters << std::endl;

        //alpaka::memcpy(queue, tkClusters.buffer(), clusters.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking


        /* Candidates*/
        size_t nCandidates = candidates.view().metadata().size()-1;
        CandidatesSoACollection tkCandidates(nCandidates, queue);
        auto CandidatesdeviceView = tkCandidates.view();
        ///if (verbose_) std::cout << "CandidatesSoACollection done " << nCandidates << std::endl;
        //- - - - - - - - - - - - - - - - - - -


        /* Geometry*/
        size_t ngeoClusters = clustergeometry.view().metadata().size()-1;
        ClusterGeometrysSoACollection tkgeoclusters(ngeoClusters, queue);
        auto deviceView = tkgeoclusters.view();
        ///if (verbose_) std::cout << "ClusterGeometrysSoACollection done " << ngeoClusters << std::endl;
        //- - - - - - - - - - - - - - - - - - -



        /* Vertices                    */
        //ZVertexSoACollection tkVertices(queue);
        //- - - - - - - - - - - - - - - - - - -

        /* SoA for the output                    */
        SiPixelDigisSoACollection tkOutputDigis(nDigis, queue);
        //SiPixelClustersSoACollection tkOutputClusters(nClusters, queue);
        //if (verbose_) std::cout << "SoA for the output done" << std::endl;

        // ------------- COPY FROM HOST TO DEVICE BUFFERS -------------------------------
        // The output SoA are initialized with the input ones (in case no cluster will be split)

        ////alpaka::memcpy(queue, tkHit.buffer(), recHits.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking

        
        alpaka::memcpy(queue, tkDigi.buffer(), digis.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking

        //alpaka::memcpy(queue, tkClusters.buffer(), clusters.buffer());
        //alpaka::memcpy(queue, tkVertices.buffer(), zVertices.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking

        alpaka::memcpy(queue, tkCandidates.buffer(), candidates.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking

        alpaka::memcpy(queue, tkgeoclusters.buffer(), clustergeometry.buffer());
        //alpaka::wait(queue);  // Ensure copy is finished before checking

        ///if (verbose_) std::cout << "Most memcpy done" << std::endl;
/*
        // Handling the per cluster calculation attributes in a struct
        std::vector<clusterProperties> gpuAlgo;
        auto clusterPropertiesHost = cms::alpakatools::make_host_buffer<clusterProperties[]>(queue, nClusters);
        auto clusterPropertiesDevice = cms::alpakatools::make_device_buffer<clusterProperties[]>(queue, nClusters);
        std::copy(gpuAlgo.begin(), gpuAlgo.end(), clusterPropertiesHost.data());
        alpaka::memcpy(queue, clusterPropertiesDevice, clusterPropertiesHost);
*/
//        alpaka::wait(queue);  // Ensure copy is finished before checking
        
        ///if (verbose_) std::cout << "All memcpy done" << std::endl;

        // Handling a global counter of the output (new) clusters (initialized to zero here)
        auto clusterCounterDevice = cms::alpakatools::make_device_buffer<uint32_t>(queue);
        alpaka::memset(queue, clusterCounterDevice, 0);

        auto pixelCounterDevice = cms::alpakatools::make_device_buffer<uint32_t>(queue);
        alpaka::memset(queue, pixelCounterDevice, 0);
        alpaka::wait(queue);  // Ensure the transfer is complete

        std::vector<uint16_t> tinyClusters;
        std::vector<uint16_t> smallClusters;
        std::vector<uint16_t> mediumClusters;
        std::vector<uint16_t> largeClusters;
        std::vector<uint16_t> heavyClusters;

        uint16_t pixelTinyThreshold = 255;
        uint16_t pixelLowThreshold = 15;
        uint16_t pixelMediumThreshold = 31;
        uint16_t pixelLargeThreshold = 127;
        uint16_t pixelHeavyThreshold = 255;

        for (size_t clusterID = 0; clusterID < clusterPixelCounts.size(); ++clusterID) {
            uint32_t pixelCount = clusterPixelCounts[clusterID];

            if (pixelCount <= pixelTinyThreshold) {
                tinyClusters.push_back(clusterID);
            }

            if (pixelCount <= pixelLowThreshold) {
                smallClusters.push_back(clusterID);
            }
            else if (pixelCount <= pixelMediumThreshold) {
                mediumClusters.push_back(clusterID);
            }
            else if (pixelCount <= pixelLargeThreshold) {
                largeClusters.push_back(clusterID);
            }
            else {
                //std::cout << "Large cluster: " << clusterID << " with " << pixelCount << " Pixels " << std::endl;
                heavyClusters.push_back(clusterID);
            }
        }

        if (!tinyClusters.empty()) {

            Splitting::runKernels<pixelTopology::Phase1>(
                  tkDigi.view(), //tkClusters.view(), 
                  tkCandidates.view(), tkgeoclusters.view(), ptMin_, deltaR_, chargeFracMin_,
                  expSizeXAtLorentzAngleIncidence_, expSizeXDeltaPerTanAlpha_, expSizeYAtNormalIncidence_,
                  centralMIPCharge_, chargePerUnit_, fractionalWidth_,
                  tkOutputDigis.view(), 
                  //tkOutputClusters.view(),
                  clusterCounterDevice.data(), pixelCounterDevice.data(),
                  vertexX, vertexY, vertexZ, vertexEta, vertexPhi,
                  verbose_, debugMode, targetDetId, targetClusterOffset,
                  tinyClusters.data(), tinyClusters.size(), pixelTinyThreshold, queue);
        }
        alpaka::wait(queue);  // Ensure the transfer is complete

/*

        if (!smallClusters.empty()) {

            Splitting::runKernels<pixelTopology::Phase1>(
                  tkDigi.view(), //tkClusters.view(), 
                  tkCandidates.view(), tkgeoclusters.view(), ptMin_, deltaR_, chargeFracMin_,
                  expSizeXAtLorentzAngleIncidence_, expSizeXDeltaPerTanAlpha_, expSizeYAtNormalIncidence_,
                  centralMIPCharge_, chargePerUnit_, fractionalWidth_,
                  tkOutputDigis.view(), 
                  //tkOutputClusters.view(),
                  clusterCounterDevice.data(), pixelCounterDevice.data(),
                  vertexX, vertexY, vertexZ, vertexEta, vertexPhi,
                  verbose_, debugMode, targetDetId, targetClusterOffset,
                  smallClusters.data(), smallClusters.size(), pixelLowThreshold, queue);
        }
        alpaka::wait(queue);  // Ensure the transfer is complete


        // Single call for medium clusters
        if (!mediumClusters.empty()) {
            Splitting::runKernels<pixelTopology::Phase1>(
                  tkDigi.view(), //tkClusters.view(), 
                  tkCandidates.view(), tkgeoclusters.view(), ptMin_, deltaR_, chargeFracMin_,
                  expSizeXAtLorentzAngleIncidence_, expSizeXDeltaPerTanAlpha_, expSizeYAtNormalIncidence_,
                  centralMIPCharge_, chargePerUnit_, fractionalWidth_,
                  tkOutputDigis.view(), 
                  //tkOutputClusters.view(),
                  clusterCounterDevice.data(), pixelCounterDevice.data(),
                  vertexX, vertexY, vertexZ, vertexEta, vertexPhi,
                  verbose_, debugMode, targetDetId, targetClusterOffset,
                  mediumClusters.data(), mediumClusters.size(), pixelMediumThreshold, queue);
        }        
        alpaka::wait(queue);  // Ensure the transfer is complete
        // Single call for large clusters
        if (!largeClusters.empty()) {
            Splitting::runKernels<pixelTopology::Phase1>(
                  tkDigi.view(), //tkClusters.view(), 
                  tkCandidates.view(), tkgeoclusters.view(), ptMin_, deltaR_, chargeFracMin_,
                  expSizeXAtLorentzAngleIncidence_, expSizeXDeltaPerTanAlpha_, expSizeYAtNormalIncidence_,
                  centralMIPCharge_, chargePerUnit_, fractionalWidth_,
                  tkOutputDigis.view(), 
                  //tkOutputClusters.view(),
                  clusterCounterDevice.data(), pixelCounterDevice.data(),
                  vertexX, vertexY, vertexZ, vertexEta, vertexPhi,
                  verbose_, debugMode, targetDetId, targetClusterOffset,
                  largeClusters.data(), largeClusters.size(), pixelLargeThreshold, queue);
        }
        alpaka::wait(queue);  // Ensure the transfer is complete

        if (!heavyClusters.empty()) {
            Splitting::runKernels<pixelTopology::Phase1>(
                  tkDigi.view(), //tkClusters.view(), 
                  tkCandidates.view(), tkgeoclusters.view(), ptMin_, deltaR_, chargeFracMin_,
                  expSizeXAtLorentzAngleIncidence_, expSizeXDeltaPerTanAlpha_, expSizeYAtNormalIncidence_,
                  centralMIPCharge_, chargePerUnit_, fractionalWidth_,
                  tkOutputDigis.view(), 
                  //tkOutputClusters.view(),
                  clusterCounterDevice.data(), pixelCounterDevice.data(),
                  vertexX, vertexY, vertexZ, vertexEta, vertexPhi,
                  verbose_, debugMode, targetDetId, targetClusterOffset,
                  heavyClusters.data(), heavyClusters.size(), pixelHeavyThreshold, queue);
        }
        alpaka::wait(queue);  // Ensure the transfer is complete
*/
        // Update from device to host
        //alpaka::memcpy(queue, gpuSharedHost, gpuSharedDevice);  // Copy device buffer to host buffer
        //alpaka::wait(queue);  // Ensure the transfer is complete
//        tkHit.updateFromDevice(queue);

        //TrackingRecHitHost<pixelTopology::Phase1> hostRecHits = cms::alpakatools::CopyToHost<TrackingRecHitDevice<pixelTopology::Phase1, Device>>::copyAsync(queue, tkHit);
        SiPixelDigisHost digisHost = cms::alpakatools::CopyToHost<SiPixelDigisDevice<Device>>::copyAsync(queue, tkDigi);
        //SiPixelClustersHost clustersHost = cms::alpakatools::CopyToHost<SiPixelClustersDevice<Device>>::copyAsync(queue, tkClusters);
//        SiPixelDigisHost outputDigisHost = cms::alpakatools::CopyToHost<SiPixelDigisDevice<Device>>::copyAsync(queue, tkOutputDigis);
//        SiPixelClustersHost outputClustersHost = cms::alpakatools::CopyToHost<SiPixelClustersDevice<Device>>::copyAsync(queue, tkOutputClusters);

        alpaka::wait(queue);
        //deviceEvent.put(std::move(digisHost));
        deviceEvent.emplace(outputdigisToken_, std::move(tkOutputDigis) );
    }
       
}


void trial::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    ////desc.add<uint32_t>("nHits", 100)->setComment("Number of hits for the test");
    ////desc.add<int32_t>("offset", 0)->setComment("Offset for hits");
    desc.add<double>("ptMin", 100.0)->setComment("Minimum pt");
    desc.add<double>("deltaR", 0.05)->setComment("Delta R");
    desc.add<double>("chargeFracMin", 2.0)->setComment("Minimum charge fraction");
    desc.add<double>("tanLorentzAngle", 0.0)->setComment("Lorentz angle");
    desc.add<double>("tanLorentzAngleBarrelLayer1", 0.0)->setComment("Lorentz angle for Barrel Layer 1");
    desc.add<double>("expSizeXAtLorentzAngleIncidence", 1.5)->setComment("Expected size X at Lorentz angle incidence");
    desc.add<double>("expSizeXDeltaPerTanAlpha", 0.0)->setComment("Expected size X delta per tan(alpha)");
    desc.add<double>("expSizeYAtNormalIncidence", 1.3)->setComment("Expected size Y at normal incidence");
    desc.add<double>("centralMIPCharge", 26000.0)->setComment("Central MIP charge");
    desc.add<double>("chargePerUnit", 2000.0)->setComment("Charge per unit");
    desc.add<double>("forceXError", 100.0)->setComment("Force X error");
    desc.add<double>("forceYError", 150.0)->setComment("Force Y error");
    desc.add<double>("fractionalWidth", 0.4)->setComment("Fractional width");
    desc.add<edm::InputTag>("candidateInput", edm::InputTag(""))->setComment("Input tag for candidate data");
    desc.add<edm::InputTag>("geometryInput", edm::InputTag(""))->setComment("Input tag for geometry data");
    desc.add<edm::InputTag>("siPixelClusters", edm::InputTag(""))->setComment("Input tag for siPixelClusters data");
    desc.add<edm::InputTag>("siPixelDigis", edm::InputTag(""))->setComment("Input tag for siPixelDigis data");
    ////desc.add<edm::InputTag>("trackingRecHits", edm::InputTag(""))->setComment("Input tag for trackingRecHits data");
    //desc.add<edm::InputTag>("zVertex", edm::InputTag(""))->setComment("Input tag for zVertex data");
    desc.add<bool>("verbose", false)->setComment("Verbose output");
    desc.add<bool>("debugMode", true);
    desc.add<int>("targetDetId");
    desc.add<int>("targetClusterOffset");
    desc.add<int>("targetEvent");    
    desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));    
    desc.add<edm::InputTag>("clusterPixelCounts", edm::InputTag("clusterPixelCounts"));    
    descriptions.add("trial", desc);
}


DEFINE_FWK_MODULE(trial);

#endif  // ALPAKA_ACC_GPU_CUDA_ENABLED
