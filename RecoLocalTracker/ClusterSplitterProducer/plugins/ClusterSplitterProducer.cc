#ifdef ALPAKA_ACC_GPU_CUDA_ENABLED

#include <memory>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/stringize.h"

#include "TFile.h"
#include "TString.h"
#include <vector>
#include <tuple>
#include <cstdlib>
#include <unistd.h>

#include <alpaka/alpaka.hpp>

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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/devices.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/ClusterGeometrySoA/interface/ClusterGeometrysSoA.h"
#include "DataFormats/ClusterGeometrySoA/interface/alpaka/ClusterGeometrysSoACollection.h"

#include "DataFormats/CandidateSoA/interface/CandidatesSoA.h"
#include "DataFormats/CandidateSoA/interface/alpaka/CandidatesSoACollection.h"

using namespace ALPAKA_ACCELERATOR_NAMESPACE;

class HelperSplitter : public edm::stream::EDProducer<> {
public:
  explicit HelperSplitter(const edm::ParameterSet&);
  ~HelperSplitter() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  const double ptMin_;
  float tanLorentzAngle_;
  float tanLorentzAngleBarrelLayer1_;  

  edm::EDGetTokenT<SiPixelClusterCollectionNew> clusterToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> candidateToken_;
  edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> const tTrackingGeom_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> const tTrackerTopo_;
  bool verbose_;

};

HelperSplitter::HelperSplitter(const edm::ParameterSet& iConfig)
    : ptMin_(iConfig.getParameter<double>("ptMin")),
      tanLorentzAngle_(iConfig.getParameter<double>("tanLorentzAngle")),
      tanLorentzAngleBarrelLayer1_(iConfig.getParameter<double>("tanLorentzAngleBarrelLayer1")),
      clusterToken_(consumes<SiPixelClusterCollectionNew>(iConfig.getParameter<edm::InputTag>("siPixelClusters"))),
      //candidateToken_(consumes<edm::View<reco::Candidate>>(edm::InputTag("Candidate"))),
      candidateToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("Candidate"))),
      tTrackingGeom_(esConsumes()),
      tTrackerTopo_(esConsumes()),
      verbose_(iConfig.getParameter<bool>("verbose"))      
        {}

HelperSplitter::~HelperSplitter() {
}

void HelperSplitter::beginStream(edm::StreamID) {
}

void HelperSplitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    printf("*********************************Starting the HelperSplitter producer.\n");
    //std::cerr << "\n\n\n*************** HELPER SPLITTER IS RUNNING ***************\n\n\n";

    // Candidate is used for retrieving the jets
    edm::Handle<edm::View<reco::Candidate>> candidatesHandle;
    iEvent.getByToken(candidateToken_, candidatesHandle);
    if (!candidatesHandle.isValid()) {
        edm::LogError("HelperSplitter") << "Could not retrieve Candidate";
        return;
    }

    // Process Candidates
    size_t nCandidates = candidatesHandle->size();
    if (verbose_) {
        std::cout << "Number of Candidates: " << nCandidates << std::endl;
    }

    // Count the number of valid candidates that pass the ptMin_ filter
    size_t validCandidatesCount = 0;
    for (const auto& candidate : *candidatesHandle) {
        if (candidate.pt() > ptMin_) {
            ++validCandidatesCount;
        }
    }
    if (verbose_) std::cout << "Number of valid Candidates: " << validCandidatesCount << std::endl;

    // Create the queue for the CPU device
    auto const& device = cms::alpakatools::devices<Platform>()[0];
    Queue queue(device);
    if (verbose_) std::cout << "Queue done" << std::endl;

    // Create the CandidateSoA on CPU
    CandidatesHost tkCandidates(nCandidates, queue);
    auto candidateView = tkCandidates.view();
    if (verbose_) std::cout << "Candidates done" << std::endl;

    // Fill the CandidateSoA
    size_t candidateIndex = 0;
    for (const auto& candidate : *candidatesHandle) {
        if (candidate.pt() > ptMin_) {  // Apply the ptMin_ filter
            candidateView.px(candidateIndex) = static_cast<float>(candidate.px());
            candidateView.py(candidateIndex) = static_cast<float>(candidate.py());
            candidateView.pz(candidateIndex) = static_cast<float>(candidate.pz());
            candidateView.pt(candidateIndex) = static_cast<float>(candidate.pt());
            candidateView.eta(candidateIndex) = static_cast<float>(candidate.eta());
            candidateView.phi(candidateIndex) = static_cast<float>(candidate.phi());
            ++candidateIndex;
        }
    }
    if (verbose_) std::cout << "Done with Candidates" << std::endl;

    // SiPixelClusters is used to get the geometry of each cluster
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> inputPixelClustersHandle;
    iEvent.getByToken(clusterToken_, inputPixelClustersHandle);
    if (!inputPixelClustersHandle.isValid()) {
        edm::LogError("HelperSplitter") << "Could not retrieve siPixelClusters.";
        return;
    }
    if (verbose_) std::cout << "siPixelClusters got it" << std::endl;


    // Process inputPixelClustersHandle
    size_t nPixelClusters = inputPixelClustersHandle->size();
    if (verbose_) {
        std::cout << "Number of Pixels: " << nPixelClusters << std::endl;
    }
    

    // Retrieve TrackerGeometry, trackerTopology from EventSetup
    const auto& trackingGeometry = iSetup.getData(tTrackingGeom_);
    const auto& trackerTopology = iSetup.getData(tTrackerTopo_);
    if (verbose_) std::cout << "TrackerGeometry/Topology got it" << std::endl;

    // Create the ClusterGeometrySoA on CPU
    ClusterGeometrysHost tkCluster(nPixelClusters, queue);
    auto clusterView = tkCluster.view();
    if (verbose_) std::cout << "Cluster done" << std::endl;

    for (auto detIt = inputPixelClustersHandle->begin(); detIt != inputPixelClustersHandle->end(); ++detIt) {
        const edmNew::DetSet<SiPixelCluster>& detset = *detIt;
        const GeomDet* det = trackingGeometry.idToDet(detset.id());
        if (!det) continue;

        const PixelTopology& topo = static_cast<const PixelTopology&>(det->topology());
        float pitchX, pitchY;
        std::tie(pitchX, pitchY) = topo.pitch();
        float thickness = det->surface().bounds().thickness();
        float tanLorentzAngle = tanLorentzAngle_;

        size_t clusterIndex = 0;
        for (const auto& cluster : detset) {
            clusterView.clusterIds(clusterIndex) = detset.id();
            clusterView.pitchX(clusterIndex) = pitchX;
            clusterView.pitchY(clusterIndex) = pitchY;
            clusterView.thickness(clusterIndex) = thickness;
            clusterView.tanLorentzAngles(clusterIndex) = tanLorentzAngle;
            ++clusterIndex;
        }
    }

    // Put the CandidateSoA and ClusterGeometrySoA into the event
    //iEvent.put(std::move(tkCandidates), "CandidateSoA");
    //iEvent.put(std::move(clusterDataSoA), "ClusterGeometrySoA");
}


void HelperSplitter::endStream() {
  edm::LogInfo("HelperSplitter") << "Processing completed.";
}

void HelperSplitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

    edm::ParameterSetDescription desc;
    desc.add<bool>("verbose", false)->setComment("Verbose output");
    desc.add<double>("ptMin", 0.5)->setComment("Minimum pt for filtering candidates");
    desc.add<double>("tanLorentzAngle", 0.1)->setComment("Lorentz angle tangent");
    desc.add<double>("tanLorentzAngleBarrelLayer1", 0.2)->setComment("Lorentz angle tangent for Barrel Layer 1");
    desc.add<edm::InputTag>("siPixelClusters", edm::InputTag("siPixelClusters"))->setComment("Collection for siPixelClusters");
    desc.add<edm::InputTag>("Candidate", edm::InputTag("Candidate"))->setComment("Candidates");
    descriptions.add("HelperSplitter", desc);
}

DEFINE_FWK_MODULE(HelperSplitter);

#endif
