#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "genmatching_utilities.h"
#include "propagate_utilities.h"

typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;

namespace MTYPE {
const char* DSA = "DSA";
const char* DGL = "DGL";
}  // namespace MTYPE

bool passTagID(const reco::Track* track, const char* mtype) {
    bool passID = false;
    if (mtype == MTYPE::DSA) {
        if (track->phi() >= -0.8 || track->phi() <= -2.1) {
            return passID;
        }
        if (abs(track->eta()) >= 0.7) {
            return passID;
        }
        if (track->pt() <= 12.5) {
            return passID;
        }
        if (track->ptError() / track->pt() >= 0.2) {
            return passID;
        }
        if (track->hitPattern().numberOfValidMuonDTHits() <= 30) {
            return passID;
        }
        if (track->normalizedChi2() >= 2) {
            return passID;
        }
        passID = true;
    } else if (mtype == MTYPE::DGL) {
        if (track->phi() >= -0.6 || track->phi() <= -2.6) {
            return passID;
        }
        if (abs(track->eta()) >= 0.9) {
            return passID;
        }
        if (track->pt() <= 20) {
            return passID;
        }
        if (track->ptError() / track->pt() >= 0.3) {
            return passID;
        }
        if (track->hitPattern().numberOfMuonHits() <= 12) {
            return passID;
        }
        if (track->hitPattern().numberOfValidStripHits() <= 5) {
            return passID;
        }
        passID = true;
    } else {
        std::cout << "Error (in passTagID): wrong muon type" << std::endl;
    }
    return passID;
}

bool passProbeID(const reco::Track* track, const TVector3& v_tag, const char* mtype) {
    bool passID = false;
    if (mtype == MTYPE::DSA) {
        if (track->hitPattern().numberOfValidMuonDTHits() +
                track->hitPattern().numberOfValidMuonCSCHits() <=
            12) {
            return passID;
        }
        if (track->pt() <= 3.5) {
            return passID;
        }
        TVector3 v_probe = TVector3();
        v_probe.SetPtEtaPhi(track->pt(), track->eta(), track->phi());
        if (v_probe.Angle(v_tag) <= 2.1) {
            return passID;
        }
        passID = true;
    } else if (mtype == MTYPE::DGL) {
        if (track->pt() <= 20) {
            return passID;
        }
        TVector3 v_probe = TVector3();
        v_probe.SetPtEtaPhi(track->pt(), track->eta(), track->phi());
        if (v_probe.Angle(v_tag) <= 2.8) {
            return passID;
        }
        passID = true;
    } else {
        std::cout << "Error (in passProbeID): wrong muon type" << std::endl;
    }
    return passID;
}

class ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
    explicit ntuplizer(const edm::ParameterSet&);
    ~ntuplizer();

    edm::ConsumesCollector iC = consumesCollector();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::ParameterSet parameters;

    bool isAOD = false;
    bool isCosmics = false;
    bool isMCSignal = false;

    // std::ofstream residualsfile;

    //
    // --- Tokens and Handles
    //

    // trigger bits
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::Handle<edm::TriggerResults> triggerBits;

    // displacedGlobalMuons (reco::Track)
    edm::EDGetTokenT<edm::View<reco::Track>> dglToken;
    edm::Handle<edm::View<reco::Track>> dgls;
    // displacedStandAloneMuons (reco::Track)
    edm::EDGetTokenT<edm::View<reco::Track>> dsaToken;
    edm::Handle<edm::View<reco::Track>> dsas;
    // displacedMuons (reco::Muon // pat::Muon)
    edm::EDGetTokenT<edm::View<reco::Muon>> dmuToken;
    edm::Handle<edm::View<reco::Muon>> dmuons;
    // prunedGenParticles (reco::GenParticle)
    edm::EDGetTokenT<edm::View<reco::GenParticle>> prunedGenToken;
    edm::Handle<edm::View<reco::GenParticle>> prunedGen;

    // Propagator
    edm::ESGetToken<Propagator, TrackingComponentsRecord> thePropAlongToken;
    edm::ESGetToken<Propagator, TrackingComponentsRecord> thePropOppositeToken;

    // Magnetic Field
    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> theMagFieldToken;
    edm::ESHandle<MagneticField> magField_;

    // Trigger tags
    std::vector<std::string> HLTPaths_;
    bool triggerPass[200] = {false};

    // Event
    Int_t event = 0;
    Int_t lumiBlock = 0;
    Int_t run = 0;

    // ----------------------------------
    // displacedMuons
    // ----------------------------------
    Int_t ndmu = 0;
    Int_t dmu_isDSA[200] = {0};
    Int_t dmu_isDGL[200] = {0};
    Int_t dmu_isDTK[200] = {0};
    Int_t dmu_isMatchesValid[200] = {0};
    Int_t dmu_numberOfMatches[200] = {0};
    Int_t dmu_numberOfChambers[200] = {0};
    Int_t dmu_numberOfChambersCSCorDT[200] = {0};
    Int_t dmu_numberOfMatchedStations[200] = {0};
    Int_t dmu_numberOfMatchedRPCLayers[200] = {0};

    // Variables by Marco
    Float_t dmu_t0_InOut[200] = {0.};
    Float_t dmu_t0_OutIn[200] = {0.};
    Float_t dmu_vertex_r[200] = {0.};
    Float_t dmu_vertex_z[200] = {0.};
    Float_t genmu_vertex_r[200] = {0.};
    Float_t genmu_vertex_z[200] = {0.};

    // Variables for gen matching
    bool dmu_hasGenMatch[200] = {false};
    Int_t dmu_genMatchedIndex[200] = {0};
    Float_t dmu_residual_r[200] = {0.};
    Float_t dmu_residual_theta[200] = {0.};
    Float_t dmu_residual_phi[200] = {0.};
    Float_t dmu_residual_p_r[200] = {0.};
    Float_t dmu_residual_p_theta[200] = {0.};
    Float_t dmu_residual_p_phi[200] = {0.};
    Float_t dmu_residual_x[200] = {0.};
    Float_t dmu_residual_y[200] = {0.};
    Float_t dmu_residual_z[200] = {0.};
    Float_t dmu_residual_p_x[200] = {0.};
    Float_t dmu_residual_p_y[200] = {0.};
    Float_t dmu_residual_p_z[200] = {0.};
    Int_t dmu_propagationSurface[200] = {0};

    Float_t dmu_dsa_pt[200] = {0.};
    Float_t dmu_dsa_eta[200] = {0.};
    Float_t dmu_dsa_phi[200] = {0.};
    Float_t dmu_dsa_ptError[200] = {0.};
    Float_t dmu_dsa_dxy[200] = {0.};
    Float_t dmu_dsa_dz[200] = {0.};
    Float_t dmu_dsa_normalizedChi2[200] = {0.};
    Float_t dmu_dsa_charge[200] = {0.};
    Int_t dmu_dsa_nMuonHits[200] = {0};
    Int_t dmu_dsa_nValidMuonHits[200] = {0};
    Int_t dmu_dsa_nValidMuonDTHits[200] = {0};
    Int_t dmu_dsa_nValidMuonCSCHits[200] = {0};
    Int_t dmu_dsa_nValidMuonRPCHits[200] = {0};
    Int_t dmu_dsa_nValidStripHits[200] = {0};
    Int_t dmu_dsa_nhits[200] = {0};
    Int_t dmu_dsa_dtStationsWithValidHits[200] = {0};
    Int_t dmu_dsa_cscStationsWithValidHits[200] = {0};
    Int_t dmu_dsa_nsegments[200] = {0};

    Float_t dmu_dgl_pt[200] = {0.};
    Float_t dmu_dgl_eta[200] = {0.};
    Float_t dmu_dgl_phi[200] = {0.};
    Float_t dmu_dgl_ptError[200] = {0.};
    Float_t dmu_dgl_dxy[200] = {0.};
    Float_t dmu_dgl_dz[200] = {0.};
    Float_t dmu_dgl_normalizedChi2[200] = {0.};
    Float_t dmu_dgl_charge[200] = {0.};
    Int_t dmu_dgl_nMuonHits[200] = {0};
    Int_t dmu_dgl_nValidMuonHits[200] = {0};
    Int_t dmu_dgl_nValidMuonDTHits[200] = {0};
    Int_t dmu_dgl_nValidMuonCSCHits[200] = {0};
    Int_t dmu_dgl_nValidMuonRPCHits[200] = {0};
    Int_t dmu_dgl_nValidStripHits[200] = {0};
    Int_t dmu_dgl_nhits[200] = {0};

    // Variables for tag and probe
    bool dmu_dsa_passTagID[200] = {false};
    bool dmu_dsa_hasProbe[200] = {false};
    Int_t dmu_dsa_probeID[200] = {0};
    Float_t dmu_dsa_cosAlpha[200] = {0.};
    bool dmu_dsa_isProbe[200] = {false};

    bool dmu_dgl_passTagID[200] = {false};
    bool dmu_dgl_hasProbe[200] = {false};
    Int_t dmu_dgl_probeID[200] = {0};
    Float_t dmu_dgl_cosAlpha[200] = {0.};
    bool dmu_dgl_isProbe[200] = {false};

    Float_t dmu_dtk_pt[200] = {0.};
    Float_t dmu_dtk_eta[200] = {0.};
    Float_t dmu_dtk_phi[200] = {0.};
    Float_t dmu_dtk_ptError[200] = {0.};
    Float_t dmu_dtk_dxy[200] = {0.};
    Float_t dmu_dtk_dz[200] = {0.};
    Float_t dmu_dtk_normalizedChi2[200] = {0.};
    Float_t dmu_dtk_charge[200] = {0.};
    Int_t dmu_dtk_nMuonHits[200] = {0};
    Int_t dmu_dtk_nValidMuonHits[200] = {0};
    Int_t dmu_dtk_nValidMuonDTHits[200] = {0};
    Int_t dmu_dtk_nValidMuonCSCHits[200] = {0};
    Int_t dmu_dtk_nValidMuonRPCHits[200] = {0};
    Int_t dmu_dtk_nValidStripHits[200] = {0};
    Int_t dmu_dtk_nhits[200] = {0};

    //
    // --- Output
    //
    std::string output_filename;
    TH1F* counts;
    TFile* file_out;
    TTree* tree_out;
};

// Constructor
ntuplizer::ntuplizer(const edm::ParameterSet& iConfig) {
    usesResource("TFileService");

    parameters = iConfig;

    // Analyzer parameters
    isAOD = parameters.getParameter<bool>("isAOD");
    isCosmics = parameters.getParameter<bool>("isCosmics");
    isMCSignal = parameters.getParameter<bool>("isMCSignal");

    counts = new TH1F("counts", "", 1, 0, 1);

    dglToken = consumes<edm::View<reco::Track>>(
        parameters.getParameter<edm::InputTag>("displacedGlobalCollection"));
    dsaToken = consumes<edm::View<reco::Track>>(
        parameters.getParameter<edm::InputTag>("displacedStandAloneCollection"));
    dmuToken = consumes<edm::View<reco::Muon>>(
        parameters.getParameter<edm::InputTag>("displacedMuonCollection"));
    if (isMCSignal) {
        prunedGenToken = consumes<edm::View<reco::GenParticle>>(
            parameters.getParameter<edm::InputTag>("prunedGenParticles"));
        thePropAlongToken = esConsumes(
            edm::ESInputTag("", parameters.getParameter<std::string>("propagatorAlong")));
        thePropOppositeToken = esConsumes(
            edm::ESInputTag("", parameters.getParameter<std::string>("propagatorOpposite")));
        theMagFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
    }

    triggerBits_ = consumes<edm::TriggerResults>(parameters.getParameter<edm::InputTag>("bits"));
}

// Destructor
ntuplizer::~ntuplizer() {}

// beginJob (Before first event)
void ntuplizer::beginJob() {
    std::cout << "Begin Job" << std::endl;

    // Open the chi2 file
    // residualsfile.open("residuals.csv");
    // residualsfile << "R, theta, phi, p_R, p_theta, p_phi, x, y, z, p_x, p_y, p_z\n";

    // Init the file and the TTree
    output_filename = parameters.getParameter<std::string>("nameOfOutput");
    file_out = new TFile(output_filename.c_str(), "RECREATE");
    tree_out = new TTree("Events", "Events");

    // Load HLT paths
    HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX3BX");
    HLTPaths_.push_back("HLT_L2Mu10_NoVertex_NoBPTX");

    // TTree branches
    tree_out->Branch("event", &event, "event/I");
    tree_out->Branch("lumiBlock", &lumiBlock, "lumiBlock/I");
    tree_out->Branch("run", &run, "run/I");

    // ----------------------------------
    // displacedMuons
    // ----------------------------------
    tree_out->Branch("ndmu", &ndmu, "ndmu/I");
    tree_out->Branch("dmu_isDSA", dmu_isDSA, "dmu_isDSA[ndmu]/I");
    tree_out->Branch("dmu_isDGL", dmu_isDGL, "dmu_isDGL[ndmu]/I");
    tree_out->Branch("dmu_isDTK", dmu_isDTK, "dmu_isDTK[ndmu]/I");
    tree_out->Branch("dmu_isMatchesValid", dmu_isMatchesValid, "dmu_isMatchesValid[ndmu]/I");
    tree_out->Branch("dmu_numberOfMatches", dmu_numberOfMatches, "dmu_numberOfMatches[ndmu]/I");
    tree_out->Branch("dmu_numberOfChambers", dmu_numberOfChambers, "dmu_numberOfChambers[ndmu]/I");
    tree_out->Branch("dmu_numberOfChambersCSCorDT", dmu_numberOfChambersCSCorDT,
                     "dmu_numberOfChambersCSCorDT[ndmu]/I");
    tree_out->Branch("dmu_numberOfMatchedStations", dmu_numberOfMatchedStations,
                     "dmu_numberOfMatchedStations[ndmu]/I");
    tree_out->Branch("dmu_numberOfMatchedRPCLayers", dmu_numberOfMatchedRPCLayers,
                     "dmu_numberOfMatchedRPCLayers[ndmu]/I");
    tree_out->Branch("dmu_t0_InOut", dmu_t0_InOut, "dmu_t0_InOut[ndmu]/F");
    tree_out->Branch("dmu_t0_OutIn", dmu_t0_OutIn, "dmu_t0_OutIn[ndmu]/F");
    tree_out->Branch("dmu_vertex_r", dmu_vertex_r, "dmu_vertex_r[ndmu]/F");
    tree_out->Branch("dmu_vertex_z", dmu_vertex_z, "dmu_vertex_z[ndmu]/F");
    tree_out->Branch("genmu_vertex_r", genmu_vertex_r, "genmu_vertex_r[ndmu]/F");
    tree_out->Branch("genmu_vertex_z", genmu_vertex_z, "genmu_vertex_z[ndmu]/F");
    // dmu_dsa
    tree_out->Branch("dmu_dsa_pt", dmu_dsa_pt, "dmu_dsa_pt[ndmu]/F");
    tree_out->Branch("dmu_dsa_eta", dmu_dsa_eta, "dmu_dsa_eta[ndmu]/F");
    tree_out->Branch("dmu_dsa_phi", dmu_dsa_phi, "dmu_dsa_phi[ndmu]/F");
    tree_out->Branch("dmu_dsa_ptError", dmu_dsa_ptError, "dmu_dsa_ptError[ndmu]/F");
    tree_out->Branch("dmu_dsa_dxy", dmu_dsa_dxy, "dmu_dsa_dxy[ndmu]/F");
    tree_out->Branch("dmu_dsa_dz", dmu_dsa_dz, "dmu_dsa_dz[ndmu]/F");
    tree_out->Branch("dmu_dsa_normalizedChi2", dmu_dsa_normalizedChi2,
                     "dmu_dsa_normalizedChi2[ndmu]/F");
    tree_out->Branch("dmu_dsa_charge", dmu_dsa_charge, "dmu_dsa_charge[ndmu]/F");
    tree_out->Branch("dmu_dsa_nMuonHits", dmu_dsa_nMuonHits, "dmu_dsa_nMuonHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nValidMuonHits", dmu_dsa_nValidMuonHits,
                     "dmu_dsa_nValidMuonHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nValidMuonDTHits", dmu_dsa_nValidMuonDTHits,
                     "dmu_dsa_nValidMuonDTHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nValidMuonCSCHits", dmu_dsa_nValidMuonCSCHits,
                     "dmu_dsa_nValidMuonCSCHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nValidMuonRPCHits", dmu_dsa_nValidMuonRPCHits,
                     "dmu_dsa_nValidMuonRPCHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nValidStripHits", dmu_dsa_nValidStripHits,
                     "dmu_dsa_nValidStripHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_nhits", dmu_dsa_nhits, "dmu_dsa_nhits[ndmu]/I");
    tree_out->Branch("dmu_dsa_dtStationsWithValidHits", dmu_dsa_dtStationsWithValidHits,
                     "dmu_dsa_dtStationsWithValidHits[ndmu]/I");
    tree_out->Branch("dmu_dsa_cscStationsWithValidHits", dmu_dsa_cscStationsWithValidHits,
                     "dmu_dsa_cscStationsWithValidHits[ndmu]/I");
    // dmu_dgl
    tree_out->Branch("dmu_dgl_pt", dmu_dgl_pt, "dmu_dgl_pt[ndmu]/F");
    tree_out->Branch("dmu_dgl_eta", dmu_dgl_eta, "dmu_dgl_eta[ndmu]/F");
    tree_out->Branch("dmu_dgl_phi", dmu_dgl_phi, "dmu_dgl_phi[ndmu]/F");
    tree_out->Branch("dmu_dgl_ptError", dmu_dgl_ptError, "dmu_dgl_ptError[ndmu]/F");
    tree_out->Branch("dmu_dgl_dxy", dmu_dgl_dxy, "dmu_dgl_dxy[ndmu]/F");
    tree_out->Branch("dmu_dgl_dz", dmu_dgl_dz, "dmu_dgl_dz[ndmu]/F");
    tree_out->Branch("dmu_dgl_normalizedChi2", dmu_dgl_normalizedChi2,
                     "dmu_dgl_normalizedChi2[ndmu]/F");
    tree_out->Branch("dmu_dgl_charge", dmu_dgl_charge, "dmu_dgl_charge[ndmu]/F");
    tree_out->Branch("dmu_dgl_nMuonHits", dmu_dgl_nMuonHits, "dmu_dgl_nMuonHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nValidMuonHits", dmu_dgl_nValidMuonHits,
                     "dmu_dgl_nValidMuonHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nValidMuonDTHits", dmu_dgl_nValidMuonDTHits,
                     "dmu_dgl_nValidMuonDTHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nValidMuonCSCHits", dmu_dgl_nValidMuonCSCHits,
                     "dmu_dgl_nValidMuonCSCHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nValidMuonRPCHits", dmu_dgl_nValidMuonRPCHits,
                     "dmu_dgl_nValidMuonRPCHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nValidStripHits", dmu_dgl_nValidStripHits,
                     "dmu_dgl_nValidStripHits[ndmu]/I");
    tree_out->Branch("dmu_dgl_nhits", dmu_dgl_nhits, "dmu_dgl_nhits[ndmu]/I");

    // Tag and probe branches
    if (isCosmics) {
        tree_out->Branch("dmu_dsa_nsegments", dmu_dsa_nsegments, "dmu_dsa_nsegments[ndmu]/I");
        tree_out->Branch("dmu_dsa_passTagID", dmu_dsa_passTagID, "dmu_dsa_passTagID[ndmu]/O");
        tree_out->Branch("dmu_dsa_hasProbe", dmu_dsa_hasProbe, "dmu_dsa_hasProbe[ndmu]/O");
        tree_out->Branch("dmu_dsa_probeID", dmu_dsa_probeID, "dmu_dsa_probeID[ndmu]/I");
        tree_out->Branch("dmu_dsa_cosAlpha", dmu_dsa_cosAlpha, "dmu_dsa_cosAlpha[ndmu]/F");
        tree_out->Branch("dmu_dsa_isProbe", dmu_dsa_isProbe, "dmu_dsa_isProbe[ndmu]/O");

        tree_out->Branch("dmu_dgl_passTagID", dmu_dgl_passTagID, "dmu_dgl_passTagID[ndmu]/O");
        tree_out->Branch("dmu_dgl_hasProbe", dmu_dgl_hasProbe, "dmu_dgl_hasProbe[ndmu]/O");
        tree_out->Branch("dmu_dgl_probeID", dmu_dgl_probeID, "dmu_dgl_probeID[ndmu]/I");
        tree_out->Branch("dmu_dgl_cosAlpha", dmu_dgl_cosAlpha, "dmu_dgl_cosAlpha[ndmu]/F");
        tree_out->Branch("dmu_dgl_isProbe", dmu_dgl_isProbe, "dmu_dgl_isProbe[ndmu]/O");
    }
    // Gen Matching branches
    if (isMCSignal) {
        tree_out->Branch("dmu_hasGenMatch", dmu_hasGenMatch, "dmu_hasGenMatch[ndmu]/O");
        tree_out->Branch("dmu_genMatchedIndex", dmu_genMatchedIndex, "dmu_genMatchedIndex[ndmu]/I");
        tree_out->Branch("dmu_residual_r", dmu_residual_r, "dmu_residual_r[ndmu]/F");
        tree_out->Branch("dmu_residual_theta", dmu_residual_theta, "dmu_residual_theta[ndmu]/F");
        tree_out->Branch("dmu_residual_phi", dmu_residual_phi, "dmu_residual_phi[ndmu]/F");
        tree_out->Branch("dmu_residual_p_r", dmu_residual_p_r, "dmu_residual_p_r[ndmu]/F");
        tree_out->Branch("dmu_residual_p_theta", dmu_residual_p_theta,
                         "dmu_residual_p_theta[ndmu]/F");
        tree_out->Branch("dmu_residual_p_phi", dmu_residual_p_phi, "dmu_residual_p_phi[ndmu]/F");
        tree_out->Branch("dmu_residual_x", dmu_residual_x, "dmu_residual_x[ndmu]/F");
        tree_out->Branch("dmu_residual_y", dmu_residual_y, "dmu_residual_y[ndmu]/F");
        tree_out->Branch("dmu_residual_z", dmu_residual_z, "dmu_residual_z[ndmu]/F");
        tree_out->Branch("dmu_residual_p_x", dmu_residual_p_x, "dmu_residual_p_x[ndmu]/F");
        tree_out->Branch("dmu_residual_p_y", dmu_residual_p_y, "dmu_residual_p_y[ndmu]/F");
        tree_out->Branch("dmu_residual_p_z", dmu_residual_p_z, "dmu_residual_p_z[ndmu]/F");
        tree_out->Branch("dmu_propagationSurface", dmu_propagationSurface,
                         "dmu_propagationSurface[ndmu]/I");
    }
    // Trigger branches
    for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
        tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
    }
}

// endJob (After event loop has finished)
void ntuplizer::endJob() {
    std::cout << "End Job" << std::endl;
    // residualsfile.close();
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();
}

// fillDescriptions
void ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// Analyze (per event)
void ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    iEvent.getByToken(dglToken, dgls);
    iEvent.getByToken(dsaToken, dsas);
    iEvent.getByToken(dmuToken, dmuons);

    // const MagneticField* magField = nullptr;
    const Propagator* propagatorAlong = nullptr;
    const Propagator* propagatorOpposite = nullptr;
    if (isMCSignal) {
        iEvent.getByToken(prunedGenToken, prunedGen);
        // propagator = &iSetup.getData(thePropToken);
        propagatorAlong = &iSetup.getData(thePropAlongToken);
        propagatorOpposite = &iSetup.getData(thePropOppositeToken);

        if (!dynamic_cast<const SteppingHelixPropagator*>(propagatorAlong) ||
            !dynamic_cast<const SteppingHelixPropagator*>(propagatorOpposite)) {
            edm::LogWarning("BadConfig") << "@SUB=CosmicGenFilterHelix::getPropagator"
                                         << "Not a SteppingHelixPropagator!";
        }
        // magField_ = iSetup.getHandle(theMagFieldToken);
        // magField = magField_.product();
    }
    iEvent.getByToken(triggerBits_, triggerBits);

    // Count number of events read
    counts->Fill(0);

    // -> Event info
    event = iEvent.id().event();
    lumiBlock = iEvent.id().luminosityBlock();
    run = iEvent.id().run();

    // ----------------------------------
    // genParticles Collection
    // ----------------------------------
    Int_t nprugenmu = 0;
    std::vector<int> goodGenMuons_indices;
    int n_goodGenMuons = 0;
    if (isMCSignal) {
        for (unsigned int i = 0; i < prunedGen->size(); i++) {
            const reco::GenParticle& p(prunedGen->at(i));
            if (abs(p.pdgId()) == 13 && p.status() == 1) {  // Check if the particle is a muon
                // Check if the muon has a mother with pdgId 1023
                if (hasMotherWithPdgId(&p, 1023)) {
                    goodGenMuons_indices.push_back(i);
                    n_goodGenMuons++;
                }
            }
            nprugenmu++;
        }
    }

    // ----------------------------------
    // displacedMuons Collection
    // ----------------------------------
    ndmu = 0;
    std::vector<const reco::Muon*> matchedMuons;
    std::map<std::pair<int, int>, float> residuals;
    std::map<std::pair<int, int>, std::pair<GlobalTrajectoryParameters, bool>> matchResults;

    for (unsigned int i = 0; i < dmuons->size(); i++) {
        // std::cout << " - - ndmu: " << ndmu << std::endl;
        const reco::Muon& dmuon(dmuons->at(i));
        dmu_isDGL[ndmu] = dmuon.isGlobalMuon();
        dmu_isDSA[ndmu] = dmuon.isStandAloneMuon();
        dmu_isDTK[ndmu] = dmuon.isTrackerMuon();
        dmu_isMatchesValid[ndmu] = dmuon.isMatchesValid();
        dmu_numberOfMatches[ndmu] = dmuon.numberOfMatches();
        dmu_numberOfChambers[ndmu] = dmuon.numberOfChambers();
        dmu_numberOfChambersCSCorDT[ndmu] = dmuon.numberOfChambersCSCorDT();
        dmu_numberOfMatchedStations[ndmu] = dmuon.numberOfMatchedStations();
        dmu_numberOfMatchedRPCLayers[ndmu] = dmuon.numberOfMatchedRPCLayers();
        dmu_t0_InOut[ndmu] = dmuon.time().timeAtIpInOut;
        dmu_t0_OutIn[ndmu] = dmuon.time().timeAtIpOutIn;
        // Initialize the residuals to 9999
        dmu_residual_r[ndmu] = 9999;
        dmu_residual_theta[ndmu] = 9999;
        dmu_residual_phi[ndmu] = 9999;
        dmu_residual_p_r[ndmu] = 9999;
        dmu_residual_p_theta[ndmu] = 9999;
        dmu_residual_p_phi[ndmu] = 9999;
        dmu_residual_x[ndmu] = 9999;
        dmu_residual_y[ndmu] = 9999;
        dmu_residual_z[ndmu] = 9999;
        dmu_residual_p_x[ndmu] = 9999;
        dmu_residual_p_y[ndmu] = 9999;
        dmu_residual_p_z[ndmu] = 9999;
        dmu_propagationSurface[ndmu] = 0;
        dmu_hasGenMatch[ndmu] = false;
        dmu_genMatchedIndex[ndmu] = -1;

        // Access the DGL track associated to the displacedMuon
        // std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
        if (dmuon.isGlobalMuon()) {
            const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
            dmu_dgl_pt[ndmu] = globalTrack->pt();
            dmu_dgl_eta[ndmu] = globalTrack->eta();
            dmu_dgl_phi[ndmu] = globalTrack->phi();
            dmu_dgl_ptError[ndmu] = globalTrack->ptError();
            dmu_dgl_dxy[ndmu] = globalTrack->dxy();
            dmu_dgl_dz[ndmu] = globalTrack->dz();
            dmu_dgl_normalizedChi2[ndmu] = globalTrack->normalizedChi2();
            dmu_dgl_charge[ndmu] = globalTrack->charge();
            dmu_dgl_nMuonHits[ndmu] = globalTrack->hitPattern().numberOfMuonHits();
            dmu_dgl_nValidMuonHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonHits();
            dmu_dgl_nValidMuonDTHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonDTHits();
            dmu_dgl_nValidMuonCSCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonCSCHits();
            dmu_dgl_nValidMuonRPCHits[ndmu] = globalTrack->hitPattern().numberOfValidMuonRPCHits();
            dmu_dgl_nValidStripHits[ndmu] = globalTrack->hitPattern().numberOfValidStripHits();
            dmu_dgl_nhits[ndmu] = globalTrack->hitPattern().numberOfValidHits();
        } else {
            dmu_dgl_pt[ndmu] = 0;
            dmu_dgl_eta[ndmu] = 0;
            dmu_dgl_phi[ndmu] = 0;
            dmu_dgl_ptError[ndmu] = 0;
            dmu_dgl_dxy[ndmu] = 0;
            dmu_dgl_dz[ndmu] = 0;
            dmu_dgl_normalizedChi2[ndmu] = 0;
            dmu_dgl_charge[ndmu] = 0;
            dmu_dgl_nMuonHits[ndmu] = 0;
            dmu_dgl_nValidMuonHits[ndmu] = 0;
            dmu_dgl_nValidMuonDTHits[ndmu] = 0;
            dmu_dgl_nValidMuonCSCHits[ndmu] = 0;
            dmu_dgl_nValidMuonRPCHits[ndmu] = 0;
            dmu_dgl_nValidStripHits[ndmu] = 0;
            dmu_dgl_nhits[ndmu] = 0;
        }

        // Access the DSA track associated to the displacedMuon
        // std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
        if (dmuon.isStandAloneMuon()) {
            const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
            dmu_dsa_pt[ndmu] = outerTrack->pt();
            dmu_dsa_eta[ndmu] = outerTrack->eta();
            dmu_dsa_phi[ndmu] = outerTrack->phi();
            dmu_dsa_ptError[ndmu] = outerTrack->ptError();
            dmu_dsa_dxy[ndmu] = outerTrack->dxy();
            dmu_dsa_dz[ndmu] = outerTrack->dz();
            dmu_dsa_normalizedChi2[ndmu] = outerTrack->normalizedChi2();
            dmu_dsa_charge[ndmu] = outerTrack->charge();
            dmu_dsa_nMuonHits[ndmu] = outerTrack->hitPattern().numberOfMuonHits();
            dmu_dsa_nValidMuonHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonHits();
            dmu_dsa_nValidMuonDTHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonDTHits();
            dmu_dsa_nValidMuonCSCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonCSCHits();
            dmu_dsa_nValidMuonRPCHits[ndmu] = outerTrack->hitPattern().numberOfValidMuonRPCHits();
            dmu_dsa_nValidStripHits[ndmu] = outerTrack->hitPattern().numberOfValidStripHits();
            dmu_dsa_nhits[ndmu] = outerTrack->hitPattern().numberOfValidHits();
            dmu_dsa_dtStationsWithValidHits[ndmu] =
                outerTrack->hitPattern().dtStationsWithValidHits();
            dmu_dsa_cscStationsWithValidHits[ndmu] =
                outerTrack->hitPattern().cscStationsWithValidHits();
            if (isAOD) {
                // Number of DT+CSC segments
                unsigned int nsegments = 0;
                for (trackingRecHit_iterator hit = outerTrack->recHitsBegin();
                     hit != outerTrack->recHitsEnd(); ++hit) {
                    if (!(*hit)->isValid()) continue;
                    DetId id = (*hit)->geographicalId();
                    if (id.det() != DetId::Muon) continue;
                    if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC) {
                        nsegments++;
                    }
                }
                dmu_dsa_nsegments[ndmu] = nsegments;
            }
        } else {
            dmu_dsa_pt[ndmu] = 0;
            dmu_dsa_eta[ndmu] = 0;
            dmu_dsa_phi[ndmu] = 0;
            dmu_dsa_ptError[ndmu] = 0;
            dmu_dsa_dxy[ndmu] = 0;
            dmu_dsa_dz[ndmu] = 0;
            dmu_dsa_normalizedChi2[ndmu] = 0;
            dmu_dsa_charge[ndmu] = 0;
            dmu_dsa_nMuonHits[ndmu] = 0;
            dmu_dsa_nValidMuonHits[ndmu] = 0;
            dmu_dsa_nValidMuonDTHits[ndmu] = 0;
            dmu_dsa_nValidMuonCSCHits[ndmu] = 0;
            dmu_dsa_nValidMuonRPCHits[ndmu] = 0;
            dmu_dsa_nValidStripHits[ndmu] = 0;
            dmu_dsa_nhits[ndmu] = 0;
            dmu_dsa_dtStationsWithValidHits[ndmu] = 0;
            dmu_dsa_cscStationsWithValidHits[ndmu] = 0;
            dmu_dsa_nsegments[ndmu] = 0;
        }
        // ----------------------------------
        // Gen Matching part
        // ----------------------------------
        if (isMCSignal) {
            const reco::Track* candidateTrack = nullptr;
            if (dmuon.isGlobalMuon()) {
                candidateTrack = (dmuon.combinedMuon()).get();
            } else if (dmuon.isStandAloneMuon()) {
                candidateTrack = (dmuon.standAloneMuon()).get();
            } else if (dmuon.isTrackerMuon()) {
                candidateTrack = (dmuon.innerTrack()).get();
            }

            TVector3 candidate_vertex = TVector3();
            candidate_vertex.SetXYZ(candidateTrack->vx(), candidateTrack->vy(),
                                    candidateTrack->vz());
            dmu_vertex_r[ndmu] = sqrt(candidate_vertex.X() * candidate_vertex.X() +
                                      candidate_vertex.Y() * candidate_vertex.Y());
            dmu_vertex_z[ndmu] = candidate_vertex.Z();
        }
        ndmu++;

        // ----------------------------------
        // Tag and probe code
        // ----------------------------------
        if (isCosmics) {
            ndmu = 0;
            for (unsigned int i = 0; i < dmuons->size(); i++) {
                const reco::Muon& dmuon(dmuons->at(i));
                // Access the DGL track associated to the displacedMuon
                // std::cout << "isGlobalMuon: " << dmuon.isGlobalMuon() << std::endl;
                if (dmuon.isGlobalMuon()) {
                    const reco::Track* globalTrack = (dmuon.combinedMuon()).get();
                    // Fill tag and probe variables
                    //   First, reset the variables
                    dmu_dgl_passTagID[ndmu] = false;
                    dmu_dgl_hasProbe[ndmu] = false;
                    dmu_dgl_probeID[ndmu] = 0;
                    dmu_dgl_cosAlpha[ndmu] = 0.;
                    // Check if muon passes tag ID
                    dmu_dgl_passTagID[ndmu] = passTagID(globalTrack, "DGL");
                    if (dmu_dgl_passTagID[ndmu]) {
                        // Search probe
                        TVector3 v_tag = TVector3();
                        v_tag.SetPtEtaPhi(globalTrack->pt(), globalTrack->eta(),
                                          globalTrack->phi());
                        const reco::Muon* muonProbeTemp =
                            nullptr;  // pointer for temporal probe (initialized to nullptr)
                        for (unsigned int j = 0; j < dmuons->size();
                             j++) {  // Loop over the rest of the muons
                            if (i == j) {
                                continue;
                            }
                            const reco::Muon& muonProbeCandidate(dmuons->at(j));
                            if (!muonProbeCandidate.isGlobalMuon()) {
                                continue;
                            }  // Get only dgls
                            const reco::Track* trackProbeCandidate =
                                (muonProbeCandidate.combinedMuon()).get();
                            if (passProbeID(trackProbeCandidate, v_tag, "DGL")) {
                                TVector3 v_probe = TVector3();
                                v_probe.SetPtEtaPhi(trackProbeCandidate->pt(),
                                                    trackProbeCandidate->eta(),
                                                    trackProbeCandidate->phi());
                                if (!dmu_dgl_hasProbe[ndmu]) {
                                    dmu_dgl_hasProbe[ndmu] = true;
                                    muonProbeTemp = &(dmuons->at(j));
                                    dmu_dgl_probeID[ndmu] = j;
                                    dmu_dgl_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
                                } else {
                                    const reco::Track* trackProbeTemp =
                                        (muonProbeTemp->combinedMuon()).get();
                                    if (trackProbeCandidate->pt() > trackProbeTemp->pt()) {
                                        muonProbeTemp = &(dmuons->at(j));
                                        dmu_dgl_probeID[ndmu] = j;
                                        dmu_dgl_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
                                    } else {
                                        std::cout << ">> Probe candidate " << j
                                                  << " has lower pt than " << dmu_dgl_probeID[ndmu]
                                                  << std::endl;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    dmu_dgl_passTagID[ndmu] = false;
                    dmu_dgl_hasProbe[ndmu] = false;
                    dmu_dgl_probeID[ndmu] = 0;
                    dmu_dgl_cosAlpha[ndmu] = 0.;
                }
                // Access the DSA track associated to the displacedMuon
                // std::cout << "isStandAloneMuon: " << dmuon.isStandAloneMuon() << std::endl;
                if (dmuon.isStandAloneMuon()) {
                    const reco::Track* outerTrack = (dmuon.standAloneMuon()).get();
                    // Fill tag and probe variables
                    //   First, reset the variables
                    dmu_dsa_passTagID[ndmu] = false;
                    dmu_dsa_hasProbe[ndmu] = false;
                    dmu_dsa_probeID[ndmu] = 0;
                    dmu_dsa_cosAlpha[ndmu] = 0.;
                    // Check if muon passes tag ID
                    dmu_dsa_passTagID[ndmu] = passTagID(outerTrack, "DSA");
                    if (dmu_dsa_passTagID[ndmu]) {
                        // Search probe
                        TVector3 v_tag = TVector3();
                        v_tag.SetPtEtaPhi(outerTrack->pt(), outerTrack->eta(), outerTrack->phi());
                        const reco::Muon* muonProbeTemp =
                            nullptr;  // pointer for temporal probe (initialized to nullptr)
                        for (unsigned int j = 0; j < dmuons->size();
                             j++) {  // Loop over the rest of the muons
                            if (i == j) {
                                continue;
                            }
                            const reco::Muon& muonProbeCandidate(dmuons->at(j));
                            if (!muonProbeCandidate.isStandAloneMuon()) {
                                continue;
                            }  // Get only dsas
                            const reco::Track* trackProbeCandidate =
                                (muonProbeCandidate.standAloneMuon()).get();
                            if (passProbeID(trackProbeCandidate, v_tag, "DSA")) {
                                TVector3 v_probe = TVector3();
                                v_probe.SetPtEtaPhi(trackProbeCandidate->pt(),
                                                    trackProbeCandidate->eta(),
                                                    trackProbeCandidate->phi());
                                if (!dmu_dsa_hasProbe[ndmu]) {
                                    dmu_dsa_hasProbe[ndmu] = true;
                                    muonProbeTemp = &(dmuons->at(j));
                                    dmu_dsa_probeID[ndmu] = j;
                                    dmu_dsa_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
                                } else {
                                    const reco::Track* trackProbeTemp =
                                        (muonProbeTemp->standAloneMuon()).get();
                                    if (trackProbeCandidate->pt() > trackProbeTemp->pt()) {
                                        muonProbeTemp = &(dmuons->at(j));
                                        dmu_dsa_probeID[ndmu] = j;
                                        dmu_dsa_cosAlpha[ndmu] = cos(v_tag.Angle(v_probe));
                                    } else {
                                        std::cout << ">> Probe candidate " << j
                                                  << " has lower pt than " << dmu_dsa_probeID[ndmu]
                                                  << std::endl;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    dmu_dsa_passTagID[ndmu] = false;
                    dmu_dsa_hasProbe[ndmu] = false;
                    dmu_dsa_probeID[ndmu] = 0;
                    dmu_dsa_cosAlpha[ndmu] = 0.;
                }
                ndmu++;
            }
            // After the tag and probe identification loop, identify the probe muons
            // Set all values of dmu_*_isProbe to false
            for (int i = 0; i < ndmu; i++) {
                dmu_dgl_isProbe[i] = false;
                dmu_dsa_isProbe[i] = false;
            }
            ndmu = 0;
            for (unsigned int i = 0; i < dmuons->size(); i++) {
                if (dmu_dgl_hasProbe[ndmu]) {
                    dmu_dgl_isProbe[dmu_dgl_probeID[ndmu]] = true;
                }
                if (dmu_dsa_hasProbe[ndmu]) {
                    dmu_dsa_isProbe[dmu_dsa_probeID[ndmu]] = true;
                }
                ndmu++;
            }
        }

        // Check if trigger fired:
        const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);
        unsigned int ipath = 0;
        for (auto path : HLTPaths_) {
            std::string path_v = path + "_v";
            // std::cout << path << "\t" << std::endl;
            bool fired = false;
            for (unsigned int itrg = 0; itrg < triggerBits->size(); ++itrg) {
                TString TrigPath = names.triggerName(itrg);
                if (!triggerBits->accept(itrg)) continue;
                if (!TrigPath.Contains(path_v)) {
                    continue;
                }
                fired = true;
            }
            triggerPass[ipath] = fired;
            ipath++;
        }
    }

    //-> Fill tree
    tree_out->Fill();
}
DEFINE_FWK_MODULE(ntuplizer);
