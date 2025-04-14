#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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
#include "Math/Error.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMatrixF.h"
#include "TTree.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryParametrization/interface/CartesianTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "genmatching_utilities.h"
#include "propagate_utilities.h"
#include "propagation_definitions.h"

class ntuplizer_gens : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
    explicit ntuplizer_gens(const edm::ParameterSet&);
    ~ntuplizer_gens();

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
    Int_t ngenmu = 0;
    Int_t genmu_isDSA[200] = {0};
    Int_t genmu_isDGL[200] = {0};
    Int_t genmu_isDTK[200] = {0};

    Float_t genmu_gen_initial_r[200] = {9999.};
    Float_t genmu_gen_initial_z[200] = {9999.};
    Float_t genmu_gen_initial_theta[200] = {9999.};
    Float_t genmu_gen_initial_phi[200] = {9999.};
    Float_t genmu_gen_final_r[200] = {9999.};
    Float_t genmu_gen_final_theta[200] = {9999.};
    Float_t genmu_gen_final_phi[200] = {9999.};
    Float_t genmu_gen_initial_p_r[200] = {9999.};
    Float_t genmu_gen_initial_p_theta[200] = {9999.};
    Float_t genmu_gen_initial_p_phi[200] = {9999.};
    Float_t genmu_gen_final_p_r[200] = {9999.};
    Float_t genmu_gen_final_p_theta[200] = {9999.};
    Float_t genmu_gen_final_p_phi[200] = {9999.};
    Int_t genmu_propagationSurface[200] = {-1};

    Float_t genmu_dsa_pt[200] = {0.};
    Float_t genmu_dsa_eta[200] = {0.};
    Float_t genmu_dsa_phi[200] = {0.};
    Float_t genmu_dsa_dxy[200] = {0.};
    Float_t genmu_dsa_dz[200] = {0.};
    Float_t genmu_dsa_normalizedChi2[200] = {0.};
    Float_t genmu_dsa_charge[200] = {0.};

    //
    // --- Output
    //
    std::string output_filename;
    TH1F* counts;
    TFile* file_out;
    TTree* tree_out;
};

// Constructor
ntuplizer_gens::ntuplizer_gens(const edm::ParameterSet& iConfig) {
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
ntuplizer_gens::~ntuplizer_gens() {}

// beginJob (Before first event)
void ntuplizer_gens::beginJob() {
    std::cout << "Begin Job" << std::endl;

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

    tree_out->Branch("ngenmu", &ngenmu, "ngenmu/I");
    tree_out->Branch("genmu_isDSA", genmu_isDSA, "genmu_isDSA[ngenmu]/I");
    tree_out->Branch("genmu_isDGL", genmu_isDGL, "genmu_isDGL[ngenmu]/I");
    tree_out->Branch("genmu_isDTK", genmu_isDTK, "genmu_isDTK[ngenmu]/I");
    tree_out->Branch("genmu_gen_initial_r", genmu_gen_initial_r, "genmu_gen_initial_r[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_z", genmu_gen_initial_z, "genmu_gen_initial_z[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_theta", genmu_gen_initial_theta,
                     "genmu_gen_initial_theta[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_phi", genmu_gen_initial_phi,
                     "genmu_gen_initial_phi[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_r", genmu_gen_final_r, "genmu_gen_final_r[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_theta", genmu_gen_final_theta,
                     "genmu_gen_final_theta[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_phi", genmu_gen_final_phi, "genmu_gen_final_phi[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_p_r", genmu_gen_initial_p_r,
                     "genmu_gen_initial_p_r[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_p_theta", genmu_gen_initial_p_theta,
                     "genmu_gen_initial_p_theta[ngenmu]/F");
    tree_out->Branch("genmu_gen_initial_p_phi", genmu_gen_initial_p_phi,
                     "genmu_gen_initial_p_phi[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_p_r", genmu_gen_final_p_r, "genmu_gen_final_p_r[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_p_theta", genmu_gen_final_p_theta,
                     "genmu_gen_final_p_theta[ngenmu]/F");
    tree_out->Branch("genmu_gen_final_p_phi", genmu_gen_final_p_phi,
                     "genmu_gen_final_p_phi[ngenmu]/F");
    tree_out->Branch("genmu_propagationSurface", genmu_propagationSurface,
                     "genmu_propagationSurface[ngenmu]/I");
    tree_out->Branch("genmu_dsa_pt", genmu_dsa_pt, "genmu_dsa_pt[ngenmu]/F");
    tree_out->Branch("genmu_dsa_eta", genmu_dsa_eta, "genmu_dsa_eta[ngenmu]/F");
    tree_out->Branch("genmu_dsa_phi", genmu_dsa_phi, "genmu_dsa_phi[ngenmu]/F");
    tree_out->Branch("genmu_dsa_charge", genmu_dsa_charge, "genmu_dsa_charge[ngenmu]/F");
    // Trigger branches
    for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
        tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
    }
}

// endJob (After event loop has finished)
void ntuplizer_gens::endJob() {
    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();
}

// fillDescriptions
void ntuplizer_gens::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// Analyze (per event)
void ntuplizer_gens::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    iEvent.getByToken(dglToken, dgls);
    iEvent.getByToken(dsaToken, dsas);
    iEvent.getByToken(dmuToken, dmuons);

    const MagneticField* magField = nullptr;
    const Propagator* propagatorAlong = nullptr;
    const Propagator* propagatorOpposite = nullptr;
    if (isMCSignal) {
        iEvent.getByToken(prunedGenToken, prunedGen);
        propagatorAlong = &iSetup.getData(thePropAlongToken);
        propagatorOpposite = &iSetup.getData(thePropOppositeToken);

        if (!dynamic_cast<const SteppingHelixPropagator*>(propagatorAlong) ||
            !dynamic_cast<const SteppingHelixPropagator*>(propagatorOpposite)) {
            edm::LogWarning("BadConfig") << "@SUB=CosmicGenFilterHelix::getPropagator"
                                         << "Not a SteppingHelixPropagator!";
        }
        magField_ = iSetup.getHandle(theMagFieldToken);
        magField = magField_.product();
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
    ngenmu = 0;
    if (isMCSignal) {
        for (unsigned int i = 0; i < prunedGen->size(); i++) {
            const reco::GenParticle& genParticle(prunedGen->at(i));
            if (abs(genParticle.pdgId()) == 13 &&
                genParticle.status() == 1) {  // Check if the particle is a muon
                // Check if the muon has a mother with pdgId 1023
                if (hasMotherWithPdgId(&genParticle, 1023)) {
                    // Initial state information in genFTS
                    GlobalTrajectoryParameters genFinalParams = GlobalTrajectoryParameters();
                    GlobalVector genMomentum(genParticle.px(), genParticle.py(), genParticle.pz());
                    // Gen vertices are in mm, make it cm
                    GlobalPoint genVertex(genParticle.vx() * 0.1, genParticle.vy() * 0.1,
                                          genParticle.vz() * 0.1);
                    int genCharge = genParticle.charge();
                    // Ntuple variables
                    genmu_gen_initial_r[ngenmu] = genVertex.perp();
                    genmu_gen_initial_z[ngenmu] = genVertex.z();
                    genmu_gen_initial_theta[ngenmu] = genVertex.theta();
                    genmu_gen_initial_phi[ngenmu] = genVertex.phi();
                    genmu_gen_initial_p_r[ngenmu] = genMomentum.perp();
                    genmu_gen_initial_p_theta[ngenmu] = genMomentum.theta();
                    genmu_gen_initial_p_phi[ngenmu] = genMomentum.phi();
                    genmu_gen_final_r[ngenmu] = 9999;
                    genmu_gen_final_theta[ngenmu] = 9999;
                    genmu_gen_final_phi[ngenmu] = 9999;
                    genmu_gen_final_p_r[ngenmu] = 9999;
                    genmu_gen_final_p_theta[ngenmu] = 9999;
                    genmu_gen_final_p_phi[ngenmu] = 9999;
                    genmu_propagationSurface[ngenmu] = -1;
                    genmu_dsa_pt[ngenmu] = genMomentum.perp();
                    genmu_dsa_eta[ngenmu] = genMomentum.eta();
                    genmu_dsa_phi[ngenmu] = genMomentum.phi();
                    genmu_dsa_charge[ngenmu] = genCharge;
                    genmu_isDSA[ngenmu] = genParticle.isStandAloneMuon();
                    genmu_isDGL[ngenmu] = genParticle.isGlobalMuon();
                    genmu_isDTK[ngenmu] = genParticle.isTrackerMuon();
                    FreeTrajectoryState genFTS(genVertex, genMomentum, genCharge, magField);
                    // final state information in genFinalParams
                    PropagationSurface genSurface = findAndPropagateToOptimalSurface(
                        genFTS, genFinalParams, magField, propagatorAlong, propagatorOpposite);

                    bool goodMatch = (static_cast<int>(genSurface.genMatchResult) > 0);
                    if (goodMatch) {
                        genFinalParams = GlobalTrajectoryParameters(genFinalParams.position(),
                                                                    genFinalParams.momentum(),
                                                                    genCharge, magField);
                        genmu_gen_final_r[ngenmu] = genFinalParams.position().perp();
                        genmu_gen_final_theta[ngenmu] = genFinalParams.position().theta();
                        genmu_gen_final_phi[ngenmu] = genFinalParams.position().phi();
                        genmu_gen_final_p_r[ngenmu] = genFinalParams.momentum().perp();
                        genmu_gen_final_p_theta[ngenmu] = genFinalParams.momentum().theta();
                        genmu_gen_final_p_phi[ngenmu] = genFinalParams.momentum().phi();
                    }
                    genmu_propagationSurface[ngenmu] =
                        static_cast<Int_t>(genSurface.genMatchResult);
                    ngenmu++;
                }
            }
        }
    }

    //-> Fill tree
    tree_out->Fill();
}
DEFINE_FWK_MODULE(ntuplizer_gens);