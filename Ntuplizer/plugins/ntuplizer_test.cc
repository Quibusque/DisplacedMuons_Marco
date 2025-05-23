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

class ntuplizer_test : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
    explicit ntuplizer_test(const edm::ParameterSet&);
    ~ntuplizer_test();

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

    // Variables for gen matching
    Int_t ngenmu = 0;
    Int_t genmu_kindOfMatching[10] = {-1};
    Int_t genmu_propagationSurface[10] = {-1};
    Int_t genmu_charge[10] = {0};
    Float_t genmu_pt[10] = {0};
    Float_t genmu_initial_r[10] = {9999.};
    Float_t genmu_initial_theta[10] = {9999.};
    Float_t genmu_initial_phi[10] = {9999.};
    Float_t genmu_final_r[10] = {9999.};
    Float_t genmu_final_theta[10] = {9999.};
    Float_t genmu_final_phi[10] = {9999.};
    Float_t genmu_initial_p_r[10] = {9999.};
    Float_t genmu_initial_p_theta[10] = {9999.};
    Float_t genmu_initial_p_phi[10] = {9999.};
    Float_t genmu_final_p_r[10] = {9999.};
    Float_t genmu_final_p_theta[10] = {9999.};
    Float_t genmu_final_p_phi[10] = {9999.};
    bool dmu_hasGenMatch[200] = {false};
    Int_t dmu_genMatchedIndex[200] = {0};
    Float_t dmu_reco_initial_r[200] = {9999.};
    Float_t dmu_reco_initial_theta[200] = {9999.};
    Float_t dmu_reco_initial_phi[200] = {9999.};
    Float_t dmu_reco_final_r[200] = {9999.};
    Float_t dmu_reco_final_theta[200] = {9999.};
    Float_t dmu_reco_final_phi[200] = {9999.};
    Float_t dmu_reco_initial_p_r[200] = {9999.};
    Float_t dmu_reco_initial_p_theta[200] = {9999.};
    Float_t dmu_reco_initial_p_phi[200] = {9999.};
    Float_t dmu_reco_final_p_r[200] = {9999.};
    Float_t dmu_reco_final_p_theta[200] = {9999.};
    Float_t dmu_reco_final_p_phi[200] = {9999.};
    Float_t dmu_reco_final_x_err[200] = {9999.};
    Float_t dmu_reco_final_y_err[200] = {9999.};
    Float_t dmu_reco_final_z_err[200] = {9999.};
    Float_t dmu_reco_final_px_err[200] = {9999.};
    Float_t dmu_reco_final_py_err[200] = {9999.};
    Float_t dmu_reco_final_pz_err[200] = {9999.};

    Int_t dmu_propagationSurface[200] = {-1};
    Float_t dmu_chi2[200] = {9999.};
    Float_t dmu_chi2_x[200] = {9999.};
    Float_t dmu_chi2_y[200] = {9999.};
    Float_t dmu_chi2_z[200] = {9999.};
    Float_t dmu_chi2_px[200] = {9999.};
    Float_t dmu_chi2_py[200] = {9999.};
    Float_t dmu_chi2_pz[200] = {9999.};

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

    //
    // --- Output
    //
    std::string output_filename;
    TH1F* counts;
    TFile* file_out;
    TTree* tree_out;
};

// Constructor
ntuplizer_test::ntuplizer_test(const edm::ParameterSet& iConfig) {
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
ntuplizer_test::~ntuplizer_test() {}

// beginJob (Before first event)
void ntuplizer_test::beginJob() {
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

    // Gen Matching branches
    if (isMCSignal) {
        tree_out->Branch("ngenmu", &ngenmu, "ngenmu/I");
        tree_out->Branch("genmu_kindOfMatching", genmu_kindOfMatching,
                    "genmu_kindOfMatching[ngenmu]/I");
        tree_out->Branch("genmu_propagationSurface", genmu_propagationSurface,
                         "genmu_propagationSurface[ngenmu]/I");
        tree_out->Branch("genmu_charge", genmu_charge, "genmu_charge[ngenmu]/I");
        tree_out->Branch("genmu_pt", genmu_pt, "genmu_pt[ngenmu]/F");
        tree_out->Branch("dmu_hasGenMatch", dmu_hasGenMatch, "dmu_hasGenMatch[ndmu]/O");
        tree_out->Branch("dmu_genMatchedIndex", dmu_genMatchedIndex, "dmu_genMatchedIndex[ndmu]/I");
        tree_out->Branch("dmu_reco_initial_r", dmu_reco_initial_r, "dmu_reco_initial_r[ndmu]/F");
        tree_out->Branch("dmu_reco_initial_theta", dmu_reco_initial_theta,
                         "dmu_reco_initial_theta[ndmu]/F");
        tree_out->Branch("dmu_reco_initial_phi", dmu_reco_initial_phi,
                         "dmu_reco_initial_phi[ndmu]/F");
        tree_out->Branch("dmu_reco_final_r", dmu_reco_final_r, "dmu_reco_final_r[ndmu]/F");
        tree_out->Branch("dmu_reco_final_theta", dmu_reco_final_theta,
                         "dmu_reco_final_theta[ndmu]/F");
        tree_out->Branch("dmu_reco_final_phi", dmu_reco_final_phi, "dmu_reco_final_phi[ndmu]/F");
        tree_out->Branch("dmu_reco_initial_p_r", dmu_reco_initial_p_r,
                         "dmu_reco_initial_p_r[ndmu]/F");
        tree_out->Branch("dmu_reco_initial_p_theta", dmu_reco_initial_p_theta,
                         "dmu_reco_initial_p_theta[ndmu]/F");
        tree_out->Branch("dmu_reco_initial_p_phi", dmu_reco_initial_p_phi,
                         "dmu_reco_initial_p_phi[ndmu]/F");
        tree_out->Branch("dmu_reco_final_p_r", dmu_reco_final_p_r, "dmu_reco_final_p_r[ndmu]/F");
        tree_out->Branch("dmu_reco_final_p_theta", dmu_reco_final_p_theta,
                         "dmu_reco_final_p_theta[ndmu]/F");
        tree_out->Branch("dmu_reco_final_p_phi", dmu_reco_final_p_phi,
                         "dmu_reco_final_p_phi[ndmu]/F");
        tree_out->Branch("genmu_initial_r", genmu_initial_r, "genmu_initial_r[ngenmu]/F");
        tree_out->Branch("genmu_initial_theta", genmu_initial_theta,
                         "genmu_initial_theta[ngenmu]/F");
        tree_out->Branch("genmu_initial_phi", genmu_initial_phi, "genmu_initial_phi[ngenmu]/F");
        tree_out->Branch("genmu_final_r", genmu_final_r, "genmu_final_r[ngenmu]/F");
        tree_out->Branch("genmu_final_theta", genmu_final_theta, "genmu_final_theta[ngenmu]/F");
        tree_out->Branch("genmu_final_phi", genmu_final_phi, "genmu_final_phi[ngenmu]/F");
        tree_out->Branch("genmu_initial_p_r", genmu_initial_p_r, "genmu_initial_p_r[ngenmu]/F");
        tree_out->Branch("genmu_initial_p_theta", genmu_initial_p_theta,
                         "genmu_initial_p_theta[ngenmu]/F");
        tree_out->Branch("genmu_initial_p_phi", genmu_initial_p_phi,
                         "genmu_initial_p_phi[ngenmu]/F");
        tree_out->Branch("genmu_final_p_r", genmu_final_p_r, "genmu_final_p_r[ngenmu]/F");
        tree_out->Branch("genmu_final_p_theta", genmu_final_p_theta,
                         "genmu_final_p_theta[ngenmu]/F");
        tree_out->Branch("genmu_final_p_phi", genmu_final_p_phi, "genmu_final_p_phi[ngenmu]/F");
        tree_out->Branch("dmu_propagationSurface", dmu_propagationSurface,
                         "dmu_propagationSurface[ndmu]/I");
        tree_out->Branch("dmu_chi2", dmu_chi2, "dmu_chi2[ndmu]/F");
        tree_out->Branch("dmu_chi2_x", dmu_chi2_x, "dmu_chi2_x[ndmu]/F");
        tree_out->Branch("dmu_chi2_y", dmu_chi2_y, "dmu_chi2_y[ndmu]/F");
        tree_out->Branch("dmu_chi2_z", dmu_chi2_z, "dmu_chi2_z[ndmu]/F");
        tree_out->Branch("dmu_chi2_px", dmu_chi2_px, "dmu_chi2_px[ndmu]/F");
        tree_out->Branch("dmu_chi2_py", dmu_chi2_py, "dmu_chi2_py[ndmu]/F");
        tree_out->Branch("dmu_chi2_pz", dmu_chi2_pz, "dmu_chi2_pz[ndmu]/F");
        tree_out->Branch("dmu_reco_final_x_err", dmu_reco_final_x_err,
                         "dmu_reco_final_x_err[ndmu]/F");
        tree_out->Branch("dmu_reco_final_y_err", dmu_reco_final_y_err,
                         "dmu_reco_final_y_err[ndmu]/F");
        tree_out->Branch("dmu_reco_final_z_err", dmu_reco_final_z_err,
                         "dmu_reco_final_z_err[ndmu]/F");
        tree_out->Branch("dmu_reco_final_px_err", dmu_reco_final_px_err,
                         "dmu_reco_final_px_err[ndmu]/F");
        tree_out->Branch("dmu_reco_final_py_err", dmu_reco_final_py_err,
                         "dmu_reco_final_py_err[ndmu]/F");
        tree_out->Branch("dmu_reco_final_pz_err", dmu_reco_final_pz_err,
                         "dmu_reco_final_pz_err[ndmu]/F");
    }
    // Trigger branches
    for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
        tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
    }
}

// endJob (After event loop has finished)
void ntuplizer_test::endJob() {
    std::cout << "End Job" << std::endl;
    file_out->cd();
    tree_out->Write();
    counts->Write();
    file_out->Close();
}

// fillDescriptions
void ntuplizer_test::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// Analyze (per event)
void ntuplizer_test::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    std::pair<PropagationSurface, GlobalTrajectoryParameters> genPropagationResults[10];
    if (isMCSignal) {
        for (unsigned int i = 0; i < prunedGen->size(); i++) {
            const reco::GenParticle& genParticle(prunedGen->at(i));
            if (abs(genParticle.pdgId()) == 13 &&
                genParticle.status() == 1) {  // Check if the particle is a muon
                // Check if the muon has a mother with pdgId 1023
                if (hasMotherWithPdgId(&genParticle, 1023)) {

                    // ---------------------------------------------------------
                    // Propagate the gens to the surface for later gen matching
                    // ---------------------------------------------------------

                    // Initial state information in genFTS, final state information in
                    // genFinalParams
                    GlobalTrajectoryParameters genFinalParams = GlobalTrajectoryParameters();
                    GlobalVector genMomentum(genParticle.px(), genParticle.py(), genParticle.pz());
                    // Gen vertices are in cm and that's okay
                    GlobalPoint genVertex(genParticle.vx(), genParticle.vy(), genParticle.vz());
                    int genCharge = genParticle.charge();
                    FreeTrajectoryState genFTS(genVertex, genMomentum, genCharge, magField);

                    PropagationSurface genSurface = findAndPropagateToOptimalSurface(
                        genFTS, genFinalParams, magField, propagatorAlong, propagatorOpposite);

                    bool goodMatch = (static_cast<int>(genSurface.genMatchResult) > 0);
                    if (goodMatch) {
                        genFinalParams = GlobalTrajectoryParameters(genFinalParams.position(),
                                                                    genFinalParams.momentum(),
                                                                    genCharge, magField);
                    }
                    genPropagationResults[ngenmu] = std::make_pair(genSurface, genFinalParams);

                    if (genSurface == PropagationConstants::GEN_OUTSIDE_CMS) {
                        genmu_kindOfMatching[ngenmu] =
                            static_cast<int>(genSurface.genMatchResult);
                    } else {
                        genmu_kindOfMatching[ngenmu] = -1;
                    }
                    genmu_propagationSurface[ngenmu] = static_cast<int>(genSurface.genMatchResult);
                    genmu_pt[ngenmu] = genParticle.pt();
                    genmu_charge[ngenmu] = genCharge;
                    genmu_initial_r[ngenmu] = genVertex.perp();
                    genmu_initial_theta[ngenmu] = genVertex.theta();
                    genmu_initial_phi[ngenmu] = genVertex.phi();
                    genmu_initial_p_r[ngenmu] = genMomentum.perp();
                    genmu_initial_p_theta[ngenmu] = genMomentum.theta();
                    genmu_initial_p_phi[ngenmu] = genMomentum.phi();

                    ngenmu++;
                }
            }
        }
    }

    // ----------------------------------
    // displacedMuons Collection
    // ----------------------------------
    ndmu = 0;
    std::map<std::pair<int, int>, std::pair<GlobalTrajectoryParameters, GlobalTrajectoryParameters>>
        propagatedTrajectories;
    TMatrixF chi2Matrix = TMatrixF(dmuons->size(), ngenmu);
    std::map<std::pair<int, int>, AlgebraicVector6> error_vectors;
    std::map<std::pair<int, int>, AlgebraicVector6> chi2_vectors;
    std::map<std::pair<int, int>, GenMatchResults> matchResults;
    for (unsigned int i = 0; i < dmuons->size(); i++) {
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
        dmu_hasGenMatch[ndmu] = false;
        dmu_propagationSurface[ndmu] = -1;
        dmu_chi2[ndmu] = 9999;
        dmu_chi2_x[ndmu] = 9999;
        dmu_chi2_y[ndmu] = 9999;
        dmu_chi2_z[ndmu] = 9999;
        dmu_chi2_px[ndmu] = 9999;
        dmu_chi2_py[ndmu] = 9999;
        dmu_chi2_pz[ndmu] = 9999;
        dmu_reco_final_x_err[ndmu] = 9999;
        dmu_reco_final_y_err[ndmu] = 9999;
        dmu_reco_final_z_err[ndmu] = 9999;
        dmu_reco_final_px_err[ndmu] = 9999;
        dmu_reco_final_py_err[ndmu] = 9999;
        dmu_reco_final_pz_err[ndmu] = 9999;

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
            TVector3 dsa_vertex = TVector3();
            TVector3 dsa_momentum = TVector3();
            dsa_vertex.SetXYZ(outerTrack->vx(), outerTrack->vy(), outerTrack->vz());
            dsa_momentum.SetXYZ(outerTrack->px(), outerTrack->py(), outerTrack->pz());
            // initial values for the gen matching
            dmu_reco_initial_r[ndmu] = dsa_vertex.Perp();
            dmu_reco_initial_theta[ndmu] = dsa_vertex.Theta();
            dmu_reco_initial_phi[ndmu] = dsa_vertex.Phi();
            dmu_reco_initial_p_r[ndmu] = dsa_momentum.Perp();
            dmu_reco_initial_p_theta[ndmu] = dsa_momentum.Theta();
            dmu_reco_initial_p_phi[ndmu] = dsa_momentum.Phi();
            dmu_reco_final_r[ndmu] = 9999;
            dmu_reco_final_theta[ndmu] = 9999;
            dmu_reco_final_phi[ndmu] = 9999;
            dmu_reco_final_p_r[ndmu] = 9999;
            dmu_reco_final_p_theta[ndmu] = 9999;
            dmu_reco_final_p_phi[ndmu] = 9999;


            for (int j = 0; j < ngenmu; j++) {
                chi2Matrix(ndmu, j) = 9999;
                matchResults[{ndmu, j}] = GenMatchResults::NONE;
                chi2_vectors[{ndmu, j}] = AlgebraicVector6(9999, 9999, 9999, 9999, 9999, 9999);
                error_vectors[{ndmu, j}] = AlgebraicVector6(9999, 9999, 9999, 9999, 9999, 9999);
            }
        } else {
            continue;
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
            } else {
                continue;
            }

            for (Int_t j = 0; j < ngenmu; j++) {
                PropagationSurface genSurface = genPropagationResults[j].first;
                GlobalTrajectoryParameters genFinalParams = genPropagationResults[j].second;
                GlobalTrajectoryParameters recoFinalParams;
                CartesianTrajectoryError recoError;

                Float_t deltaR_threshold = 1000;

                GenMatchResults matchResult = matchRecoTrackToGenSurface(
                    genSurface, candidateTrack, magField, propagatorAlong, propagatorOpposite,
                    genFinalParams, recoFinalParams, recoError, deltaR_threshold);

                propagatedTrajectories[{ndmu, j}] = {genFinalParams, recoFinalParams};
                matchResults[{ndmu, j}] = matchResult;
                bool goodMatch = (static_cast<int>(matchResult) > 0);
                AlgebraicVector6 chi2_vector =
                    calculateChi2Vector(genFinalParams, recoFinalParams, recoError);

                AlgebraicVector6 error_vector = AlgebraicVector6(
                    recoError.matrix()(0, 0), recoError.matrix()(1, 1), recoError.matrix()(2, 2),
                    recoError.matrix()(3, 3), recoError.matrix()(4, 4), recoError.matrix()(5, 5));
                error_vectors[{ndmu, j}] = error_vector;
                Float_t chi2 = 0;
                for (int k = 0; k < 6; k++) {
                    chi2 += chi2_vector[k];
                    // fill chi2_vectors
                    chi2_vectors[{ndmu, j}][k] = chi2_vector[k];
                }
                chi2Matrix(ndmu, j) = (goodMatch) ? chi2 : 9999;
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
        ndmu++;
    }
    // end of loop over reco displacedMuons

    //  -----------------------------------------------------
    //  Gen Matching continued MUST BE OUTSIDE LOOP ON RECOS
    //  -----------------------------------------------------

    // dmuons->size() is larger than ndmu
    // resize the chi2Matrix to ndmu,ngenmu size
    chi2Matrix.ResizeTo(ndmu, ngenmu);
    TMatrixF boolMatrix = TMatrixF(200, 200);
    markUniqueBestMatches(chi2Matrix, boolMatrix);
    // additional safety step: in principle 9999 could be the
    // minimal value, but it should not be considered
    for (int i = 0; i < ndmu; i++) {
        for (int j = 0; j < ngenmu; j++) {
            if (chi2Matrix(i, j) == 9999) {
                boolMatrix(i, j) = 0;
            }
        }
    }
    if (ndmu == 0) {
        for (int j = 0; j < ngenmu; j++) {
            if (genmu_kindOfMatching[j] != static_cast<Int_t>(GenMatchResults::GEN_OUTSIDE_CMS)) {
                genmu_kindOfMatching[j] = static_cast<Int_t>(GenMatchResults::ZERO_RECOS_FOUND);
            }
        }
    }

    // -------------------------------------
    // Assign residuals and gen match info
    // -------------------------------------
    for (Int_t j = 0; j < ngenmu; j++) {
        for (int i = 0; i < ndmu; i++) {
            if (boolMatrix(i, j) == 1.0) {
                dmu_hasGenMatch[i] = true;
                dmu_genMatchedIndex[i] = j;
                GenMatchResults matchResult = matchResults[{i, j}];
                dmu_propagationSurface[i] = static_cast<Int_t>(matchResult);
                genmu_kindOfMatching[j] = static_cast<Int_t>(matchResult);

                dmu_chi2[i] = chi2Matrix(i, j);
                dmu_chi2_x[i] = chi2_vectors[{i, j}][0];
                dmu_chi2_y[i] = chi2_vectors[{i, j}][1];
                dmu_chi2_z[i] = chi2_vectors[{i, j}][2];
                dmu_chi2_px[i] = chi2_vectors[{i, j}][3];
                dmu_chi2_py[i] = chi2_vectors[{i, j}][4];
                dmu_chi2_pz[i] = chi2_vectors[{i, j}][5];

                dmu_reco_final_x_err[i] = error_vectors[{i, j}][0];
                dmu_reco_final_y_err[i] = error_vectors[{i, j}][1];
                dmu_reco_final_z_err[i] = error_vectors[{i, j}][2];
                dmu_reco_final_px_err[i] = error_vectors[{i, j}][3];
                dmu_reco_final_py_err[i] = error_vectors[{i, j}][4];
                dmu_reco_final_pz_err[i] = error_vectors[{i, j}][5];

                GlobalTrajectoryParameters genFinalParams = propagatedTrajectories[{i, j}].first;
                GlobalTrajectoryParameters recoFinalParams = propagatedTrajectories[{i, j}].second;
                GlobalPoint recoFinalVertex = recoFinalParams.position();
                GlobalVector recoFinalMomentum = recoFinalParams.momentum();
                GlobalPoint genFinalVertex = genFinalParams.position();
                GlobalVector genFinalMomentum = genFinalParams.momentum();

                dmu_reco_final_r[i] = recoFinalVertex.perp();
                dmu_reco_final_theta[i] = recoFinalVertex.theta();
                dmu_reco_final_phi[i] = recoFinalVertex.phi();
                dmu_reco_final_p_r[i] = recoFinalMomentum.perp();
                dmu_reco_final_p_theta[i] = recoFinalMomentum.theta();
                dmu_reco_final_p_phi[i] = recoFinalMomentum.phi();

                genmu_final_r[j] = genFinalVertex.perp();
                genmu_final_theta[j] = genFinalVertex.theta();
                genmu_final_phi[j] = genFinalVertex.phi();
                genmu_final_p_r[j] = genFinalMomentum.perp();
                genmu_final_p_theta[j] = genFinalMomentum.theta();
                genmu_final_p_phi[j] = genFinalMomentum.phi();
                break;
            }
        }
        // fill genmu_kindOfMatching for the failed matches
        if (genmu_kindOfMatching[j] == -1 && ndmu > 0 &&
            (static_cast<Int_t>(matchResults[{0, j}]) < 0)) {
            Int_t first = static_cast<Int_t>(matchResults[{0, j}]);
            if (first == -1) {
                break;
            }
            bool all_same = true;
            for (Int_t i = 1; i < ndmu; i++) {
                if (static_cast<Int_t>(matchResults[{i, j}]) != first) {
                    all_same = false;
                    break;
                }
            }
            if (all_same) {
                genmu_kindOfMatching[j] = first;
            }
        }
    }
    //-> Fill tree
    tree_out->Fill();
}
DEFINE_FWK_MODULE(ntuplizer_test);