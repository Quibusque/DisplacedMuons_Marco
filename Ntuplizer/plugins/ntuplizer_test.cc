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
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
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

typedef std::pair<TrajectoryStateOnSurface, double> TsosPath;

bool hasMotherWithPdgId_test(const reco::Candidate* particle, int pdgId) {
    // Loop on mothers, if any, and return true if a mother with the specified pdgId is found
    for (size_t i = 0; i < particle->numberOfMothers(); i++) {
        const reco::Candidate* mother = particle->mother(i);
        if (mother->pdgId() == pdgId || hasMotherWithPdgId_test(mother, pdgId)) {
            return true;
        }
    }
    return false;
}

bool propagateToSurface_test(float radius, float minZ, float maxZ, const FreeTrajectoryState& fts,
                             const Propagator* propagator, TsosPath& tsosPath) {
    const Surface::RotationType dummyRot;
    Cylinder::CylinderPointer theTargetCylinder =
        Cylinder::build(Surface::PositionType(0., 0., 0.), dummyRot, radius);
    Plane::PlanePointer theTargetPlaneMin =
        Plane::build(Surface::PositionType(0., 0., minZ), dummyRot);
    Plane::PlanePointer theTargetPlaneMax =
        Plane::build(Surface::PositionType(0., 0., maxZ), dummyRot);

    bool goodResult = true;
    tsosPath = propagator->propagateWithPath(fts, *theTargetCylinder);
    if (!tsosPath.first.isValid()) {
        goodResult = false;
    } else if (tsosPath.first.globalPosition().z() < theTargetPlaneMin->position().z()) {
        tsosPath = propagator->propagateWithPath(fts, *theTargetPlaneMin);
        if (!tsosPath.first.isValid() ||
            tsosPath.first.globalPosition().perp() > theTargetCylinder->radius()) {
            goodResult = false;
        }
    } else if (tsosPath.first.globalPosition().z() > theTargetPlaneMax->position().z()) {
        tsosPath = propagator->propagateWithPath(fts, *theTargetPlaneMax);
        if (!tsosPath.first.isValid() ||
            tsosPath.first.globalPosition().perp() > theTargetCylinder->radius()) {
            goodResult = false;
        }
    }
    return goodResult;
}

std::pair<GlobalTrajectoryParameters, bool> myGenMatching_test(const reco::GenParticle& genParticle,
                                                               const reco::Track* recoTrack,
                                                               const MagneticField* magField,
                                                               const Propagator* propagator) {
    // The Gen Matching is done by propagating the particles to a common surface and doing a check
    // based on the chi-square of the difference between the momenta

    // First barrel muon chambers for the radius and first CSC station for the z
    Float_t radius = 420.0;
    Float_t minZ = -700.0;
    Float_t maxZ = 700.0;

    // Initial states
    GlobalPoint genVertex(genParticle.vx(), genParticle.vy(), genParticle.vz());
    GlobalVector genMomentum(genParticle.px(), genParticle.py(), genParticle.pz());
    int genCharge = genParticle.charge();
    FreeTrajectoryState genFTS(genVertex, genMomentum, genCharge, magField);
    TsosPath genTsosPath;

    GlobalPoint recoVertex(recoTrack->vx(), recoTrack->vy(), recoTrack->vz());
    GlobalVector recoMomentum(recoTrack->px(), recoTrack->py(), recoTrack->pz());
    int recoCharge = recoTrack->charge();
    FreeTrajectoryState recoFTS(recoVertex, recoMomentum, recoCharge, magField);
    TsosPath recoTsosPath;

    // Propagate the particles to the target surface
    bool genPropagationGood =
        propagateToSurface_test(radius, minZ, maxZ, genFTS, propagator, genTsosPath);
    bool recoPropagationGood =
        propagateToSurface_test(radius, minZ, maxZ, recoFTS, propagator, recoTsosPath);

    if (!genPropagationGood || !recoPropagationGood)
        return std::make_pair(GlobalTrajectoryParameters(), false);

    // Calculate the residuals
    GlobalPoint genPositionPropagated = genTsosPath.first.globalPosition();
    GlobalVector genMomentumPropagated = genTsosPath.first.globalMomentum();
    GlobalPoint recoPositionPropagated = recoTsosPath.first.globalPosition();
    GlobalVector recoMomentumPropagated = recoTsosPath.first.globalMomentum();

    GlobalPoint residualPosition =
        GlobalPoint(genPositionPropagated.x() - recoPositionPropagated.x(),
                    genPositionPropagated.y() - recoPositionPropagated.y(),
                    genPositionPropagated.z() - recoPositionPropagated.z());

    GlobalVector residualMomentum = genMomentumPropagated - recoMomentumPropagated;

    GlobalTrajectoryParameters residualParams(residualPosition, residualMomentum, genCharge,
                                              magField);

    return std::make_pair(residualParams, true);
}

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

    std::ofstream chi2file;

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
    edm::ESGetToken<Propagator, TrackingComponentsRecord> thePropToken;

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
        thePropToken =
            esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("propagator")));
        theMagFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
    }

    triggerBits_ = consumes<edm::TriggerResults>(parameters.getParameter<edm::InputTag>("bits"));
}

// Destructor
ntuplizer_test::~ntuplizer_test() {}

// beginJob (Before first event)
void ntuplizer_test::beginJob() {
    std::cout << "Begin Job" << std::endl;

    // Open the chi2 file
    chi2file.open("chi2.txt");
    chi2file << "deltaR_residual, residual, match?\n";

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

    // Gen Matching branches
    if (isMCSignal) {
        tree_out->Branch("dmu_hasGenMatch", dmu_hasGenMatch, "dmu_hasGenMatch[ndmu]/O");
        tree_out->Branch("dmu_genMatchedIndex", dmu_genMatchedIndex, "dmu_genMatchedIndex[ndmu]/I");
    }
    // Trigger branches
    for (unsigned int ihlt = 0; ihlt < HLTPaths_.size(); ihlt++) {
        tree_out->Branch(TString(HLTPaths_[ihlt]), &triggerPass[ihlt]);
    }
}

// endJob (After event loop has finished)
void ntuplizer_test::endJob() {
    std::cout << "End Job" << std::endl;
    chi2file.close();
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
    const Propagator* propagator = nullptr;
    if (isMCSignal) {
        iEvent.getByToken(prunedGenToken, prunedGen);
        propagator = &iSetup.getData(thePropToken);
        if (!dynamic_cast<const SteppingHelixPropagator*>(propagator)) {
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
    Int_t nprugenmu = 0;
    int goodGenMuons_indices[10] = {-1};
    int n_goodGenMuons = 0;
    if (isMCSignal) {
        for (unsigned int i = 0; i < prunedGen->size(); i++) {
            const reco::GenParticle& p(prunedGen->at(i));
            if (abs(p.pdgId()) == 13 && p.status() == 1) {  // Check if the particle is a muon
                // Check if the muon has a mother with pdgId 1023
                if (hasMotherWithPdgId_test(&p, 1023)) {
                    goodGenMuons_indices[n_goodGenMuons] = i;
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
    std::vector<Float_t> deltaR_residuals;
    std::vector<Float_t> residuals;
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
            // if (isAOD) {
            //     // Number of DT+CSC segments
            //     unsigned int nsegments = 0;
            //     for (trackingRecHit_iterator hit = outerTrack->recHitsBegin();
            //          hit != outerTrack->recHitsEnd(); ++hit) {
            //         if (!(*hit)->isValid()) continue;
            //         DetId id = (*hit)->geographicalId();
            //         if (id.det() != DetId::Muon) continue;
            //         if (id.subdetId() == MuonSubdetId::DT || id.subdetId() == MuonSubdetId::CSC)
            //         {
            //             nsegments++;
            //         }
            //     }
            //     dmu_dsa_nsegments[ndmu] = nsegments;
            // }
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
            dmu_hasGenMatch[ndmu] = false;
            dmu_genMatchedIndex[ndmu] = -1;

            const reco::Track* candidateTrack = nullptr;
            if (dmuon.isGlobalMuon()) {
                continue;
            } else if (dmuon.isStandAloneMuon()) {
                if (dmu_dsa_dxy[ndmu] > 0.3) {
                    continue;
                }
                candidateTrack = (dmuon.standAloneMuon()).get();
            } else if (dmuon.isTrackerMuon()) {
                continue;
            }

            TVector3 candidate_vertex = TVector3();
            candidate_vertex.SetXYZ(candidateTrack->vx(), candidateTrack->vy(),
                                    candidateTrack->vz());
            dmu_vertex_r[ndmu] = sqrt(candidate_vertex.X() * candidate_vertex.X() +
                                      candidate_vertex.Y() * candidate_vertex.Y());
            dmu_vertex_z[ndmu] = candidate_vertex.Z();

            for (Int_t j = 0; j < n_goodGenMuons; j++) {
                const reco::GenParticle& p(prunedGen->at(goodGenMuons_indices[j]));
                TVector3 gen_vertex = TVector3();
                gen_vertex.SetXYZ(p.vx(), p.vy(), p.vz());
                std::pair<GlobalTrajectoryParameters, bool> genMatchResult =
                    myGenMatching_test(p, candidateTrack, magField, propagator);
                GlobalPoint residualPoint = genMatchResult.first.position();
                GlobalVector residualMomentum = genMatchResult.first.momentum();
                bool matchGood = genMatchResult.second;
                // Residual is the magnitude of the point summed to momentum
                Float_t residual = matchGood ? residualPoint.mag() + residualMomentum.mag() : -1;
                Float_t deltaR_residual = reco::deltaR(*candidateTrack, p);
                deltaR_residuals.push_back(deltaR_residual);
                residuals.push_back(residual);
            }
        }
        ndmu++;

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

        //-> Fill tree
        tree_out->Fill();
    }

    /// Check the residuals vectors
    // Print each pair of values to the file
    for (size_t i = 0; i < deltaR_residuals.size(); ++i) {
        chi2file << deltaR_residuals[i] << ", " << residuals[i] << '\n';
    }

    // Check that the index of the minimum residual is the same as the index of the minimum
    // deltaR_residual
    auto min_deltaR_it = std::min_element(deltaR_residuals.begin(), deltaR_residuals.end());
    auto min_residual_it = std::min_element(residuals.begin(), residuals.end());

    size_t min_deltaR_index = std::distance(deltaR_residuals.begin(), min_deltaR_it);
    size_t min_residual_index = std::distance(residuals.begin(), min_residual_it);

    chi2file << "Min index match: " << (min_deltaR_index == min_residual_index) << '\n';

    // // Debugging: just check if this is reasonable
    // if (matchedMuons.size() == 2 &&
    //     matchedMuons[0]->charge() != matchedMuons[1]->charge()) {
    //     TLorentzVector muon1, muon2;
    //     muon1.SetPtEtaPhiM(matchedMuons[0]->pt(), matchedMuons[0]->eta(),
    //                         matchedMuons[0]->phi(), 0.105);
    //     muon2.SetPtEtaPhiM(matchedMuons[1]->pt(), matchedMuons[1]->eta(),
    //                         matchedMuons[1]->phi(), 0.105);
    //     double invariantMass = (muon1 + muon2).M();
    //     std::cout << "Invariant mass of the two gen-matched muons: " << invariantMass
    //                 << std::endl;
    // }
}
DEFINE_FWK_MODULE(ntuplizer_test);
