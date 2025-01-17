import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('MC_Ntuples.root'),
    isAOD                         = cms.bool(False),
    isCosmics                     = cms.bool(False),
    isMCSignal                    = cms.bool(True),
    EventInfo                     = cms.InputTag("generator"),
    RunInfo                       = cms.InputTag("generator"),
    BeamSpot                      = cms.InputTag("offlineBeamSpot"),
    displacedGlobalCollection     = cms.InputTag("displacedGlobalMuons"),
    displacedStandAloneCollection = cms.InputTag("displacedStandAloneMuons"),
    displacedMuonCollection       = cms.InputTag("slimmedDisplacedMuons"),
    prunedGenParticles            = cms.InputTag("prunedGenParticles"),

    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)
