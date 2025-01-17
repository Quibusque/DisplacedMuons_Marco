import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('ntuples.root'),
    isAOD                         = cms.bool(False),
    isCosmics                     = cms.bool(True),
    isMCSignal                    = cms.bool(False),
    EventInfo                     = cms.InputTag("generator"),
    RunInfo                       = cms.InputTag("generator"),
    BeamSpot                      = cms.InputTag("offlineBeamSpot"),
    displacedGlobalCollection     = cms.InputTag("displacedGlobalMuons"),
    displacedStandAloneCollection = cms.InputTag("displacedStandAloneMuons"),
    displacedMuonCollection       = cms.InputTag("slimmedDisplacedMuons"),

    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)
