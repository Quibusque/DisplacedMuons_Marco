import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer_test',
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
    propagatorAlong               = cms.string('SteppingHelixPropagatorAlong'),
    propagatorOpposite            = cms.string('SteppingHelixPropagatorOpposite'),
    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)

SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    useTuningForL2Speed = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    endcapShiftInZNeg = cms.double(0.0),
    SetVBFPointer = cms.bool(False),
    AssumeNoMaterial = cms.bool(False),
    endcapShiftInZPos = cms.double(0.0),
    useInTeslaFromMagField = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    useEndcapShiftsInZ = cms.bool(False),
    sendLogWarning = cms.bool(False),
    useMatVolumes = cms.bool(True),
    debug = cms.bool(False),
    #This sort of works but assumes a measurement at propagation origin  
    ApplyRadX0Correction = cms.bool(True),
    useMagVolumes = cms.bool(True),
    returnTangentPlane = cms.bool(True)
)

SteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ComponentName = cms.string('SteppingHelixPropagatorOpposite'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('oppositeToMomentum'),
    useTuningForL2Speed = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    endcapShiftInZNeg = cms.double(0.0),
    SetVBFPointer = cms.bool(False),
    AssumeNoMaterial = cms.bool(False),
    endcapShiftInZPos = cms.double(0.0),
    useInTeslaFromMagField = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    useEndcapShiftsInZ = cms.bool(False),
    sendLogWarning = cms.bool(False),
    useMatVolumes = cms.bool(True),
    debug = cms.bool(False),
    #This sort of works but assumes a measurement at propagation origin  
    ApplyRadX0Correction = cms.bool(True),
    useMagVolumes = cms.bool(True),
    returnTangentPlane = cms.bool(True)
)
