import FWCore.ParameterSet.Config as cms
import os

#################################  CONSTANTS  #################################################
ERA = 'C'
ReReco = True
###############################################################################################


process = cms.Process("demo")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

# Debug printout and summary.
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  # Set up multi-threaded run. Must be consistent with config.JobType.numCores in crab_cfg.py.
  #numberOfThreads=cms.untracked.uint32(8)
)

from Configuration.AlCa.GlobalTag import GlobalTag

# Select number of events to be processed
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
listOfFiles = ['']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
if ERA in 'ABCD' and not ReReco: gTag = '124X_dataRun3_PromptAnalysis_v1'
if ERA in 'ABCDE' and ReReco:    gTag = '124X_dataRun3_v15'
if ERA in 'FG' and not ReReco:   gTag = '124X_dataRun3_PromptAnalysis_v2'
process.GlobalTag = GlobalTag(process.GlobalTag, gTag)

## Define the process to run 
## 
process.load("DisplacedMuons-FrameWork.Ntuplizer.Cosmics_ntuples_MiniAOD_cfi")

process.p = cms.EndPath(process.ntuples)
