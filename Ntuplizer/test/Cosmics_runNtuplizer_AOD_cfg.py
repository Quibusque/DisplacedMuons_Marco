import FWCore.ParameterSet.Config as cms
import os

#################################  CONSTANTS  #################################################
ERA = 'F'
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
if ERA in 'ABCD': gTag = '124X_dataRun3_PromptAnalysis_v1'
if ERA in 'EFG':  gTag = '124X_dataRun3_Prompt_v10'
process.GlobalTag = GlobalTag(process.GlobalTag, gTag)

## Define the process to run 
## 
process.load("DisplacedMuons-FrameWork.Ntuplizer.Cosmics_ntuples_AOD_cfi")

process.p = cms.EndPath(process.ntuples)
