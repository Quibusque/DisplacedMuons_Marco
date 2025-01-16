import FWCore.ParameterSet.Config as cms
import os

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
nEvents = 462173
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
listOfFiles = ['/store/data/Run2023C/NoBPTX/MINIAOD/22Sep2023_v2-v1/50000/00f7e24b-cddd-484c-93d7-6434d6726085.root',
               '/store/data/Run2023C/NoBPTX/MINIAOD/22Sep2023_v2-v1/50000/0107ed00-aa45-410b-9603-9dac4fd79bb2.root']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
gTag = '130X_mcRun3_2023_realistic_v14'
process.GlobalTag = GlobalTag(process.GlobalTag, gTag)

## Define the process to run 
## 
process.load("DisplacedMuons-FrameWork.Ntuplizer.CosmicsData_ntuples_MiniAOD_cfi")

process.p = cms.EndPath(process.ntuples)
