import FWCore.ParameterSet.Config as cms
import os
import re
import argparse


#Argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-input_dir', type=str, help='Directory that contains all the root files to be processed')
parser.add_argument('-out_file', type=str, help='Output file name')
args = parser.parse_args()

main_dir = "/eos/home-m/mcrucian/datasets/displaced_muons/"
my_dir = os.path.join(main_dir, args.input_dir)


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
listOfFiles = [
    os.path.join(root, file)
    for root, dirs, files in os.walk(my_dir)
    for file in files
    if file.endswith(".root")
]
listOfFiles = ['file:'+f for f in files]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
gTag = '130X_mcRun3_2023_realistic_v14'
process.GlobalTag = GlobalTag(process.GlobalTag, gTag)

## Define the process to run 
## 
process.load("DisplacedMuons-FrameWork.Ntuplizer.MC_ntuples_MiniAOD_cfi")

process.ntuples.nameOfOutput = args.out_file

process.p = cms.EndPath(process.ntuples)
