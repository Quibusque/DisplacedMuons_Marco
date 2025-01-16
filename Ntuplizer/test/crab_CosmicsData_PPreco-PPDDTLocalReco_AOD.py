import CRABClient
from CRABClient.UserUtilities import config, getLumiListInValidFiles
#from FWCore.DataStructs.LumiList import LumiList
from FWCore.PythonUtilities.LumiList import LumiList

config = config()

# General
config.General.workArea = 'crab_projects'
config.General.requestName = 'CosmicsPPreco-PPDDTLocalReco_AOD-Ntuples'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CosmicsData_PPreco-PPDDTLocalReco_runNtuplizer_cfg.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['ntuples.root']

# Data
config.Data.inputDataset = '/Cosmics/fernance-CosmicPPreco-PPDTLocalReco-7307c75268af57e7de289394389e476f/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/rlopezru/Cosmics2024/'
config.Data.publication = False
config.Data.outputDatasetTag = 'CosmicsPPreco-PPDDTLocalReco_Ntuples'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
