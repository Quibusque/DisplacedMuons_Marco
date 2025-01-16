import CRABClient
from CRABClient.UserUtilities import config, getLumiListInValidFiles
#from FWCore.DataStructs.LumiList import LumiList
from FWCore.PythonUtilities.LumiList import LumiList

config = config()

# General
config.General.workArea = 'crab_projects'
config.General.requestName = 'CosmicsPPreco-CosmicDTLocalReco_AOD-Ntuples'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'CosmicsData_PPreco-CosmicDTLocalReco_runNtuplizer_cfg.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['ntuples.root']

# Data
config.Data.inputDataset = '/Cosmics/fernance-CosmicPPreco-CosmicDTLocalReco-1ff72641f835c7d0ae35987c50731b57/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/rlopezru/Cosmics2024/'
config.Data.publication = False
config.Data.outputDatasetTag = 'CosmicsPPreco-CosmicDTLocalReco_Ntuples'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
