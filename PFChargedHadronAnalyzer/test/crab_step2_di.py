from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'step2_dipigun01'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'step2_dipigun.py'
config.JobType.maxMemoryMB = 2500
#config.JobType.numCores = 4

#config.Data.inputDataset = '/SinglePiPt100Eta1p6_2p8/PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'step2_dipigun01'
config.Data.ignoreLocality = True
#config.Data.totalUnits = 100000
config.Data.userInputFiles = open('step1_dipigun01_0p5_0p2_list').readlines()+open('step1_dipigun01_0p2_0p1_list').readlines()+open('step1_dipigun01_0p1_0p05_list').readlines()+open('step1_dipigun01_0p05_0p02_list').readlines()+open('step1_dipigun01_0p02_0p01_list').readlines()+open('step1_dipigun01_0p01_0p00_list').readlines()
#Data.outputPrimaryDataset

config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.ignoreGlobalBlacklist = True
config.Site.whitelist = ['T2_US_*']
#config.Site.whitelist = ['T3_US_FNALLPC']
#config.Site.blacklist = ['T2_US_Florida']
