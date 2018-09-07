from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

# config.General.requestName = 'step1_dipigun01_0p5_0p2'
# config.General.requestName = 'step1_dipigun01_0p2_0p1'
#config.General.requestName = 'step1_dipigun01_0p1_0p05'
# config.General.requestName = 'step1_dipigun01_0p05_0p02'
# config.General.requestName = 'step1_dipigun01_0p02_0p01'
config.General.requestName = 'step1_dipigun01_0p01_0p00'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'step1_dipigun.py'
config.JobType.maxMemoryMB = 2000
config.JobType.numCores = 4

#config.Data.inputDataset = '/SinglePiPt100Eta1p6_2p8/PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 750
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
# config.Data.outputDatasetTag  = 'step1_dipigun01_0p5_0p2'
# config.Data.outputDatasetTag  = 'step1_dipigun01_0p2_0p1'
#config.Data.outputDatasetTag  = 'step1_dipigun01_0p1_0p05'
# config.Data.outputDatasetTag  = 'step1_dipigun01_0p05_0p02'
# config.Data.outputDatasetTag  = 'step1_dipigun01_0p02_0p01'
config.Data.outputDatasetTag  = 'step1_dipigun01_0p01_0p00'
config.Data.ignoreLocality = True
config.Data.totalUnits = 150000

config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.ignoreGlobalBlacklist = True
config.Site.whitelist = ['T2_US_*']
#config.Site.whitelist = ['T3_US_FNALLPC']
#config.Site.blacklist = ['T2_US_Florida']
