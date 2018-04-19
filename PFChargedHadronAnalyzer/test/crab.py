from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'rereco_HGCALfix06'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'rereco.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDataset = '/SinglePiPt100Eta1p6_2p8/PhaseIITDRFall17DR-noPUFEVT_93X_upgrade2023_realistic_v2-v1/GEN-SIM-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'rereco_HGCALfix06'
config.Data.ignoreLocality = True

config.Site.storageSite = 'T3_US_Brown'
config.Site.ignoreGlobalBlacklist = True
#config.Site.whitelist = ['T2_US_*']
config.Site.blacklist = ['T2_US_Florida']
