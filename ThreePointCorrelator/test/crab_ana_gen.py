### this is an example for running on RECO
### the name must be changed crab.cfg for actual running

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
outputName = 'CME_QvsdEta_pPb_EPOS_GEN_v1'
config.General.requestName = outputName
config.General.workArea = outputName
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'threepointcorrelatorgen_cfg.py'
config.Data.inputDBS = 'phys03'
#EPOS:
config.Data.inputDataset = '/ReggeGribovPartonMCfix_EposLHC_5TeV_pPb/davidlw-RecoSkim_ReTracking_v4_5M-5cde49c8740ff28f897f533d05a99dbc/USER'
#HIJING:
#config.Data.inputDataset = '/Hijing_PPb502_MinimumBias/davidlw-RecoSkim_ReTracking_v4_10M-5cde49c8740ff28f897f533d05a99dbc/USER'

config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Site.storageSite = 'T2_US_MIT'
