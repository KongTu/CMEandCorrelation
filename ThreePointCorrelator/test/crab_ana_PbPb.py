### this is an example for running on RECO
### the name must be changed crab.cfg for actual running

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.requestName = 'CME_QvsdEta_PbPb_50_100_v2'
config.General.workArea = 'CME_QvsdEta_PbPb_50_100_v2'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'threepointcorrelatoretagap_cfg.py'

config.Data.inputDBS = 'phys03'
#MB
config.Data.inputDataset = '/HIMinBiasUPC/davidlw-PR2011_MBPPRereco_TRKANASKIM_v6-ccf03100d177f42de0f9cdc7627799d3/USER'

config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Site.storageSite = 'T2_US_MIT'
