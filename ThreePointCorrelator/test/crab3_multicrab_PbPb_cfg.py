from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatoretagap_cfi")

ntrkRange = [0,35,60,90,120,150,185,220,260]

outputName = "multicrab_CME_QvsdEta_PbPb_50_100_v1"

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'threepointcorrelatoretagap_cfg.py'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 15
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = outputName

config.Site.storageSite = 'T2_US_MIT'

if __name__ == '__main__':
   from CRABAPI.RawCommand import crabCommand
   from CRABClient.ClientExceptions import ClientException
   from httplib import HTTPException

   config.General.workArea = outputName

   def submit(config):
      try:
           crabCommand('submit', config = config)
      except HTTPException as hte:
           print "Failed submitting task: %s" % (hte.headers)
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)
   
   for num in range(0,8):
		
	print 'double check that multiplicity range is fram %r to %r' % (ntrkRange[paths],ntrkRange[paths+1])
      		
	process.ana.Nmin = ntrkRange[paths]
	process.ana.Nmax = ntrkRange[paths+1]
       	RequestName = outputName + "_" + str(num)
       	DataSetName = '/HIMinBiasUPC/davidlw-PR2011_MBPPRereco_TRKANASKIM_v6-ccf03100d177f42de0f9cdc7627799d3/USER'
       	config.General.requestName = RequestName
       	config.Data.inputDataset = DataSetName
       	submit(config)

# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
