from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatoretagap_cfi")

ntrkRange = [0,35,60,90,120]

outputName = "multicrab_CME_QvsdEta_pPb_MB_v2"

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
   
   sampleName = ["/PAMinBiasUPC/davidlw-PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18-25c9a89be536a41c8ccb3c75e9fd9358/USER",
		 "/PAMinBiasUPC/davidlw-PA2013_FlowCorr_PromptReco_MB_Gplus_ReTracking_v18-25c9a89be536a41c8ccb3c75e9fd9358/USER",
		 "/PAMinBiasUPC/davidlw-PA2013_FlowCorr_PromptReco_MB_Gplus_Reverse_ReTracking_v18-25c9a89be536a41c8ccb3c75e9fd9358/USER"
		]
   
   beam = [False,False,True]

   for paths in range(0,4):
   	for num in range(0,3):
		
		print 'double check that multiplicity range is fram %r to %r' $ (ntrkRange[paths],ntrkRange[paths+1])
       		print 'double check that with %r sample the reverse beam option is %r ' % (num, beam[num])
      		
		process.ana.Nmin = ntrkRange[paths]
		process.ana.Nmax = ntrkRange[paths+1]
		process.ana.reverseBeam = beam[num]       
       		RequestName = outputName + "_" + str(paths) + "_"  + str(num)
       		DataSetName = sampleName[num]
       		config.General.requestName = RequestName
       		config.Data.inputDataset = DataSetName
       		submit(config)

# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
