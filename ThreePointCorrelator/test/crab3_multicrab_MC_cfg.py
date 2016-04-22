from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import FWCore.ParameterSet.Config as cms
#load the cfi file and rewrite cross section parameter each time:
process = cms.Process('Demo')
process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatoretagap_cfi")

outputName = "multicrab_CME_QvsdEta_pPb_EPOS_v36"

config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'threepointcorrelatoretagap_cfg.py'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
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
   
   sampleName = ['/EPOS_GEN-SIM_1M_batch12/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch13/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch14/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch15/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch16/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch17/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch18/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch19/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch2/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch20/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch21/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch22/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch23/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch24/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch25/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch3/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch4/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch5/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch6/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch7/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch8/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch9/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/EPOS_GEN-SIM_1M_batch10/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
	         '/EPOS_GEN-SIM_1M_batch11/davidlw-RecoSkim_ReTracking_v4-5cde49c8740ff28f897f533d05a99dbc/USER',
		 '/ReggeGribovPartonMCfix_EposLHC_5TeV_pPb/davidlw-RecoSkim_ReTracking_v4_5M-5cde49c8740ff28f897f533d05a99dbc/USER'
		]
   
   for num in range(0,26):
		
       	print 'double check this is the %r sample of EPOS with multiplicity 120 to 150' %(num)
      		
       	RequestName = outputName + "_" + str(num)
       	DataSetName = sampleName[num]
       	config.General.requestName = RequestName
       	config.Data.inputDataset = DataSetName
       	submit(config)

# python crab3_ppTrackingAnalyzer.py to execute 
# ./multicrab -c status -w crab_projects/ to check status 
