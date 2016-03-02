import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.Timing = cms.Service("Timing")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*',
                          'HLT_PAPixelTracks_Multiplicity130_v*',
                          'HLT_PAPixelTracks_Multiplicity160_v*'
                          #'HLT_PAPixelTracks_Multiplicity190_v*'
                          #'HLT_PAPixelTracks_Multiplicity220_v*'
]

process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_53_LV6::All'

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    #nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_BPd.root'    
#'root://xrootd-cms.infn.it//store/user/davidlw/HIMinBiasUPC/PR2011_MBPPRereco_TRKANASKIM_v6/ccf03100d177f42de0f9cdc7627799d3/pPb_HM_1000_1_lq0.root'
'root://xrootd-cms.infn.it//store/user/davidlw/HIMinBiasUPC/Skim_rereco_MB_pixeltracks_final_v2/9c1b4b9b6b9ff3e493a474ba7d01bc76/hiMB_1000_1_WLn.root'
)
)

process.ana = cms.EDAnalyzer('ThreePointCorrelatorEtaGap',
						  vertexSrc = cms.string('hiSelectedVertex'),
                          			  trackSrc = cms.InputTag('hiGeneralAndPixelTracks'),
						  towerSrc = cms.InputTag('towerMaker'),
					      	  #centralitySrc = cms.InputTag('hiCentrality'),

						  Nmin = cms.untracked.int32(1),
						  Nmax = cms.untracked.int32(10000)
)



process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))
process.p = cms.Path(   #process.hltHM*
                        #process.centralityBin*
			process.ana
                                    )
