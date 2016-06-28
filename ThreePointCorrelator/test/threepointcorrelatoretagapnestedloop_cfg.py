import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.Timing = cms.Service("Timing")
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

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

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_BPd.root'    
#'root://xrootd-cms.infn.it//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_YyR.root'
#'root://cmsxrootd.fnal.gov//store/user/davidlw/ReggeGribovPartonMCfix_EposLHC_5TeV_pPb/RecoSkim_ReTracking_v4_5M/5cde49c8740ff28f897f533d05a99dbc/pPb_MB_100_2_s6H.root'
#'root://xrootd-cms.infn.it//store/user/davidlw/Hijing_PPb502_MinimumBias/RecoSkim_ReTracking_v4_10M/5cde49c8740ff28f897f533d05a99dbc/pPb_MB_100_1_yjL.root'
'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_HM_1000_1_BPd.root'
#'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_MB_100_2_s6H.root'
#'root://cmsxrootd.fnal.gov//store/user/davidlw/PAMinBiasUPC/PA2013_FlowCorr_PromptReco_MB_Gplus_Rereco_ReTracking_v18/25c9a89be536a41c8ccb3c75e9fd9358/pPb_HM_1000_1_Bgt.root'
)
)

process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatoretagap_cfi")

#define the cuts
process.ana_nested.Nmin = 185
process.ana_nested.Nmax = 220
process.ana_nested.doEffCorrection = True

process.TFileService = cms.Service("TFileService",fileName = cms.string("test_nested.root"))
process.p = cms.Path(   process.hltHM*
			process.ana_nested
                                    )
