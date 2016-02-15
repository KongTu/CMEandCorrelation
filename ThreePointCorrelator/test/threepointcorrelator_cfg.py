import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*'
                          #'HLT_PAPixelTracks_Multiplicity130_v*',
                          #'HLT_PAPixelTracks_Multiplicity160_v*',
                          #'HLT_PAPixelTracks_Multiplicity190_v*'
                          #'HLT_PAPixelTracks_Multiplicity220_v*'
]

process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
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
        'root://xrootd3.cmsaf.mit.edu//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Rereco_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_YQB.root'
    )
)

process.ana = cms.EDAnalyzer('ThreePointCorrelator',
						  vertexSrc = cms.string('offlinePrimaryVertices'),
                          trackSrc = cms.InputTag('generalTracks')
)


process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))
process.p = cms.Path(   process.hltHM*
                        process.ana
                                    )