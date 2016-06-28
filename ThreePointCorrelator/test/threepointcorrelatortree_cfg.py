import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.Timing = cms.Service("Timing")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'STARTHI53_V17::All'
#globaltag for pPb data
#process.GlobalTag.globaltag = 'GR_P_V43F::All'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_ReTracking_v18/28b2b9cce04ec3f20baeb96fbd2295a8/pPb_HM_1000_1_YyR.root'
#'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_HM_1000_1_BPd.root'
'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_MB_100_2_s6H.root'
)
)
process.load("HLTrigger.HLTanalyzers.HLTBitAnalyser_cfi")
process.hltanalysis = cms.EDAnalyzer("HLTBitAnalyzer",
    HLTProcessName = cms.string('HLT'),
                                     RunParameters = cms.PSet(
        HistogramFile = cms.untracked.string('hltbitanalysis.root'),
 	Monte         = cms.bool(True)
    
    ),
    UseTFileService = cms.untracked.bool(True),
    mctruth                         = cms.InputTag("genParticles::SIM"),
    genEventInfo                    = cms.InputTag("generator::SIM"),
    OfflinePrimaryVertices0     = cms.InputTag('offlinePrimaryVertices'),
    simhits                         = cms.InputTag("g4SimHits"),
    hltresults = cms.InputTag("TriggerResults","","HLT"),
    l1GctHFBitCounts = cms.InputTag("gctDigis"),
    l1GctHFRingSums = cms.InputTag("gctDigis"),
    l1GtObjectMapRecord = cms.InputTag("L1GtObjectMap","","RECO"),
    l1GtReadoutRecord = cms.InputTag("gtDigis","","RECO"),
    l1extramc = cms.string('l1extraParticles'),
    l1extramu = cms.string('l1extraParticles')
)

process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatortree_cfi")
#define the cuts
process.TFileService = cms.Service("TFileService",fileName = cms.string("test_tree.root"))
process.p = cms.Path( process.ana_tree )
#process.hltbitanalysis = cms.EndPath(process.hltanalysis)
