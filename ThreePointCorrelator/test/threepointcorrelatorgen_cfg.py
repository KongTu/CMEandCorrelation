import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )

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
#'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_HM_1000_1_BPd.root'
'file:/afs/cern.ch/work/z/ztu/CME/CMSSW_5_3_20/src/CMEandCorrelation/ThreePointCorrelator/test/pPb_MB_100_2_s6H.root'
)
)

process.load("CMEandCorrelation.ThreePointCorrelator.threepointcorrelatorgen_cfi")

#define the cuts
process.ana.Nmin = 40
process.ana.Nmax = 10000

process.TFileService = cms.Service("TFileService",fileName = cms.string("test.root"))
process.p = cms.Path( process.ana )
                                    
