import FWCore.ParameterSet.Config as cms

				
ana_tree = cms.EDAnalyzer('ThreePointCorrelatorTree',
                                                  vertexSrc = cms.string('offlinePrimaryVertices'),
                                                  trackSrc = cms.InputTag('generalTracks'),
                                                  towerSrc = cms.InputTag('towerMaker'),
                                                  useCentrality = cms.untracked.bool(False)

)
