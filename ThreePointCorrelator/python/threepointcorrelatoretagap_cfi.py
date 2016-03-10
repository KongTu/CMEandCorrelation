import FWCore.ParameterSet.Config as cms

				
ana = cms.EDAnalyzer('ThreePointCorrelatorEtaGap',
                                                  vertexSrc = cms.string('offlinePrimaryVertices'),
                                                  trackSrc = cms.InputTag('generalTracks'),
                                                  towerSrc = cms.InputTag('towerMaker'),
                                                  offlineDCA = cms.untracked.double(3.0),
                                                  offlineptErr = cms.untracked.double(0.1),
                                                  useCentrality = cms.untracked.bool(False),
                                                  reverseBeam = cms.untracked.bool(False),
                                                  Nmin = cms.untracked.int32(185),
                                                  Nmax = cms.untracked.int32(220),
                                                  vzLow = cms.untracked.double(-15.0),
                                                  vzHigh = cms.untracked.double(+15.0),
                                                  etaLowHF = cms.untracked.double(4.4),
                                                  etaHighHF = cms.untracked.double(5.0),
                                                  etaBins = cms.untracked.vdouble(-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,
                                                                                  -1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
                                                                                  -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,
                                                                                  1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4),
                                                  dEtaBins = cms.untracked.vdouble(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                                                                                   1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,
                                                                                   2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8)

)
