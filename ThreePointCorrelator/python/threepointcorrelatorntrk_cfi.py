import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltHM.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*',
                  'HLT_PAPixelTracks_Multiplicity130_v*',
                  'HLT_PAPixelTracks_Multiplicity160_v*'
                  #'HLT_PAPixelTracks_Multiplicity190_v*'
                  #'HLT_PAPixelTracks_Multiplicity220_v*'
]

hltHM.andOr = cms.bool(True)
hltHM.throw = cms.bool(False)				


ana = cms.EDAnalyzer('ThreePointCorrelatorNtrk',
                                                  vertexSrc = cms.string('offlinePrimaryVertices'),
                                                  trackSrc = cms.InputTag('generalTracks'),
                                                  towerSrc = cms.InputTag('towerMaker'),
                                                  offlineDCA = cms.untracked.double(3.0),
                                                  offlineptErr = cms.untracked.double(0.1),
                                                  useCentrality = cms.untracked.bool(False),
                                                  reverseBeam = cms.untracked.bool(False),
                                                  messAcceptance = cms.untracked.bool(False),
					          useBothSide = cms.untracked.bool(False),
						  doEffCorrection = cms.untracked.bool(False),
						  holesize = cms.untracked.double(0.0),
						  Nmin = cms.untracked.int32(1),
                                                  Nmax = cms.untracked.int32(10000),
                                                  vzLow = cms.untracked.double(-15.0),
                                                  vzHigh = cms.untracked.double(+15.0),
                                                  etaLowHF = cms.untracked.double(4.4),
                                                  etaHighHF = cms.untracked.double(5.0),
					          ptLow = cms.untracked.double(0.3),
                                                  ptHigh = cms.untracked.double(3.0),		
						  holeLeft = cms.untracked.double(0.0),
                                                  holeRight = cms.untracked.double(0.6),
						  ntrkBins = cms.untracked.vint32(10,20,30,40,50,60,70,80,90,100,110,120,135,150,165,180,200,220,240,260,300)
)
