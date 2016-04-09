// -*- C++ -*-
//
// Package:    TrackAnalyzer
// Class:      TrackAnalyzer
//
/**\class TrackAnalyzer TrackAnalyzer.cc MitHig/TrackAnalyzer/src/TrackAnalyzer.cc

   Description: <one line class summary>

   Implementation:
   Prepare the Track Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie Lee
// Updated: Frank Ma, Matt Nguyen
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: TrackAnalyzer.cc,v 1.55 2013/06/11 20:58:09 yjlee Exp $

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <functional>

// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
// #include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
// #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
// #include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
// #include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
// #include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class decleration
//

#define PI 3.14159265358979

#define MAXTRACKS 60000
#define MAXVTX 100
#define MAXQUAL 5
#define MAXMATCH 5

struct TrackEvent{

  // event information
  int nRun;
  int nEv;
  int nLumi;
  int nBX;
  int N; // multiplicity variable

  int nTrk;
  int nVtx;
  int nParticle;

};

class TrackAnalyzer : public edm::EDAnalyzer {

public:
  explicit TrackAnalyzer(const edm::ParameterSet&);
  ~TrackAnalyzer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // ----------member data ---------------------------

  bool doTrack_;
  bool doTrackExtra_;
  bool doSimTrack_;
  bool doSimVertex_;
  bool fillSimTrack_;
  bool doPFMatching_;
  bool useQuality_;
  bool doDeDx_;
  bool doDebug_;
  bool doMVA_;
  // bool associateChi2_;
  bool doHighestPtVertex_;
  bool doTrackVtxWImpPar_;

  double trackPtMin_;
  double trackVtxMaxDistance_;
  std::vector<std::string> qualityStrings_;
  std::string qualityString_;

  double simTrackPtMin_;
  bool fiducialCut_; 
  edm::InputTag trackSrc_;
  std::string mvaSrc_;
  edm::InputTag particleSrc_;
  edm::InputTag tpFakeSrc_;
  edm::InputTag tpEffSrc_;
  edm::InputTag pfCandSrc_;
  edm::InputTag DeDxSrc_;
  edm::InputTag associatorMap_;

  vector<string> vertexSrc_;
  edm::InputTag simVertexSrc_;

  const TrackerGeometry* geo_;
  edm::Service<TFileService> fs;
  edm::ESHandle < ParticleDataTable > pdt;
  edm::Handle<TrackingParticleCollection> trackingParticles;

  edm::InputTag beamSpotProducer_;

  // Root object
  TTree* trackTree_;

  TrackEvent pev_;


};

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)

{

  doTrack_             = iConfig.getUntrackedParameter<bool>  ("doTrack",true);

  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");

  particleSrc_ = iConfig.getParameter<edm::InputTag>("particleSrc");

  vertexSrc_ = iConfig.getParameter<vector<string> >("vertexSrc");


}

//--------------------------------------------------------------------------------------------------
TrackAnalyzer::~TrackAnalyzer()
{
}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get tracker geometry
  //  cout <<"StartFill"<<endl;

  edm::ESHandle<TrackerGeometry> tGeo;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeo);
  geo_ = tGeo.product();
  iSetup.getData(pdt);

  //  cout <<"Got data"<<endl;
  pev_.nEv = (int)iEvent.id().event();
  pev_.nRun = (int)iEvent.id().run();
  pev_.nLumi = (int)iEvent.luminosityBlock();
  pev_.nBX = (int)iEvent.bunchCrossing();
  pev_.N = 0;

  // pev_.nv = 0;
  pev_.nParticle = 0;
  pev_.nTrk = 0;

  //cout <<"Fill Vtx"<<endl;

  //cout <<"Fill Tracks"<<endl;
  if (doTrack_) fillTracks(iEvent, iSetup);
  //cout <<"Tracks filled!"<<endl;
  //cout <<"SimTracks filled!"<<endl;
  trackTree_->Fill();
  //cout <<"Tree filled!"<<endl;

}

//--------------------------------------------------------------------------------------------------
void
TrackAnalyzer::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<vector<Track> > etracks;
  iEvent.getByLabel(trackSrc_, etracks);

}
// ------------ method called once each job just before starting event loop  ------------
void
TrackAnalyzer::beginJob()
{

  trackTree_ = fs->make<TTree>("trackTree","v1");

  // event
  trackTree_->Branch("nEv",&pev_.nEv,"nEv/I");
  trackTree_->Branch("nLumi",&pev_.nLumi,"nLumi/I");
  trackTree_->Branch("nBX",&pev_.nBX,"nBX/I");
  trackTree_->Branch("nRun",&pev_.nRun,"nRun/I");
  trackTree_->Branch("N",&pev_.N,"N/I");


}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);