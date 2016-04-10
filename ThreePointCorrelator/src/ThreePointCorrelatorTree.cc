// -*- C++ -*-
//
// Package:    ThreePointCorrelatorTree
// Class:      ThreePointCorrelatorTree
//
/**\class ThreePointCorrelatorTree ThreePointCorrelatorTree.cc MitHig/ThreePointCorrelatorTree/src/ThreePointCorrelatorTree.cc

   Description: <one line class summary>

   Implementation:
   Prepare the Track Tree for analysis
*/
//
// Original Author:  Yilmaz Yetkin, Yen-Jie Lee
// Updated: Frank Ma, Matt Nguyen
//         Created:  Tue Sep 30 15:14:28 CEST 2008
// $Id: ThreePointCorrelatorTree.cc,v 1.55 2013/06/11 20:58:09 yjlee Exp $

// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
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
#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
  int nTower;
  int cbin;

  double vtxZ;
  double trkEta[MAXTRACKS];
  double trkPhi[MAXTRACKS];
  double trkPt[MAXTRACKS];
  double trkPtError[MAXTRACKS];
  double trkDCAxy[MAXTRACKS];
  double trkDCAz[MAXTRACKS];

  double towEta[MAXTRACKS];
  double towPhi[MAXTRACKS];
  double towEt[MAXTRACKS];
  double towEnergy[MAXTRACKS];

};

class ThreePointCorrelatorTree : public edm::EDAnalyzer {

public:
  explicit ThreePointCorrelatorTree(const edm::ParameterSet&);
  ~ThreePointCorrelatorTree();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  // ----------member data ---------------------------

  bool useCentrality_;

  edm::InputTag trackSrc_;
  edm::InputTag towerSrc_;
  string vertexSrc_;

  // Root object
  TTree* trackTree_;

  TrackEvent pev_;


};

//--------------------------------------------------------------------------------------------------
ThreePointCorrelatorTree::ThreePointCorrelatorTree(const edm::ParameterSet& iConfig)

{

  useCentrality_       = iConfig.getUntrackedParameter<bool>  ("doCentrality",false);
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  vertexSrc_ = iConfig.getParameter<string>("vertexSrc");

}

//--------------------------------------------------------------------------------------------------
ThreePointCorrelatorTree::~ThreePointCorrelatorTree()
{
}

//--------------------------------------------------------------------------------------------------
void
ThreePointCorrelatorTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  pev_.nEv = (int)iEvent.id().event();
  pev_.nRun = (int)iEvent.id().run();
  pev_.nLumi = (int)iEvent.luminosityBlock();
  pev_.nBX = (int)iEvent.bunchCrossing();

  fillTracks(iEvent, iSetup);
  trackTree_->Fill();

}

//--------------------------------------------------------------------------------------------------
void
ThreePointCorrelatorTree::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

  pev_.vtxZ = (double) bestvz;

  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

  if( useCentrality_ ){

    CentralityProvider * centProvider = 0;
    if (!centProvider) centProvider = new CentralityProvider(iSetup);
    centProvider->newEvent(iEvent,iSetup);
    int hiBin = centProvider->getBin();
    pev_.cbin = (int) hiBin;

  }

  pev_.nTrk = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

      const reco::Track & trk = (*tracks)[it];

      math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

      double dzvtx = trk.dz(bestvtx);
      double dxyvtx = trk.dxy(bestvtx);
      double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
      double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
      if(!trk.quality(reco::TrackBase::highPurity)) continue;

      pev_.trkEta[pev_.nTrk] = trk.eta();
      pev_.trkPhi[pev_.nTrk] = trk.phi();
      pev_.trkPt[pev_.nTrk] = trk.pt();
      pev_.trkDCAz[pev_.nTrk] = fabs(dzvtx/dzerror);
      pev_.trkDCAxy[pev_.nTrk] = fabs(dxyvtx/dxyerror);
      pev_.trkPtError[pev_.nTrk] = fabs(trk.ptError())/trk.pt();

      pev_.nTrk++;//this is not Ntrkoffline, it needs to give a cut

  } 

   pev_.nTower = 0;
   for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
        double energy = hit.emEnergy() + hit.hadEnergy();

        pev_.towEta[pev_.nTower] = caloEta;
        pev_.towPhi[pev_.nTower] = caloPhi;
        pev_.towEt[pev_.nTower] = w;
        pev_.towEnergy[pev_.nTower] = energy;

        pev_.nTower++;
  }

}
// ------------ method called once each job just before starting event loop  ------------
void
ThreePointCorrelatorTree::beginJob()
{  
  edm::Service<TFileService> fs;

  trackTree_ = fs->make<TTree>("trackTree","v1");

  // event
  trackTree_->Branch("nEv",&pev_.nEv,"nEv/I");
  trackTree_->Branch("nLumi",&pev_.nLumi,"nLumi/I");
  trackTree_->Branch("nBX",&pev_.nBX,"nBX/I");
  trackTree_->Branch("nRun",&pev_.nRun,"nRun/I");

  if( useCentrality_ )  trackTree_->Branch("cbin", &pev_.cbin,"cbin/I");

  trackTree_->Branch("nTower", &pev_.nTower,"nTower/I");
  trackTree_->Branch("nTrk", &pev_.nTrk,"nTrk/I");
  trackTree_->Branch("vtxZ", &pev_.vtxZ,"vtxZ/D");
  trackTree_->Branch("trkEta", &pev_.trkEta,"trkEta[nTrk]/D");
  trackTree_->Branch("trkPhi", &pev_.trkPhi,"trkPhi[nTrk]/D");
  trackTree_->Branch("trkPt", &pev_.trkPt,"trkPt[nTrk]/D");
  trackTree_->Branch("trkPtError", &pev_.trkPtError,"trkPtError[nTrk]/D");
  trackTree_->Branch("trkDCAz", &pev_.trkDCAz,"trkDCAz[nTrk]/D");
  trackTree_->Branch("trkDCAxy", &pev_.trkDCAxy,"trkDCAxy[nTrk]/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void
ThreePointCorrelatorTree::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorTree);