// -*- C++ -*-
//
// Package:    ThreePointCorrelatorNtrk
// Class:      ThreePointCorrelatorNtrk
// 
/**\class ThreePointCorrelatorNtrk ThreePointCorrelatorNtrk.cc CMEandCorrelation/ThreePointCorrelatorNtrk/src/ThreePointCorrelatorNtrk.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu,42 1-015,,
//         Created:  Mon Feb 15 10:38:19 CET 2016
// $Id$
//
//

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
//
// class declaration
//

class ThreePointCorrelatorNtrk : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorNtrk(const edm::ParameterSet&);
      ~ThreePointCorrelatorNtrk();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      TH1D* Ntrk;
      TH1D* evtWeightedQp3;
      TH1D* evtWeight;
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;
      
      TH1D* averageCosPlus;
      TH1D* averageSinPlus;
      TH1D* averageCosMinus;
      TH1D* averageSinMinus;
      
      //culmulants: with [0] denote plusHF, [1] minusHF, [2] bothHF;
      TH2D* QvsNtrkPlusPlus[3];
      TH2D* QvsNtrkMinusMinus[3];
      TH2D* QvsNtrkPlusMinus[3];
      TH2D* QvsNtrkMinusPlus[3];

      //two and single particle sum to correct the acceptance effect
      //HF:
      TH1D* HFcosSum[2];
      TH1D* HFsinSum[2];
      TH1D* weightSum[2]; 

      TH1D* TRKcosPlusSum;
      TH1D* TRKsinPlusSum;
      TH1D* TRK2cosPlusSum;
      TH1D* TRK2sinPlusSum;

      TH1D* TRKcosMinusSum;
      TH1D* TRKsinMinusSum;
      TH1D* TRK2cosMinusSum;
      TH1D* TRK2sinMinusSum;

      TH1D* cbinHist;

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;

      double offlineptErr_;
      double offlineDCA_;

      bool useCentrality_;
      bool useBothSide_;

      std::vector<double> ntrkBins_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ThreePointCorrelatorNtrk::ThreePointCorrelatorNtrk(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  
  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");
  
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  useBothSide_ = iConfig.getUntrackedParameter<bool>("useBothSide");

  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);

  ntrkBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ntrkBins");

}


ThreePointCorrelatorNtrk::~ThreePointCorrelatorNtrk()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//


//cos(p1-p2)
double get2RealduP(double cos1, double cos2, double sin1, double sin2){

  double t1 = cos1*cos2;
  double t2 = sin1*sin2;

  return t1+t2;
}

    
// ------------ method called for each event  ------------
void
ThreePointCorrelatorNtrk::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
  
  //first selection; vertices
  if(bestvz < -15.0 || bestvz > 15.0) return;
  
  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);
  
  if( useCentrality_ ){

    CentralityProvider * centProvider = 0;
    if (!centProvider) centProvider = new CentralityProvider(iSetup);
    centProvider->newEvent(iEvent,iSetup);
    int hiBin = centProvider->getBin();
    cbinHist->Fill( hiBin );
    if( hiBin < Nmin_ || hiBin >= Nmax_ ) return;

  }

// initialize Qcos and Qsin and Q2cos and Q2sin
// each array has 2 entries corresponding to one positive particle and one negative particle.

  double Qcos[2];
  double Qsin[2];
  double Q2cos[2];
  double Q2sin[2];
  int Qcounts[2];

  for(int sign = 0; sign < 2; sign++){

      Qcounts[sign] = 0;
      Qcos[sign] = 0.;
      Qsin[sign] = 0.;
      Q2cos[sign] = 0.;
      Q2sin[sign] = 0.;
    }

//this is for v2 of c particle, doesn't matter sign. 
  double QcosTRK = 0.;
  double QsinTRK = 0.;
  int QcountsTrk = 0;

  int nTracks = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];

     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if ( fabs(trk.eta()) > 2.4 || trk.pt() < 0.4  ) continue;
        nTracks++;   

        //for c particle v2;
        QcosTRK = QcosTRK + cos( 2*trk.phi() );
        QsinTRK = QsinTRK + sin( 2*trk.phi() );
        QcountsTrk++;

        if( trk.charge() == 1 ){

          Qcos[0] += cos( trk.phi() );
          Qsin[0] += sin( trk.phi() );

          Q2cos[0] += cos( 2*trk.phi() );
          Q2sin[0] += sin( 2*trk.phi() );

          Qcounts[0]++;

          TRKcosPlusSum->Fill( cos( trk.phi() ) );
          TRKsinPlusSum->Fill( sin( trk.phi() ) );

          TRK2cosPlusSum->Fill( cos( 2*trk.phi() ) );
          TRK2sinPlusSum->Fill( sin( 2*trk.phi() ) );

        }
        else if( trk.charge() == -1 ){

          Qcos[1] += cos( trk.phi() );
          Qsin[1] += sin( trk.phi() );

          Q2cos[1] += cos( 2*trk.phi() );
          Q2sin[1] += sin( 2*trk.phi() );

          Qcounts[1]++;

          TRKcosMinusSum->Fill( cos( trk.phi() ) );
          TRKsinMinusSum->Fill( sin( trk.phi() ) );

          TRK2cosMinusSum->Fill( cos( 2*trk.phi() ) );
          TRK2sinMinusSum->Fill( sin( 2*trk.phi() ) );

        }
        else{
          continue;
        }
  } 

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

//loop over calo towers (HF)

  double HFqVcos = 0.;
  double HFqVsin = 0.;
  int HFcounts = 0;
  double ETT = 0.;

  double HFqVcosPlus = 0.;
  double HFqVsinPlus = 0.;
  int HFminusCounts = 0;
  double ETTplus = 0.;

  double HFqVcosMinus = 0.;
  double HFqVsinMinus = 0.;
  int HFplusCounts = 0;
  double ETTminus = 0.;

  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );

        if( fabs(caloEta) > etaLowHF_ && fabs(caloEta) < etaHighHF_ ){

          HFqVcos = HFqVcos + w*cos( 2*caloPhi );
          HFqVsin = HFqVsin + w*sin( 2*caloPhi );
          
          HFcounts++;
          ETT += w;
        }
        if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

          HFcosSum[1]->Fill( w*cos( 2*caloPhi ) );
          HFsinSum[1]->Fill( w*sin( 2*caloPhi ) );
          weightSum[1]->Fill( w );

          HFqVcosMinus = HFqVcosMinus + w*cos( 2*caloPhi );
          HFqVsinMinus = HFqVsinMinus + w*sin( 2*caloPhi );
         
          HFminusCounts++; 
          ETTminus += w;

        }
        if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){

          HFcosSum[0]->Fill( w*cos( 2*caloPhi ) );
          HFsinSum[0]->Fill( w*sin( 2*caloPhi ) );
          weightSum[0]->Fill( w );
          
          HFqVcosPlus = HFqVcosPlus + w*cos( 2*caloPhi );
          HFqVsinPlus = HFqVsinPlus + w*sin( 2*caloPhi );
          
          HFplusCounts++;
          ETTplus += w;
        }
  }

//weight by ET() and renormalize by ETT in order to have dimensionless 

  if( HFplusCounts == 0 || HFminusCounts == 0 ) return;

//self normalize the Qvectors from HF:
  HFqVcosPlus = HFqVcosPlus/ETTplus;
  HFqVsinPlus = HFqVsinPlus/ETTplus;

  HFqVcosMinus = HFqVcosMinus/ETTminus;
  HFqVsinMinus = HFqVsinMinus/ETTminus;

  HFqVcos = HFqVcos/ETT;
  HFqVsin = HFqVsin/ETT;

  QcosTRK = QcosTRK/QcountsTrk;
  QsinTRK = QsinTRK/QcountsTrk;

  double QaQc = get2RealduP(HFqVcosMinus, QcosTRK, HFqVsinMinus, QsinTRK );
  double QaQb = get2RealduP(HFqVcosMinus, HFqVcosPlus, HFqVsinMinus, HFqVsinPlus);
  double QcQb = get2RealduP(QcosTRK, HFqVcosPlus, QsinTRK, HFqVsinPlus);

  c2_ac->Fill( QaQc );
  c2_cb->Fill( QcQb  );
  c2_ab->Fill( QaQb );

  double W2 = ETTminus*ETTplus;
  evtWeight->Fill( W2 );
  evtWeightedQp3->Fill( W2*QaQb );

//3 point correlator

//like-sign correlator:
  for(int type = 0; type < 3; type++){

    double tempHFcos = 0.0;
    double tempHFsin = 0.0;

    if( type == 0 ) {
      tempHFcos = HFqVcosPlus;
      tempHFsin = HFqVsinPlus;
    }
    else if ( type == 1 ){
      tempHFcos = HFqVcosMinus;
      tempHFsin = HFqVsinMinus;
    }
    else if ( type == 3 ){
      tempHFcos = HFqVcos;
      tempHFsin = HFqVsin;
    }
    else{
      return;
    }

    for(int sign = 0; sign < 2; sign++){

      double realPart_like = Qcos[sign]*Qcos[sign] - Qsin[sign]*Qsin[sign];
      double imagPart_like = 2*Qcos[sign]*Qsin[sign];

      realPart_like = realPart_like - Q2cos[sign];
      imagPart_like = imagPart_like - Q2sin[sign];

      realPart_like = realPart_like/(Qcounts[sign]*(Qcounts[sign] - 1));
      imagPart_like = imagPart_like/(Qcounts[sign]*(Qcounts[sign] - 1));

      double fQ = get2RealduP(realPart_like, tempHFcos, imagPart_like, tempHFsin);

      if( sign == 0 ) QvsNtrkPlusPlus[type]->Fill( nTracks, fQ );
      else if (sign == 1) QvsNtrkMinusMinus[type]->Fill( nTracks, fQ );

    }

  //un-like sign correlator:
    double realPart_unlike = Qcos[0]*Qcos[1] - Qsin[0]*Qsin[1];
    double imagPart_unlike = Qcos[0]*Qsin[1] + Qcos[1]*Qsin[0];

    realPart_unlike = realPart_unlike/(Qcounts[0]*Qcounts[1]);
    imagPart_unlike = imagPart_unlike/(Qcounts[0]*Qcounts[1]);

    double fQ_unlike = get2RealduP(realPart_unlike, tempHFcos, imagPart_unlike, tempHFsin);

    QvsNtrkPlusMinus[type]->Fill( nTracks, fQ_unlike );  
      
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorNtrk::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);

  double ntrkBinsFill[100];
  const int nNtrkBins = ntrkBins_.size() - 1;
  for(unsigned num = 0; num < ntrkBins_.size(); num++ ){

    ntrkBinsFill[num] = ntrkBins_[num];
  }

//HF:

  evtWeight = fs->make<TH1D>("evtWeight",";evtWeight", 10000000,0,5000);
  evtWeightedQp3 = fs->make<TH1D>("evtWeightedQp3",";evtWeightedQp3", 1000000,-50,50);
  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 10000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 10000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 10000,-1,1);
  
  for(int sign = 0; sign < 2; sign++){

    HFsinSum[sign] = fs->make<TH1D>(Form("HFsinSum_%d", sign), ";HFsinSum", 2000, -1.0, 1.0 );
    HFcosSum[sign] = fs->make<TH1D>(Form("HFcosSum_%d", sign), ";HFcosSum", 2000, -1.0, 1.0 );
    weightSum[sign] = fs->make<TH1D>(Form("weightSum_%d", sign), ";weightSum", 3000, 0, 150 );

  }
//TRK:
    for(int type = 0; type < 3; type++ ){
     
      QvsNtrkPlusPlus[type] = fs->make<TH2D>(Form("QvsNtrkPlusPlus_%d", type), ";N^{offline}_{trk};Q_{#phi_{1,+}}Q_{#phi_{2,+}}Q^{*}_{2#phi_{3}}", nNtrkBins, ntrkBinsFill, 20000,-0.1,0.1 );
      QvsNtrkMinusMinus[type] = fs->make<TH2D>(Form("QvsNtrkMinusMinus_%d", type), ";N^{offline}_{trk};Q_{#phi_{1,-}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", nNtrkBins, ntrkBinsFill, 20000,-0.1,0.1 );
      QvsNtrkPlusMinus[type] = fs->make<TH2D>(Form("QvsNtrkPlusMinus_%d", type), ";N^{offline}_{trk};Q_{#phi_{1,+}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", nNtrkBins, ntrkBinsFill, 20000,-0.1,0.1 );
      QvsNtrkMinusPlus[type] = fs->make<TH2D>(Form("QvsNtrkMinusPlus_%d", type), ";N^{offline}_{trk};Q_{#phi_{1,-}}Q_{#phi_{2,+}}Q^{*}_{2#phi_{3}}", nNtrkBins, ntrkBinsFill, 20000,-0.1,0.1 );

    }

    TRKcosPlusSum = fs->make<TH1D>("TRKcosPlusSum", ";TRKcosPlusSum", 20000, -1.0, 1.0 );
    TRKsinPlusSum = fs->make<TH1D>("TRKsinPlusSum", ";TRKsinPlusSum", 20000, -1.0, 1.0 );
    TRKcosMinusSum = fs->make<TH1D>("TRKcosMinusSum", ";TRKcosMinusSum", 20000, -1.0, 1.0 );
    TRKsinMinusSum = fs->make<TH1D>("TRKsinMinusSum", ";TRKsinMinusSum", 20000, -1.0, 1.0 );

    TRK2cosPlusSum = fs->make<TH1D>("TRK2cosPlusSum", ";TRK2cosPlusSum", 20000, -1.0, 1.0 );
    TRK2sinPlusSum = fs->make<TH1D>("TRK2sinPlusSum", ";TRK2sinPlusSum", 20000, -1.0, 1.0 );
    TRK2cosMinusSum = fs->make<TH1D>("TRK2cosMinusSum", ";TRK2cosMinusSum", 20000, -1.0, 1.0 );
    TRK2sinMinusSum = fs->make<TH1D>("TRK2sinMinusSum", ";TRK2sinMinusSum", 20000, -1.0, 1.0 );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorNtrk::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorNtrk::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorNtrk::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorNtrk::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorNtrk::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorNtrk::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorNtrk);
