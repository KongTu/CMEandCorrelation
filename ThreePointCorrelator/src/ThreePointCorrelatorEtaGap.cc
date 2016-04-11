// -*- C++ -*-
//
// Package:    ThreePointCorrelatorEtaGap
// Class:      ThreePointCorrelatorEtaGap
// 
/**\class ThreePointCorrelatorEtaGap ThreePointCorrelatorEtaGap.cc CMEandCorrelation/ThreePointCorrelatorEtaGap/src/ThreePointCorrelatorEtaGap.cc

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
//
// class declaration
//

class ThreePointCorrelatorEtaGap : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorEtaGap(const edm::ParameterSet&);
      ~ThreePointCorrelatorEtaGap();

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
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* aveQ3[2][2];//calculate the correction on v2

//end v2

      TH1D* QvsdEta[48][3];
      TH1D* XY_real[48][3];TH1D* XY_imag[48][3];
      TH1D* XZ_real[48][3];TH1D* XZ_imag[48][3];
      TH1D* YZ_real[48][3];TH1D* YZ_imag[48][3];
      TH1D* X_real[48][3]; TH1D* X_imag[48][3];
      TH1D* Y_real[48][3]; TH1D* Y_imag[48][3];
      TH1D* Z_real[48][3]; TH1D* Z_imag[48][3];

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double holesize_;

      bool useCentrality_;
      bool reverseBeam_;
      bool messAcceptance_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;

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
ThreePointCorrelatorEtaGap::ThreePointCorrelatorEtaGap(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  
  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");
  
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  reverseBeam_ = iConfig.getUntrackedParameter<bool>("reverseBeam");
  messAcceptance_ = iConfig.getUntrackedParameter<bool>("messAcceptance");

  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);

  holesize_ = iConfig.getUntrackedParameter<double>("holesize");

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");

}


ThreePointCorrelatorEtaGap::~ThreePointCorrelatorEtaGap()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
double get3RealDup(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;

}
double get3ImagDup(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double get2RealDup(double R1, double R2, double I1, double I2){

  double t1 = R1*R2;
  double t2 = I1*I2;

  return t1-t2;
}
double get2ImagDup(double R1, double R2, double I1, double I2){

  double t1 = R1*I2;
  double t2 = R2*I1;

  return t1+t2;
}    
// ------------ method called for each event  ------------
void
ThreePointCorrelatorEtaGap::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  if(bestvz < vzLow_ || bestvz > vzHigh_ ) return;
  
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

  const int NetaBins = etaBins_.size() - 1 ;
  const int NdEtaBins = dEtaBins_.size() - 1;


// initialize Qcos and Qsin

  double QcosTRK = 0.;
  double QsinTRK = 0.;
  int QcountsTrk = 0;

  double Q1[NetaBins][2][2];
  int Q1_count[NetaBins][2];

  for(int i = 0; i < NetaBins; i++){
    for(int j = 0; j < 2; j++){
      Q1_count[i][j] = 0;
      for(int k = 0; k < 2; k++){
        Q1[i][j][k] = 0.0;
      }
    }
  }

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
        if(fabs(trk.eta()) > 2.4 || trk.pt() < 0.4) continue;
        if( messAcceptance_ ) {
          if( ( trk.phi() < (0.0 + holesize_) && trk.phi() > (0.0 - holesize_) )    ||
              ( trk.phi() < (2.09 + holesize_) && trk.phi() > (2.09 - holesize_) )  ||
              ( trk.phi() < (-2.09 + holesize_) && trk.phi() > (-2.09 - holesize_) )     ) continue;
        }
        nTracks++;  

        trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off

        QcosTRK += cos( 2*trk.phi() );
        QsinTRK += sin( 2*trk.phi() );
        QcountsTrk++;

        for(int eta = 0; eta < NetaBins; eta++){
          if( trk.eta() > etaBins_[eta] && trk.eta() < etaBins_[eta+1] ){

            if( trk.charge() == 1){
              Q1[eta][0][0] += cos( trk.phi() );
              Q1[eta][0][1] += sin( trk.phi() );
              Q1_count[eta][0]++;
            }
            else if( trk.charge() == -1){
              Q1[eta][1][0] += cos( trk.phi() );
              Q1[eta][1][1] += sin( trk.phi() );
              Q1_count[eta][1]++;
            }
          }
        }     
  } 

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

//loop over calo towers (HF)

  double Q3[2][2];
  double ETT[2];

  for(int i = 0; i < 2; i++){
    ETT[i] = 0.;
    for(int j = 0; j < 2; j++){
      Q3[i][j] = 0.;
    }
  }

  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
        if( messAcceptance_ ){ 
          if( caloPhi < -1.5 ) continue;
        }

        hfPhi->Fill( caloPhi );//make sure if messAcceptance is on or off

        if( reverseBeam_ ) caloEta = -hit.eta();
        if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
          
            Q3[0][0] += w*cos( -2*caloPhi );
            Q3[0][1] += w*sin( -2*caloPhi );
            ETT[0] += w;
        }
        else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

            Q3[1][0] += w*cos( -2*caloPhi );
            Q3[1][1] += w*sin( -2*caloPhi );
            ETT[1] += w;

        }
        else{continue;}
  }

  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){
    
      if( ieta == jeta ) continue;

      double deltaEta = fabs(etaBins_[jeta] - etaBins_[ieta]);

      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBins_[deta] && deltaEta < dEtaBins_[deta+1] ){
          for(int sign = 0; sign < 2; sign++ ){
            if( Q1_count[ieta][sign] == 0 || Q1_count[jeta][sign] == 0 || ETT[0] == 0.0 ) continue; //USE HF plus first

            double Q_real = get3RealDup(Q1[ieta][sign][0]/Q1_count[ieta][sign],Q1[jeta][sign][0]/Q1_count[jeta][sign],Q3[0][0]/ETT[0], Q1[ieta][sign][1]/Q1_count[ieta][sign], Q1[jeta][sign][1]/Q1_count[jeta][sign], Q3[0][1]/ETT[0]);
            QvsdEta[deta][sign]->Fill( Q_real );  

            double XY_real_temp = get2RealDup(Q1[ieta][sign][0], Q1[jeta][sign][0], Q1[ieta][sign][1], Q1[jeta][sign][1]);
            double XY_imag_temp = get2ImagDup(Q1[ieta][sign][0], Q1[jeta][sign][0], Q1[ieta][sign][1], Q1[jeta][sign][1]);
            
            double XZ_real_temp = get2RealDup(Q1[ieta][sign][0], Q3[0][0], Q1[ieta][sign][1], Q3[0][1]);
            double XZ_imag_temp = get2ImagDup(Q1[ieta][sign][0], Q3[0][0], Q1[ieta][sign][1], Q3[0][1]);

            double YZ_real_temp = get2RealDup(Q1[jeta][sign][0], Q3[0][0], Q1[jeta][sign][1], Q3[0][1]);
            double YZ_imag_temp = get2ImagDup(Q1[jeta][sign][0], Q3[0][0], Q1[jeta][sign][1], Q3[0][1]);

            XY_real[deta][sign]->Fill( XY_real_temp/(Q1_count[ieta][sign]*Q1_count[jeta][sign]), Q1_count[ieta][sign]*Q1_count[jeta][sign] );
            XY_imag[deta][sign]->Fill( XY_imag_temp/(Q1_count[ieta][sign]*Q1_count[jeta][sign]), Q1_count[ieta][sign]*Q1_count[jeta][sign] );
            
            XZ_real[deta][sign]->Fill( XZ_real_temp/(Q1_count[ieta][sign]*ETT[0]), Q1_count[ieta][sign]*ETT[0] );
            XZ_imag[deta][sign]->Fill( XZ_imag_temp/(Q1_count[ieta][sign]*ETT[0]), Q1_count[ieta][sign]*ETT[0] );

            YZ_real[deta][sign]->Fill( YZ_real_temp/(Q1_count[jeta][sign]*ETT[0]), Q1_count[jeta][sign]*ETT[0] );
            YZ_imag[deta][sign]->Fill( YZ_imag_temp/(Q1_count[jeta][sign]*ETT[0]), Q1_count[jeta][sign]*ETT[0] );

            double X_real_temp = Q1[ieta][sign][0]; double X_imag_temp = Q1[ieta][sign][1]; 
            double Y_real_temp = Q1[jeta][sign][0]; double Y_imag_temp = Q1[jeta][sign][1]; 
            double Z_real_temp = Q3[0][0];          double Z_imag_temp = Q3[0][1]; 

            X_real[deta][sign]->Fill( X_real_temp/Q1_count[ieta][sign], Q1_count[ieta][sign]);    
            Y_real[deta][sign]->Fill( Y_real_temp/Q1_count[jeta][sign], Q1_count[jeta][sign]);    
            Z_real[deta][sign]->Fill( Z_real_temp/ETT[0], ETT[0]);  
          
            X_imag[deta][sign]->Fill( X_imag_temp/Q1_count[ieta][sign], Q1_count[ieta][sign]);    
            Y_imag[deta][sign]->Fill( Y_imag_temp/Q1_count[jeta][sign], Q1_count[jeta][sign]);    
            Z_imag[deta][sign]->Fill( Z_imag_temp/ETT[0], ETT[0]);

          }
        }
      }
    }
  }

/*
calculate v2 using 3 sub-events method:
 */

  aveQ3[0][0]->Fill( Q3[0][0]/ETT[0] );
  aveQ3[0][1]->Fill( Q3[0][1]/ETT[0] );
  
  aveQ3[1][0]->Fill( Q3[1][0]/ETT[1] );
  aveQ3[1][1]->Fill( Q3[1][1]/ETT[1] );

  QcosTRK = QcosTRK/QcountsTrk;
  QsinTRK = QsinTRK/QcountsTrk;

  double QaQc = get2RealDup(Q3[1][0]/ETT[1], QcosTRK/QcountsTrk, Q3[1][1]/ETT[1], QsinTRK/QcountsTrk );
  double QaQb = get2RealDup(Q3[1][0]/ETT[1], Q3[0][0]/ETT[0], Q3[1][1]/ETT[1], -Q3[0][1]/ETT[0]);//an extra minus sign 
  double QcQb = get2RealDup(QcosTRK/QcountsTrk, Q3[0][0]/ETT[0], QsinTRK/QcountsTrk, Q3[0][1]/ETT[0]);

  c2_ac->Fill( QaQc );
  c2_cb->Fill( QcQb  );
  c2_ab->Fill( QaQb );

}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorEtaGap::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);

  const int NdEtaBins = dEtaBins_.size() - 1;
//HF:
  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 20000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 20000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 20000,-1,1);

  for(int i = 0; i < 2; i++ ){
      for(int j = 0; j < 2; j++){

        aveQ3[i][j] = fs->make<TH1D>(Form("aveQ3_%d_%d",i, j), ";aveQ3", 20000, -1.0, 1.0);
      }
  }

  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 2; sign++){
      
      QvsdEta[deta][sign] = fs->make<TH1D>(Form("QvsdEta_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      
      XY_real[deta][sign] = fs->make<TH1D>(Form("XY_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      XZ_real[deta][sign] = fs->make<TH1D>(Form("XZ_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      YZ_real[deta][sign] = fs->make<TH1D>(Form("YZ_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      
      XY_imag[deta][sign] = fs->make<TH1D>(Form("XY_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      XZ_imag[deta][sign] = fs->make<TH1D>(Form("XZ_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      YZ_imag[deta][sign] = fs->make<TH1D>(Form("YZ_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      
      X_real[deta][sign] = fs->make<TH1D>(Form("X_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      Y_real[deta][sign] = fs->make<TH1D>(Form("Y_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      Z_real[deta][sign] = fs->make<TH1D>(Form("Z_real_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      
      X_imag[deta][sign] = fs->make<TH1D>(Form("X_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      Y_imag[deta][sign] = fs->make<TH1D>(Form("Y_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);
      Z_imag[deta][sign] = fs->make<TH1D>(Form("Z_imag_%d_%d",deta,sign), "", 20000,-1.0,1.0);

    }
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorEtaGap::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorEtaGap::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorEtaGap::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorEtaGap::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorEtaGap::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorEtaGap::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorEtaGap);
