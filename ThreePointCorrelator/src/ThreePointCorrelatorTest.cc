// -*- C++ -*-
//
// Package:    ThreePointCorrelatorTest
// Class:      ThreePointCorrelatorTest
// 
/**\class ThreePointCorrelatorTest ThreePointCorrelatorTest.cc CMEandCorrelation/ThreePointCorrelatorTest/src/ThreePointCorrelatorTest.cc

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
#include <stdlib.h>


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

class ThreePointCorrelatorTest : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorTest(const edm::ParameterSet&);
      ~ThreePointCorrelatorTest();

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

      TH1D* averageCosHF[2];//calculate the correction on v2
      TH1D* averageSinHF[2];
//end v2

      //culmulants: 
      TH1D* QvsdEta[2];


      TH1D* XY[2];
      TH1D* XZ[2];
      TH1D* YZ[2];

      TH1D* XYcount;
      TH1D* XZcount;
      TH1D* YZcount;

      TH1D* X[2];
      TH1D* Y[2];
      TH1D* Z[2];

      TH1D* Xcount;
      TH1D* Ycount;
      TH1D* Zcount;

      //double particle product, corresponds to the first 3 terms, <QaQb><Qc>, <QaQc><Qb>, and <QbQc>
      //like/unlike sign, real/imaginary
      TH2D* QaQbvsdEta[3][2];
      TH2D* NaNbvsdEta[3];//count is the same for real and imaginary

      //type, like/unlike sign, real/imaginary
      TH2D* QaQcvsdEta[3][3][2];
      TH2D* QbvsdEta[3][3][2];
      TH2D* NaNcvsdEta[3][3];

      //like/unlike sign, real/imaginary corresponds to <Qa><Qb>
      TH2D* QaSinglevsdEta[3][2];
      TH2D* QbSinglevsdEta[3][2];

      //single particle sum, corresponds to the last term in correction, <Qc>
      TH1D* HFcosSum[2];
      TH1D* HFsinSum[2];
      TH1D* weightSum[2];

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
ThreePointCorrelatorTest::ThreePointCorrelatorTest(const edm::ParameterSet& iConfig)

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


ThreePointCorrelatorTest::~ThreePointCorrelatorTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double get3Real(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;

}
double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double get2Real(double R1, double R2, double I1, double I2){

  double t1 = R1*R2;
  double t2 = I1*I2;

  return t1-t2;
}
double get2Imag(double R1, double R2, double I1, double I2){

  double t1 = R1*I2;
  double t2 = R2*I1;

  return t1+t2;
}
    
// ------------ method called for each event  ------------
void
ThreePointCorrelatorTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  //const int binSize_ = etaBins_.size() - 1 ;

// initialize Qcos and Qsin
  double Qcos[2];
  double Qsin[2];
  int Qcounts[2];
  for(int eta = 0; eta < 2; eta++){
      Qcounts[eta] = 0;
      Qcos[eta] = 0.;
      Qsin[eta] = 0.;
  }



  double x[1000];
  double y[1000];
  double z[1000];
  double HFqVcos = 0.;
  double HFqVsin = 0.;
  int HFcounts = 0;
  for(int i = 0; i < 100; i++){

    x[i]= fRand(-1.5,3.14);
    Qcos[0] += cos(x[i]);
    Qsin[0] += sin(x[i]);
    Qcounts[0]++;
    trkPhi->Fill(x[i]);
  }

  for(int i = 0; i < 200; i++){

    y[i]= fRand(-1.5,3.14);
    Qcos[1] += cos(y[i]);
    Qsin[1] += sin(y[i]);
    Qcounts[1]++;


  }
  for(int i = 0; i < 50; i++){

    z[i]= fRand(-1.5,3.14);
    HFqVcos += cos( -2*z[i] );
    HFqVsin += sin( -2*z[i] );
    HFcounts++;

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
          if( trk.phi() < -1.5 ) continue;
        }
        nTracks++; 

        // trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off

        // if( trk.eta() > -2.4 && trk.eta() < -2.0 ){

        //         Qcos[0] += cos( trk.phi() );
        //         Qsin[0] += sin( trk.phi() );
        //         Qcounts[0]++;

        // }
        // else if( trk.eta() > 2.0 && trk.eta() < 2.4 ){
             
        //      if(trk.charge() == 1){

        //         Qcos[1] += cos( trk.phi() );
        //         Qsin[1] += sin( trk.phi() );
        //         Qcounts[1]++;
        //      }
        // }
        // else if( trk.eta() > -1.0 && trk.eta() < 1.0 ){

        //      if(trk.charge() == 1 ){

        //       HFqVcos += cos( -2*trk.phi() );
        //       HFqVsin += sin( -2*trk.phi() );
        //       HFcounts++;
        //      }

        // }
        // else {continue;}
  }



  //   for(unsigned i = 0; i < towers->size(); ++i){

  //       const CaloTower & hit= (*towers)[i];

  //       double caloEta = hit.eta();
  //       double caloPhi = hit.phi();
  //       if( messAcceptance_ ){ 
  //         if( caloPhi < -1.5 ) continue;
  //       }

  //       caloPhi = fRand(-1.5,3.14);

  //       hfPhi->Fill( caloPhi );//make sure if messAcceptance is on or off

  //       if( fabs(caloEta) > etaLowHF_ && fabs(caloEta) < etaHighHF_ ){

  //         HFqVcos += cos( -2*caloPhi );
  //         HFqVsin += sin( -2*caloPhi );          
  //         HFcounts++;
  //       }
  // }

  if( nTracks <= Nmin_ || nTracks > Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

  double XY_real = get2Real(Qcos[0], Qcos[1], Qsin[0], Qsin[1]);
  double XY_imag = get2Imag(Qcos[0], Qcos[1], Qsin[0], Qsin[1]); 

  double XZ_real = get2Real(Qcos[0], HFqVcos, Qsin[0], HFqVsin);
  double XZ_imag = get2Imag(Qcos[0], HFqVcos, Qsin[0], HFqVsin);

  double YZ_real = get2Real(Qcos[1], HFqVcos, Qsin[1], HFqVsin);
  double YZ_imag = get2Imag(Qcos[1], HFqVcos, Qsin[1], HFqVsin);

  int XY_count = Qcounts[0]*Qcounts[1];
  int XZ_count = Qcounts[0]*HFcounts;
  int YZ_count = Qcounts[1]*HFcounts;

  double X_real = Qcos[0];
  double Y_real = Qcos[1];
  double Z_real = HFqVcos;

  double X_imag = Qsin[0];
  double Y_imag = Qsin[1];
  double Z_imag = HFqVsin;

  int X_count = Qcounts[0];
  int Y_count = Qcounts[1];
  int Z_count = HFcounts;

  XY[0]->Fill( XY_real ); XY[1]->Fill( XY_imag );
  XZ[0]->Fill( XZ_real ); XZ[1]->Fill( XZ_imag );
  YZ[0]->Fill( YZ_real ); YZ[1]->Fill( YZ_imag );

  XYcount->Fill( XY_count ); XZcount->Fill( XZ_count ); YZcount->Fill( YZ_count ); 
   
  X[0]->Fill( X_real ); X[1]->Fill( X_imag );
  Y[0]->Fill( Y_real ); Y[1]->Fill( Y_imag );
  Z[0]->Fill( Z_real ); Z[1]->Fill( Z_imag );

  Xcount->Fill( X_count ); Ycount->Fill( Y_count ); Zcount->Fill( Z_count );

  //self normalize the Qvectors from HF:
  HFqVcos = HFqVcos/HFcounts;
  HFqVsin = HFqVsin/HFcounts;

//3 point correlator
//self normalize the Qvectors;
  for(int eta = 0; eta < 2; eta++){

    if( Qcounts[eta] == 0 ) continue;

    Qcos[eta] = Qcos[eta]/Qcounts[eta];
    Qsin[eta] = Qsin[eta]/Qcounts[eta];

  }

  double real_totalQplusplus = get3Real(Qcos[0],Qcos[1], HFqVcos, Qsin[0], Qsin[1], HFqVsin ); 
  double imag_totalQplusplus = get3Imag(Qcos[0],Qcos[1], HFqVcos, Qsin[0], Qsin[1], HFqVsin );

  QvsdEta[0]->Fill( real_totalQplusplus );
  QvsdEta[1]->Fill( imag_totalQplusplus );

 
}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorTest::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);

  for(int real = 0; real < 2; real++ ){
      QvsdEta[real] = fs->make<TH1D>(Form("QvsdEta_%d",real),";Q_{#phi_{1}}Q_{#phi_{2}}Q^{*}_{2#phi_{3}}", 4000,-2.0-0.000005,2.0-0.000050 );
  }
  
  for(int real = 0; real < 2; real++){

    XY[real] = fs->make<TH1D>(Form("XY_%d", real), ";XY", 2000000, -100000, 100000);
    XZ[real] = fs->make<TH1D>(Form("XZ_%d", real), ";XZ", 2000000, -100000, 100000);
    YZ[real] = fs->make<TH1D>(Form("YZ_%d", real), ";YZ", 2000000, -100000, 100000);

    X[real] = fs->make<TH1D>(Form("X_%d", real), ";X", 2000000, -100, 100);
    Y[real] = fs->make<TH1D>(Form("Y_%d", real), ";Y", 2000000, -100, 100);
    Z[real] = fs->make<TH1D>(Form("Z_%d", real), ";Z", 2000000, -100, 100);

  }

  XYcount = fs->make<TH1D>("XYcount", ";XYcount", 100000, 0,100000);
  XZcount = fs->make<TH1D>("XZcount", ";XZcount", 100000, 0,100000);
  YZcount = fs->make<TH1D>("YZcount", ";YZcount", 100000, 0,100000);

  Xcount = fs->make<TH1D>("Xcount", ";Xcount", 2000, 0,2000);
  Ycount = fs->make<TH1D>("Ycount", ";Ycount", 2000, 0,2000);
  Zcount = fs->make<TH1D>("Zcount", ";Zcount", 2000, 0,2000);


}
// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorTest::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorTest);
