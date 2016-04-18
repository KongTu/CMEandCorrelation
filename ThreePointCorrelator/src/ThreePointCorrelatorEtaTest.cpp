// -*- C++ -*-
//
// Package:    ThreePointCorrelatorEtaTest
// Class:      ThreePointCorrelatorEtaTest
// 
/**\class ThreePointCorrelatorEtaTest ThreePointCorrelatorEtaTest.cc CMEandCorrelation/ThreePointCorrelatorEtaTest/src/ThreePointCorrelatorEtaTest.cc

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

#include "CMEandCorrelation/ThreePointCorrelator/interface/ThreePointCorrelatorBase.h"

ThreePointCorrelatorEtaTest::ThreePointCorrelatorEtaTest(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  
  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");
  
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  useBothSide_ = iConfig.getUntrackedParameter<bool>("useBothSide");  
  reverseBeam_ = iConfig.getUntrackedParameter<bool>("reverseBeam");
  messAcceptance_ = iConfig.getUntrackedParameter<bool>("messAcceptance");

  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");
  ptLow_ = iConfig.getUntrackedParameter<double>("ptLow");
  ptHigh_ = iConfig.getUntrackedParameter<double>("ptHigh");
  holeLeft_ = iConfig.getUntrackedParameter<double>("holeLeft");
  holeRight_ = iConfig.getUntrackedParameter<double>("holeRight");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);
  
  holesize_ = iConfig.getUntrackedParameter<double>("holesize");

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");

}


ThreePointCorrelatorEtaTest::~ThreePointCorrelatorEtaTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
ThreePointCorrelatorEtaTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

// initialize Qcos and Qsin

  double QcosTRK = 0.;
  double QsinTRK = 0.;
  double QcountsTrk = 0;

  double Q1[NetaBins][2][2];
  double Q1_count[NetaBins][2];

  double Q2[NetaBins][2][2];
  double Q2_count[NetaBins][2];

  for(int i = 0; i < NetaBins; i++){
    for(int j = 0; j < 2; j++){
      Q1_count[i][j] = 0.0;
      Q2_count[i][j] = 0.0;
      for(int k = 0; k < 2; k++){
        Q1[i][j][k] = 0.0;
        Q2[i][j][k] = 0.0;
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
        double nhits = trk.numberOfValidHits();
        double chi2n = trk.normalizedChi2();
        double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
        chi2n = chi2n/nlayers;

        double weight = 1.0;

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if(fabs(trk.eta()) < 2.4 && trk.pt() > 0.4 ){nTracks++;}// NtrkOffline        
        if(fabs(trk.eta()) > 2.4 || trk.pt() < ptLow_ || trk.pt() > ptHigh_) continue;
        if(chi2n > offlineChi2_) continue;
        if(nhits < offlinenhits_) continue;
        if( messAcceptance_ ) { if( trk.phi() < holeRight_ && trk.phi() > holeLeft_ ) continue;}
        if( doEffCorrection_ ){ weight = 1.0/effTable->GetBinContent( effTable->FindBin(trk.eta(), trk.pt()) );}
       
        trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off

        QcosTRK += weight*cos( 2*trk.phi() );
        QsinTRK += weight*sin( 2*trk.phi() );
        QcountsTrk += weight;

        for(int eta = 0; eta < NetaBins; eta++){
          if( trk.eta() > etaBins_[eta] && trk.eta() < etaBins_[eta+1] ){

            if( trk.charge() == 1){
              Q1[eta][0][0] += weight*cos( trk.phi() );
              Q1[eta][0][1] += weight*sin( trk.phi() );
              Q1_count[eta][0] += weight;

              Q2[eta][0][0] += weight*cos( 2*trk.phi() );
              Q2[eta][0][1] += weight*sin( 2*trk.phi() );
              Q2_count[eta][0] += weight;

            }
            else if( trk.charge() == -1){
              Q1[eta][1][0] += weight*cos( trk.phi() );
              Q1[eta][1][1] += weight*sin( trk.phi() );
              Q1_count[eta][1] += weight;

              Q2[eta][1][0] += weight*cos( 2*trk.phi() );
              Q2[eta][1][1] += weight*sin( 2*trk.phi() );
              Q2_count[eta][1] += weight;

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

  int HFside = 2;
  if( useBothSide_ ) HFside = 1;

  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
        
        if( reverseBeam_ ) caloEta = -hit.eta();
        if( messAcceptance_ ){if( caloPhi < holeRight_ && caloPhi > holeLeft_ ) continue;} hfPhi->Fill( caloPhi );//make sure if messAcceptance is on or off
        
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
    if( ieta > etaBins_[ieta] && ieta < etaBins_[ieta+1] ){  
      for(int HF = 0; HF < HFside; HF++){
        for(int sign = 0; sign < 2; sign++){
          if( Q1_count[ieta][sign] == 0.0 || ETT[HF] == 0.0 ) continue;

            double Q_real = get3RealOverlap(Q1[ieta][sign][0], Q2[ieta][sign][0], Q3[HF][0], Q1[ieta][sign][1], Q2[ieta][sign][1], Q3[HF][1], Q1_count[ieta][sign], ETT[HF] );
            QvsdEta[ieta][sign][HF]->Fill( Q_real, Q1_count[ieta][sign]*(Q1_count[ieta][sign]-1)*ETT[HF] );
        
          } 
        if( Q1_count[ieta][0] == 0 || Q1_count[ieta][1] == 0 || ETT[HF] == 0.0 ) continue;

        double Q_real = get3Real(Q1[ieta][0][0]/Q1_count[ieta][0],Q1[jeta][1][0]/Q1_count[jeta][1],Q3[HF][0]/ETT[HF], Q1[ieta][0][1]/Q1_count[ieta][0], Q1[jeta][1][1]/Q1_count[jeta][1], Q3[HF][1]/ETT[HF]);
        QvsdEta[ieta][2][HF]->Fill( Q_real, Q1_count[ieta][0]*Q1_count[jeta][1]*ETT[HF] );  
      } 
    }
  }

/*
calculate v2 using 3 sub-events method:
 */

  aveQ3[0][0]->Fill( Q3[0][0]/ETT[0], ETT[0] );//HF+ cos
  aveQ3[0][1]->Fill( Q3[0][1]/ETT[0], ETT[0] );//HF+ sin
  
  aveQ3[1][0]->Fill( Q3[1][0]/ETT[1], ETT[1] );//HF- cos
  aveQ3[1][1]->Fill( Q3[1][1]/ETT[1], ETT[1] );//HF- sin

  QcosTRK = QcosTRK/QcountsTrk;
  QsinTRK = QsinTRK/QcountsTrk;

  double QaQc = get2Real(Q3[1][0]/ETT[1], QcosTRK/QcountsTrk, Q3[1][1]/ETT[1], QsinTRK/QcountsTrk );
  double QaQb = get2Real(Q3[1][0]/ETT[1], Q3[0][0]/ETT[0], Q3[1][1]/ETT[1], -Q3[0][1]/ETT[0]);//an extra minus sign 
  double QcQb = get2Real(QcosTRK/QcountsTrk, Q3[0][0]/ETT[0], QsinTRK/QcountsTrk, Q3[0][1]/ETT[0]);

  c2_ac->Fill( QaQc, ETT[1]*QcountsTrk );
  c2_cb->Fill( QcQb, ETT[1]*ETT[0]  );
  c2_ab->Fill( QaQb, ETT[0]*QcountsTrk );

}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorEtaTest::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  edm::FileInPath fip1("CMEandCorrelation/ThreePointCorrelator/data/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
  TFile f1(fip1.fullPath().c_str(),"READ");
  effTable = (TH2D*)f1.Get("rTotalEff3D");

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);

  const int NetaBins = etaBins_.size() - 1 ;
  int HFside = 2;
  if( useBothSide_ ) HFside = 1;
//HF:
  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 20000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 20000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 20000,-1,1);

  for(int i = 0; i < 2; i++ ){
      for(int j = 0; j < 2; j++){

        aveQ3[i][j] = fs->make<TH1D>(Form("aveQ3_%d_%d",i, j), ";aveQ3", 20000, -1.0, 1.0);
      }
  }

  for(int eta = 0; eta < NetaBins; eta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < HFside; HF++){       
        QvsdEta[eta][sign][HF] = fs->make<TH1D>(Form("QvsdEta_%d_%d_%d",eta,sign,HF), "", 20000,-1.0-0.00005, 1.0-0.00005);
      }
    }
  }
}


double 
ThreePointCorrelatorEtaTest::get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}

double ThreePointCorrelatorEtaTest::get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = N1*(N1-1)*N3;

      if( N == 0.0 ){return 0.0;}
      else{return (t1-t2)/N;}

}
double ThreePointCorrelatorEtaTest::get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double ThreePointCorrelatorEtaTest::get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = N1*(N1-1)*N3;

      if( N == 0.0 ){ return 0.0;}
      else{return (t1+t2)/N;}

}

double ThreePointCorrelatorEtaTest::get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double ThreePointCorrelatorEtaTest::get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double ThreePointCorrelatorEtaTest::get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double ThreePointCorrelatorEtaTest::get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
} 
// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorEtaTest::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorEtaTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorEtaTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorEtaTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorEtaTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorEtaTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorEtaTest);