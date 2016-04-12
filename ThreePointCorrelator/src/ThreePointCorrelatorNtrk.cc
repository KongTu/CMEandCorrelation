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

#include "CMEandCorrelation/ThreePointCorrelator/interface/ThreePointCorrelatorBase.h"


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
  reverseBeam_ = iConfig.getUntrackedParameter<bool>("reverseBeam");
  messAcceptance_ = iConfig.getUntrackedParameter<bool>("messAcceptance");

  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);

  holesize_ = iConfig.getUntrackedParameter<double>("holesize");

  ntrkBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ntrkBins");

}


ThreePointCorrelatorNtrk::~ThreePointCorrelatorNtrk()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

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

  const int nNtrkBins = ntrkBins_.size() - 1;

  double QcosTRK = 0.;
  double QsinTRK = 0.;
  int QcountsTrk = 0;

  double Q1[2][2];
  double Q2[2][2];
  double Q3[2][2];
  double Q1_count[2]; 
  double ETT[2];

  for(int i = 0; i < 2; i++){

    Q1_count[i] = 0.;
    ETT[i] = 0.;
    for(int j = 0; j < 2; j++){

      Q1[i][j] = 0.;
      Q2[i][j] = 0.;
      Q3[i][j] = 0.;
    }
  }

  int HFside = 2;
  if( useCentrality_ ) HFside = 1;

//track loop
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
        if( messAcceptance_ ) { if( trk.phi() < 0.6 && trk.phi() > 0.0 ) continue;}trkPhi->Fill( trk.phi() );
        nTracks++;   

        //for c particle v2;
        QcosTRK += cos( 2*trk.phi() );
        QsinTRK += sin( 2*trk.phi() );
        QcountsTrk++;

        if( trk.charge() == 1 ){

          Q1[0][0] += cos( trk.phi() );
          Q1[0][1] += sin( trk.phi() );
          Q1_count[0]++;

          Q2[0][0] += cos( 2*trk.phi() );
          Q2[0][1] += sin( 2*trk.phi() );

        }
        else if( trk.charge() == -1 ){

          Q1[1][0] += cos( trk.phi() );
          Q1[1][1] += sin( trk.phi() );
          Q1_count[1]++;

          Q2[1][0] += cos( 2*trk.phi() );
          Q2[1][1] += sin( 2*trk.phi() );

        }
        else{
          continue;
        }
  } 

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

//loop over calo towers (HF)

  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
        if( reverseBeam_ ) caloEta = -hit.eta();
        if( messAcceptance_ ){if( caloPhi < 0.6 && caloPhi > 0.0 ) continue;} hfPhi->Fill( caloPhi );//make sure if messAcceptance is on or off
        
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

  for(int trk = 0; trk < nNtrkBins; trk++){
    if( nTracks >= ntrkBins_[trk] && nTracks < ntrkBins_[trk+1] ){
        for(int HF = 0; HF < HFside; HF++){
          for(int sign = 0; sign < 2; sign++){
            if( Q1_count[sign] == 0 || ETT[HF] == 0 ) continue;
              double Q_real = 0.;
              Q_real = get3RealOverlap(Q1[sign][0], Q2[sign][0], Q3[HF][0], Q1[sign][1], Q2[sign][1], Q3[HF][1], Q1_count[sign],ETT[HF]);
              QvsNtrk[trk][sign][HF]->Fill( Q_real );

              double XY_real_temp = get2RealOverlap(Q1[sign][0],Q2[sign][0],Q1[sign][1],Q2[sign][1]);
              double XZ_real_temp = get2Real(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);
              double YZ_real_temp = get2Real(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);

              double XY_imag_temp = get2ImagOverlap(Q1[sign][0],Q2[sign][0],Q1[sign][1],Q2[sign][1]);
              double XZ_imag_temp = get2Imag(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);
              double YZ_imag_temp = get2Imag(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);

              double X_real_temp = Q1[sign][0];
              double Y_real_temp = Q1[sign][0];
              double Z_real_temp = Q3[HF][0];

              double X_imag_temp = Q1[sign][1];
              double Y_imag_temp = Q1[sign][1];
              double Z_imag_temp = Q3[HF][1];

              XY_real[trk][sign][HF]->Fill( XY_real_temp/(Q1_count[sign]*(Q1_count[sign]-1) ), (Q1_count[sign]*(Q1_count[sign]-1) ) );
              XZ_real[trk][sign][HF]->Fill( XZ_real_temp/(Q1_count[sign]*ETT[HF]), (Q1_count[sign]*ETT[HF]));
              YZ_real[trk][sign][HF]->Fill( YZ_real_temp/(Q1_count[sign]*ETT[HF]), (Q1_count[sign]*ETT[HF]));

              XY_imag[trk][sign][HF]->Fill( XY_imag_temp/(Q1_count[sign]*(Q1_count[sign]-1) ), (Q1_count[sign]*(Q1_count[sign]-1) ) );
              XZ_imag[trk][sign][HF]->Fill( XZ_imag_temp/(Q1_count[sign]*ETT[HF]), (Q1_count[sign]*ETT[HF]));
              YZ_imag[trk][sign][HF]->Fill( YZ_imag_temp/(Q1_count[sign]*ETT[HF]), (Q1_count[sign]*ETT[HF]));

              X_real[trk][sign][HF]->Fill( X_real_temp/Q1_count[sign], Q1_count[sign]);
              Y_real[trk][sign][HF]->Fill( Y_real_temp/Q1_count[sign], Q1_count[sign]);
              Z_real[trk][sign][HF]->Fill( Z_real_temp/ETT[HF], ETT[HF]);
              
              X_imag[trk][sign][HF]->Fill( X_imag_temp/Q1_count[sign], Q1_count[sign]);
              Y_imag[trk][sign][HF]->Fill( Y_imag_temp/Q1_count[sign], Q1_count[sign]);
              Z_imag[trk][sign][HF]->Fill( Z_imag_temp/ETT[HF], ETT[HF]);
            } 
            
            if( Q1_count[0] == 0 || Q1_count[1] == 0 || ETT[HF] == 0 ) continue;
            double Q_real = 0.;
            Q_real = get3Real(Q1[0][0]/Q1_count[0], Q1[1][0]/Q1_count[1], Q3[HF][0]/ETT[HF], Q1[0][1]/Q1_count[0], Q2[1][1]/Q1_count[1], Q3[HF][1]/ETT[HF]);
            QvsNtrk[trk][2][HF]->Fill( Q_real );

            double XY_real_temp = get2Real(Q1[0][0],Q1[1][0],Q1[0][1],Q1[1][1]);
            double XZ_real_temp = get2Real(Q1[0][0],Q3[HF][0],Q1[0][1],Q3[HF][1]);
            double YZ_real_temp = get2Real(Q1[1][0],Q3[HF][0],Q1[1][1],Q3[HF][1]);

            double XY_imag_temp = get2Imag(Q1[0][0],Q1[1][0],Q1[0][1],Q1[1][1]);
            double XZ_imag_temp = get2Imag(Q1[0][0],Q3[HF][0],Q1[0][1],Q3[HF][1]);
            double YZ_imag_temp = get2Imag(Q1[1][0],Q3[HF][0],Q1[1][1],Q3[HF][1]);

            double X_real_temp = Q1[0][0];
            double Y_real_temp = Q1[1][0];
            double Z_real_temp = Q3[HF][0];

            double X_imag_temp = Q1[0][1];
            double Y_imag_temp = Q1[1][1];
            double Z_imag_temp = Q3[HF][1];

            XY_real[trk][2][HF]->Fill( XY_real_temp/(Q1_count[0]*Q1_count[1]), (Q1_count[0]*Q1_count[1])  );
            XZ_real[trk][2][HF]->Fill( XZ_real_temp/(Q1_count[0]*ETT[HF]), (Q1_count[0]*ETT[HF]));
            YZ_real[trk][2][HF]->Fill( YZ_real_temp/(Q1_count[1]*ETT[HF]), (Q1_count[1]*ETT[HF]));

            XY_imag[trk][2][HF]->Fill( XY_imag_temp/(Q1_count[0]*Q1_count[1]), (Q1_count[0]*Q1_count[1])  );
            XZ_imag[trk][2][HF]->Fill( XZ_imag_temp/(Q1_count[0]*ETT[HF]), (Q1_count[0]*ETT[HF]));
            YZ_imag[trk][2][HF]->Fill( YZ_imag_temp/(Q1_count[1]*ETT[HF]), (Q1_count[1]*ETT[HF]));

            X_real[trk][2][HF]->Fill( X_real_temp/Q1_count[0], Q1_count[0]);
            Y_real[trk][2][HF]->Fill( Y_real_temp/Q1_count[1], Q1_count[1]);
            Z_real[trk][2][HF]->Fill( Z_real_temp/ETT[HF], ETT[HF]);
            
            X_imag[trk][2][HF]->Fill( X_imag_temp/Q1_count[0], Q1_count[0]);
            Y_imag[trk][2][HF]->Fill( Y_imag_temp/Q1_count[1], Q1_count[1]);
            Z_imag[trk][2][HF]->Fill( Z_imag_temp/ETT[HF], ETT[HF]);
      }       
    }
  }


/*
calculate v2 using 3 sub-events method:
 */

  aveQ3[0][0]->Fill( Q3[0][0]/ETT[0] );//HF+ cos
  aveQ3[0][1]->Fill( Q3[0][1]/ETT[0] );//HF+ sin
  
  aveQ3[1][0]->Fill( Q3[1][0]/ETT[1] );//HF- cos
  aveQ3[1][1]->Fill( Q3[1][1]/ETT[1] );//HF- sin

  QcosTRK = QcosTRK/QcountsTrk;
  QsinTRK = QsinTRK/QcountsTrk;

  double QaQc = get2Real(Q3[1][0]/ETT[1], QcosTRK/QcountsTrk, Q3[1][1]/ETT[1], QsinTRK/QcountsTrk );
  double QaQb = get2Real(Q3[1][0]/ETT[1], Q3[0][0]/ETT[0], Q3[1][1]/ETT[1], -Q3[0][1]/ETT[0]);//an extra minus sign 
  double QcQb = get2Real(QcosTRK/QcountsTrk, Q3[0][0]/ETT[0], QsinTRK/QcountsTrk, Q3[0][1]/ETT[0]);

  c2_ac->Fill( QaQc );
  c2_cb->Fill( QcQb  );
  c2_ab->Fill( QaQb );

}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorNtrk::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();
  
  const int nNtrkBins = ntrkBins_.size() - 1;

  int HFside = 2;
  if( useCentrality_ ) HFside = 1;

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 10000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 10000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 10000,-1,1);
  
  for(int i = 0; i < 2; i++ ){
      for(int j = 0; j < 2; j++){

        aveQ3[i][j] = fs->make<TH1D>(Form("aveQ3_%d_%d",i, j), ";aveQ3", 20000, -1.0, 1.0);
      }
  }
  for(int trk = 0; trk < nNtrkBins; trk++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < HFside; HF++){       
        QvsNtrk[trk][sign][HF] = fs->make<TH1D>(Form("QvsNtrk_%d_%d_%d",trk,sign,HF), "", 20000,-1.0-0.00005, 1.0-0.00005);
        
        XY_real[trk][sign][HF] = fs->make<TH1D>(Form("XY_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        XZ_real[trk][sign][HF] = fs->make<TH1D>(Form("XZ_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        YZ_real[trk][sign][HF] = fs->make<TH1D>(Form("YZ_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        
        XY_imag[trk][sign][HF] = fs->make<TH1D>(Form("XY_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        XZ_imag[trk][sign][HF] = fs->make<TH1D>(Form("XZ_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        YZ_imag[trk][sign][HF] = fs->make<TH1D>(Form("YZ_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        
        X_real[trk][sign][HF] = fs->make<TH1D>(Form("X_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        Y_real[trk][sign][HF] = fs->make<TH1D>(Form("Y_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        Z_real[trk][sign][HF] = fs->make<TH1D>(Form("Z_real_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        
        X_imag[trk][sign][HF] = fs->make<TH1D>(Form("X_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        Y_imag[trk][sign][HF] = fs->make<TH1D>(Form("Y_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
        Z_imag[trk][sign][HF] = fs->make<TH1D>(Form("Z_imag_%d_%d_%d",trk,sign,HF), "", 20000,-1.0,1.0);
      }
    }
  }

}
//
double 
ThreePointCorrelatorNtrk::get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}

double ThreePointCorrelatorNtrk::get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = N1*(N1-1)*N3;

      return (t1-t2)/N;
}
double ThreePointCorrelatorNtrk::get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double ThreePointCorrelatorNtrk::get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = N1*(N1-1)*N3;

      return (t1+t2)/N;
}

double ThreePointCorrelatorNtrk::get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double ThreePointCorrelatorNtrk::get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double ThreePointCorrelatorNtrk::get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double ThreePointCorrelatorNtrk::get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
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
