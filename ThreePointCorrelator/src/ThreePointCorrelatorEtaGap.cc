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

#include "CMEandCorrelation/ThreePointCorrelator/interface/ThreePointCorrelatorBase.h"

ThreePointCorrelatorEtaGap::ThreePointCorrelatorEtaGap(const edm::ParameterSet& iConfig)

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
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  do3pTracker_ = iConfig.getUntrackedParameter<bool>("do3pTracker");
  
  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");
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
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");

}


ThreePointCorrelatorEtaGap::~ThreePointCorrelatorEtaGap()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

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
  if( fabs(bestvz) < vzLow_ || fabs(bestvz) > vzHigh_ ) return;

  vtxZ->Fill( bestvz );
  
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
  double dEtaBinsArray[100];

  for(unsigned i = 0; i < dEtaBins_.size(); i++){

    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
  }


// initialize Qcos and Qsin

  double QcosTRK = 0.;
  double QsinTRK = 0.;
  double QcountsTrk = 0;

  double Q1[NetaBins][2][2];
  double Q1_count[NetaBins][2];

  double Q2[NetaBins][2][2];
  double Q2_count[NetaBins][2];

  double P1[NetaBins][2][2];
  double P2[NetaBins][2][2];

  for(int i = 0; i < NetaBins; i++){
    for(int j = 0; j < 2; j++){
      Q1_count[i][j] = 0.0;
      Q2_count[i][j] = 0.0;
      for(int k = 0; k < 2; k++){
        Q1[i][j][k] = 0.0;
        Q2[i][j][k] = 0.0;

        P1[i][j][k] = 0.0;
        P2[i][j][k] = 0.0;
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
        if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
        if(fabs(dzvtx/dzerror) > 3.0) continue;
        if(fabs(dxyvtx/dxyerror) > 3.0) continue;
        if(fabs(trk.eta()) < 2.4 && trk.pt() > 0.4 ){nTracks++;}// NtrkOffline

  }

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];

     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        double nlayers = trk.hitPattern().pixelLayersWithMeasurement();//only pixel layers
        double chi2n = trk.normalizedChi2();
        double nlayersTracker = trk.hitPattern().trackerLayersWithMeasurement();
        chi2n = chi2n/nlayersTracker;
        double trkEta = trk.eta();

        double weight = 1.0;
        if( doEffCorrection_ ) { weight = 1.0/effTable->GetBinContent( effTable->FindBin(trk.eta(), trk.pt()) );}

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if(nlayers <= 0 ) continue;
        if(chi2n > offlineChi2_ ) continue;
        if(fabs(trk.eta()) > etaTracker_ || trk.pt() < ptLow_ || trk.pt() > ptHigh_) continue;
        if( messAcceptance_ ) { if( trk.phi() < holeRight_ && trk.phi() > holeLeft_ ) continue;}
        if( reverseBeam_ ) {trkEta = -trk.eta();}

        trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off
        trkPt->Fill( trk.pt(), weight);//single particle closure
        trk_eta->Fill( trk.eta(), weight);

        QcosTRK += weight*cos( 2*trk.phi() );
        QsinTRK += weight*sin( 2*trk.phi() );
        QcountsTrk += weight;

        for(int eta = 0; eta < NetaBins; eta++){
          if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

            if( trk.charge() == 1){

              Q1[eta][0][0] += weight*cos( trk.phi() );
              Q1[eta][0][1] += weight*sin( trk.phi() );

              P1[eta][0][0] += weight*cos( trk.phi() );//the same as Q1, just for consistency
              P1[eta][0][1] += weight*sin( trk.phi() );

              Q1_count[eta][0] += weight;

              Q2[eta][0][0] += weight*weight*cos( 2*trk.phi() );
              Q2[eta][0][1] += weight*weight*sin( 2*trk.phi() );

              P2[eta][0][0] += weight*cos( -trk.phi() );
              P2[eta][0][1] += weight*sin( -trk.phi() );

              Q2_count[eta][0] += weight*weight;

            }
            else if( trk.charge() == -1){
              
              Q1[eta][1][0] += weight*cos( trk.phi() );
              Q1[eta][1][1] += weight*sin( trk.phi() );

              P1[eta][1][0] += weight*cos( trk.phi() );//the same as Q1, just for consistency
              P1[eta][1][1] += weight*sin( trk.phi() );
            
              Q1_count[eta][1] += weight;

              Q2[eta][1][0] += weight*weight*cos( 2*trk.phi() );
              Q2[eta][1][1] += weight*weight*sin( 2*trk.phi() );

              P2[eta][1][0] += weight*cos( -trk.phi() );
              P2[eta][1][1] += weight*sin( -trk.phi() );
              
              Q2_count[eta][1] += weight*weight;

            }
          }
        }     
  } 

//loop over calo towers (HF)

  double Q3[2][2];
  double ETT[2];

  for(int i = 0; i < 2; i++){
    ETT[i] = 0.;
    for(int j = 0; j < 2; j++){
      Q3[i][j] = 0.;
    }
  }


  if( do3pTracker_ ){

    for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];

     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        double nlayers = trk.hitPattern().pixelLayersWithMeasurement();//only pixel layers
        double trkEta = trk.eta();

        double weight = 1.0;

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if( trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
        if(nlayers <= 0 ) continue;
        if( messAcceptance_ ) { if( trk.phi() < holeRight_ && trk.phi() > holeLeft_ ) continue;}
        if( doEffCorrection_ ) { weight = 1.0/effTable->GetBinContent( effTable->FindBin(trk.eta(), trk.pt()) );}
        if( reverseBeam_ ) {trkEta = -trk.eta();}

        if( trkEta > etaTracker_ && trkEta < 2.4 ){

          Q3[0][0] += weight*cos( -2*trk.phi() );
          Q3[0][1] += weight*sin( -2*trk.phi() );
          ETT[0] += weight;
        }
        else if( trkEta > -2.4 && trkEta < -etaTracker_ ){
              
          Q3[1][0] += weight*cos( -2*trk.phi() );
          Q3[1][1] += weight*sin( -2*trk.phi() );
          ETT[1] += weight;
        }
        else{continue;}
      }

  }
  else{

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
  }

//randomly rotate the reaction plane:

  double num1 = fRand(-3.14,3.14);
  double num2 = fRand(-3.14,3.14);

  double constant1 = Q3[0][0]*Q3[0][0] + Q3[0][1]*Q3[0][1];
  double d1 = 1 + tan(num1)*tan(num1);

  double x1 = sqrt( constant1/d1 );
  double y1 = tan(num1) * x1;

  double constant2 = Q3[1][0]*Q3[1][0] + Q3[1][1]*Q3[1][1];
  double d2 = 1 + tan(num2)*tan(num2);

  double x2 = sqrt( constant2/d2 );
  double y2 = tan(num2) * x2;

  Q3[0][0] = x1;
  Q3[0][1] = y1;

  Q3[1][0] = x2;
  Q3[1][1] = y2;

//2p correlators
  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){

      double deltaEta = fabs(etaBins_[ieta] - etaBins_[jeta]);
      
      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){
          if( deta == 0){
            for(int sign = 0; sign < 2; sign++){
              if( Q1_count[ieta][sign] == 0.0 ) continue;

              delEta2p[sign]->Fill( deltaEta );
              double P_real = get2Real(P1[ieta][sign][0], P2[ieta][sign][0], P1[ieta][sign][1], P2[ieta][sign][1]);
              double P_real_count = Q1_count[ieta][sign]*Q1_count[ieta][sign] - Q2_count[ieta][sign];
              P_real = (P_real - Q2_count[ieta][sign])/P_real_count; //for COS(P1-P2) needs to minus the N.

              PvsdEta[deta][sign]->Fill( P_real, P_real_count  );
            }

            if( Q1_count[ieta][0] == 0.0 || Q1_count[ieta][1] == 0.0 ) continue;

              delEta2p[2]->Fill( deltaEta );
              double P_real = get2Real(P1[ieta][0][0], P2[ieta][1][0], P1[ieta][0][1], P2[ieta][1][1]);
              double P_real_count = Q1_count[ieta][0]*Q1_count[ieta][1];
              P_real = P_real/P_real_count;

              PvsdEta[deta][2]->Fill( P_real, P_real_count  );

          }
          else{
            for(int sign = 0; sign < 2; sign++){
              if( Q1_count[ieta][sign] == 0.0 || Q1_count[jeta][sign] == 0.0 ) continue;
             
              delEta2p[sign]->Fill( deltaEta );
              double P_real = get2Real(P1[ieta][sign][0], P2[jeta][sign][0], P1[ieta][sign][1], P2[jeta][sign][1]);
              double P_real_count = Q1_count[ieta][sign]*Q1_count[jeta][sign];
              P_real = P_real/P_real_count; 

              PvsdEta[deta][sign]->Fill( P_real, P_real_count  );
            }

            if( Q1_count[ieta][0] == 0.0 || Q1_count[jeta][1] == 0.0 ) continue;

              delEta2p[2]->Fill( deltaEta );
              double P_real = get2Real(P1[ieta][0][0], P2[jeta][1][0], P1[ieta][0][1], P2[jeta][1][1]);
              double P_real_count = Q1_count[ieta][0]*Q1_count[jeta][1];
              P_real = P_real/P_real_count;

              PvsdEta[deta][2]->Fill( P_real, P_real_count  );

          }
        }
      }
    }
  }

  int HFside = 2;
  if( useBothSide_ ) HFside = 1;

//3p correlators
  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){
    
      double deltaEta = fabs(etaBins_[ieta] - etaBins_[jeta]);

      for(int deta = 0; deta < NdEtaBins; deta++){
        if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){
          if( deta == 0 ){
            for(int HF = 0; HF < HFside; HF++){
              for(int sign = 0; sign < 2; sign++){

                if( Q1_count[ieta][sign] == 0.0 || ETT[HF] == 0.0 ) continue;

                delEta3p[sign][HF]->Fill( deltaEta );
                double Q_real = get3RealOverlap(Q1[ieta][sign][0], Q2[ieta][sign][0], Q3[HF][0], Q1[ieta][sign][1], Q2[ieta][sign][1], Q3[HF][1], Q1_count[ieta][sign], Q2_count[ieta][sign], ETT[HF] );
                QvsdEta[deta][sign][HF]->Fill( Q_real, (Q1_count[ieta][sign]*Q1_count[ieta][sign] - Q2_count[ieta][sign])*ETT[HF] );
              }

              if( Q1_count[ieta][0] == 0 || Q1_count[ieta][1] == 0 || ETT[HF] == 0.0 ) continue;

              delEta3p[2][HF]->Fill( deltaEta );
              double Q_real = get3Real(Q1[ieta][0][0]/Q1_count[ieta][0],Q1[jeta][1][0]/Q1_count[jeta][1],Q3[HF][0]/ETT[HF], Q1[ieta][0][1]/Q1_count[ieta][0], Q1[jeta][1][1]/Q1_count[jeta][1], Q3[HF][1]/ETT[HF]);
              QvsdEta[deta][2][HF]->Fill( Q_real, Q1_count[ieta][0]*Q1_count[jeta][1]*ETT[HF] );  
            }
          }
          else{
            for(int HF = 0; HF < HFside; HF++){
              for(int sign = 0; sign < 2; sign++ ){
              
                if( Q1_count[ieta][sign] == 0.0 || Q1_count[jeta][sign] == 0.0 || ETT[HF] == 0.0 ) continue; 

                delEta3p[sign][HF]->Fill( deltaEta );
                double Q_real = get3Real(Q1[ieta][sign][0]/Q1_count[ieta][sign],Q1[jeta][sign][0]/Q1_count[jeta][sign],Q3[HF][0]/ETT[HF], Q1[ieta][sign][1]/Q1_count[ieta][sign], Q1[jeta][sign][1]/Q1_count[jeta][sign], Q3[HF][1]/ETT[HF]);
                QvsdEta[deta][sign][HF]->Fill( Q_real, Q1_count[ieta][sign]*Q1_count[jeta][sign]*ETT[HF] );  

              }

              if( Q1_count[ieta][0] == 0.0 || Q1_count[jeta][1] == 0.0 || ETT[HF] == 0.0 ) continue;
                
                delEta3p[2][HF]->Fill( deltaEta );
                double Q_real = get3Real(Q1[ieta][0][0]/Q1_count[ieta][0],Q1[jeta][1][0]/Q1_count[jeta][1],Q3[HF][0]/ETT[HF], Q1[ieta][0][1]/Q1_count[ieta][0], Q1[jeta][1][1]/Q1_count[jeta][1], Q3[HF][1]/ETT[HF]);
                QvsdEta[deta][2][HF]->Fill( Q_real, Q1_count[ieta][0]*Q1_count[jeta][1]*ETT[HF] );  

            }
          }
        }
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

  double QaQc = get2Real(Q3[1][0]/ETT[1], QcosTRK/QcountsTrk, Q3[1][1]/ETT[1], QsinTRK/QcountsTrk );
  double QcQb = get2Real(QcosTRK/QcountsTrk, Q3[0][0]/ETT[0], QsinTRK/QcountsTrk, Q3[0][1]/ETT[0]);
  double QaQb = get2Real(Q3[1][0]/ETT[1], Q3[0][0]/ETT[0], Q3[1][1]/ETT[1], -Q3[0][1]/ETT[0]);//an extra minus sign 

  c2_ac->Fill( QaQc, ETT[1]*QcountsTrk );
  c2_cb->Fill( QcQb, ETT[0]*QcountsTrk  );
  c2_ab->Fill( QaQb, ETT[1]*ETT[0] );

}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorEtaGap::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH2D::SetDefaultSumw2();

  const int NdEtaBins = dEtaBins_.size() - 1;
  double dEtaBinsArray[100];

  for(unsigned i = 0; i < dEtaBins_.size(); i++){

    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
  }
  const int NetaBins = etaBins_.size() - 1;
  double etaBinsArray[100];

  for(unsigned i = 0; i < etaBins_.size(); i++){

    etaBinsArray[i] = etaBins_[i];
  }
  const int Nptbins = ptBins_.size() - 1;
  double ptBinsArray[100];

  for(unsigned i = 0; i < ptBins_.size(); i++){

    ptBinsArray[i] = ptBins_[i];
  }
  int HFside = 2;
  if( useBothSide_ ) HFside = 1;

  // edm::FileInPath fip1("CMEandCorrelation/ThreePointCorrelator/data/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
  // TFile f1(fip1.fullPath().c_str(),"READ");
  // effTable = (TH2D*)f1.Get("rTotalEff3D");

  edm::FileInPath fip1("CMEandCorrelation/ThreePointCorrelator/data/EPOS_eff.root");  
  TFile f1(fip1.fullPath().c_str(),"READ");
  effTable = (TH2D*)f1.Get("recoHist");

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", NetaBins, etaBinsArray);

  for(int sign = 0; sign < 3; sign++){

    delEta2p[sign] = fs->make<TH1D>(Form("delEta2p_%d",sign),";#Delta#eta", NdEtaBins, dEtaBinsArray);

    for(int HF = 0; HF < HFside; HF++){
        delEta3p[sign][HF] = fs->make<TH1D>(Form("delEta3p_%d_%d",sign,HF),";#Delta#eta", NdEtaBins, dEtaBinsArray);

    }
  }
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
    for(int sign = 0; sign < 3; sign++){

      PvsdEta[deta][sign] = fs->make<TH1D>(Form("PvsdEta_%d_%d", deta, sign), "", 20000, -1.0-0.00005, 1.0-0.00005);

      for(int HF = 0; HF < HFside; HF++){       
        QvsdEta[deta][sign][HF] = fs->make<TH1D>(Form("QvsdEta_%d_%d_%d",deta,sign,HF), "", 20000,-1.0-0.00005, 1.0-0.00005);
      }
    }
  }
}


double 
ThreePointCorrelatorEtaGap::get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}

double ThreePointCorrelatorEtaGap::get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = (N1*N1-N2)*N3;

      if( N == 0.0 ){return 0.0;}
      else{return (t1-t2)/N;}

}
double ThreePointCorrelatorEtaGap::get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double ThreePointCorrelatorEtaGap::get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = (N1*N1-N2)*N3;

      if( N == 0.0 ){ return 0.0;}
      else{return (t1+t2)/N;}

}

double ThreePointCorrelatorEtaGap::get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double ThreePointCorrelatorEtaGap::get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double ThreePointCorrelatorEtaGap::get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double ThreePointCorrelatorEtaGap::get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
} 

double ThreePointCorrelatorEtaGap::fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
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