// -*- C++ -*-
//
// Package:    ThreePointCorrelatorEtaGapNestedLoop
// Class:      ThreePointCorrelatorEtaGapNestedLoop
// 
/**\class ThreePointCorrelatorEtaGapNestedLoop ThreePointCorrelatorEtaGapNestedLoop.cc CMEandCorrelation/ThreePointCorrelatorEtaGapNestedLoop/src/ThreePointCorrelatorEtaGapNestedLoop.cc

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

ThreePointCorrelatorEtaGapNestedLoop::ThreePointCorrelatorEtaGapNestedLoop(const edm::ParameterSet& iConfig)

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


ThreePointCorrelatorEtaGapNestedLoop::~ThreePointCorrelatorEtaGapNestedLoop()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
ThreePointCorrelatorEtaGapNestedLoop::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  const int NdEtaBins = dEtaBins_.size() - 1;
  double dEtaBinsArray[100];

  for(unsigned i = 0; i < dEtaBins_.size(); i++){

    dEtaBinsArray[i] = dEtaBins_[i]-0.0001;
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
        if(fabs(trk.eta()) < 2.4 && trk.pt() > 0.4 ){nTracks++;}// NtrkOffline        

  } 

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;
  
  Ntrk->Fill(nTracks);

  double QcosTRK = 0.;
  double QsinTRK = 0.;
  double QcountsTrk = 0;

  double Q3[2][2];
  double ETT[2];

  for(int i = 0; i < 2; i++){
    ETT[i] = 0.;
    for(int j = 0; j < 2; j++){
      Q3[i][j] = 0.;
    }
  }

  double real_term[48][3][2];
  double Npairs[48][3][2];

  for(int eta = 0; eta < NdEtaBins; eta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){
        
        real_term[eta][sign][HF] = 0.;
        Npairs[eta][sign][HF] = 0.;
      }
    } 
  }

  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];

     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        double nlayers = trk.hitPattern().pixelLayersWithMeasurement();//only pixel layers

        double weight = 1.0;
        if( doEffCorrection_ ) { weight = 1.0/effTable->GetBinContent( effTable->FindBin(trk.eta(), trk.pt()) );}

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if(fabs(trk.eta()) > 2.4 || trk.pt() < ptLow_ || trk.pt() > ptHigh_) continue;
        if(nlayers <= 0 ) continue;

        trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off
        trkPt->Fill( trk.pt(), weight);//single particle closure
        trk_eta->Fill( trk.eta(), weight);

        if( !do3pTracker_ ){
          QcosTRK += weight*cos( 2*trk.phi() );
          QsinTRK += weight*sin( 2*trk.phi() );
          QcountsTrk += weight;
        }

        for(unsigned jt = 0; jt < tracks->size(); jt++){

          const reco::Track & trk1 = (*tracks)[jt];

          math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

          double dzvtx = trk1.dz(bestvtx);
          double dxyvtx = trk1.dxy(bestvtx);
          double dzerror = sqrt(trk1.dzError()*trk1.dzError()+bestvzError*bestvzError);
          double dxyerror = sqrt(trk1.d0Error()*trk1.d0Error()+bestvxError*bestvyError);
          double nlayers = trk1.hitPattern().pixelLayersWithMeasurement();//only pixel layers

          double weight1 = 1.0;
          if( doEffCorrection_ ) { weight1 = 1.0/effTable->GetBinContent( effTable->FindBin(trk1.eta(), trk1.pt()) );}

          if(!trk1.quality(reco::TrackBase::highPurity)) continue;
          if(fabs(trk1.ptError())/trk1.pt() > offlineptErr_ ) continue;
          if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
          if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
          if(fabs(trk1.eta()) > 2.4 || trk1.pt() < ptLow_ || trk1.pt() > ptHigh_) continue;
          if(nlayers <= 0 ) continue;

          if(it == jt ) continue;

          if( do3pTracker_ ){

            for(unsigned kt = 0; kt < tracks->size(); kt++){

            const reco::Track & trk2 = (*tracks)[kt];

            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

            double dzvtx = trk2.dz(bestvtx);
            double dxyvtx = trk2.dxy(bestvtx);
            double dzerror = sqrt(trk2.dzError()*trk2.dzError()+bestvzError*bestvzError);
            double dxyerror = sqrt(trk2.d0Error()*trk2.d0Error()+bestvxError*bestvyError);
            double nlayers = trk2.hitPattern().pixelLayersWithMeasurement();//only pixel layers

            double weight2 = 1.0;
            if( doEffCorrection_ ) { weight1 = 1.0/effTable->GetBinContent( effTable->FindBin(trk2.eta(), trk2.pt()) );}

            if(!trk2.quality(reco::TrackBase::highPurity)) continue;
            if(fabs(trk2.ptError())/trk2.pt() > offlineptErr_ ) continue;
            if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
            if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
            if(fabs(trk2.eta()) > 2.4) continue;
            if(nlayers <= 0 ) continue;

            double totalWeight = weight*weight1*weight2;

            if( fabs(trk2.eta()- trk.eta()) < 2.0 || fabs(trk2.eta() - trk1.eta()) < 2.0 ) continue;
            
            QcosTRK += weight2*cos( 2*trk2.phi() );
            QsinTRK += weight2*sin( 2*trk2.phi() );
            QcountsTrk += weight2;

              double deltaEta = fabs( trk.eta() - trk1.eta() );
              for(int deta = 0; deta < NdEtaBins; deta++){
                if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){

                    if( trk.charge() == 1 && trk1.charge() == 1){
                      
                      real_term[deta][0][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*trk2.phi() );
                      Npairs[deta][0][0] += totalWeight;
                    }
                    if( trk.charge() == -1 && trk1.charge() == -1){
                      real_term[deta][1][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*trk2.phi() );
                      Npairs[deta][1][0] += totalWeight;
                    }
                    if( trk.charge() == 1 && trk1.charge() == -1 ){

                      real_term[deta][2][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*trk2.phi() );
                      Npairs[deta][2][0] += totalWeight;
                    }
                    
                }
              }  
            }
          }
          else{
            for(unsigned i = 0; i < towers->size(); ++i){

              const CaloTower & hit= (*towers)[i];

              double caloEta = hit.eta();
              double caloPhi = hit.phi();
              double caloEt  = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );

              double totalWeight = weight*weight1*caloEt;
              
              if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
                double deltaEta = fabs( trk.eta() - trk1.eta() );
                for(int deta = 0; deta < NdEtaBins; deta++){
                  if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){

                      if( trk.charge() == 1 && trk1.charge() == 1){
                        
                        real_term[deta][0][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][0][0] += totalWeight;
                      }
                      if( trk.charge() == -1 && trk1.charge() == -1){
                        real_term[deta][1][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][1][0] += totalWeight;
                      }
                      if( trk.charge() == 1 && trk1.charge() == -1 ){

                        real_term[deta][2][0] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][2][0] += totalWeight;
                      }
                      
                  }
                }  
              }
              else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){
                double deltaEta = fabs( trk.eta() - trk1.eta() );
                for(int deta = 0; deta < NdEtaBins; deta++){
                  if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta+1] ){

                      if( trk.charge() == 1 && trk1.charge() == 1){
                        
                        real_term[deta][0][1] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][0][1] += totalWeight;
                      }
                      if( trk.charge() == -1 && trk1.charge() == -1){
                        real_term[deta][1][1] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][1][1] += totalWeight;
                      }
                      if( trk.charge() == 1 && trk1.charge() == -1 ){

                        real_term[deta][2][1] += totalWeight*cos( trk.phi() + trk1.phi() - 2*caloPhi );
                        Npairs[deta][2][1] += totalWeight;
                      }
                  }
                }  
              }
              else{continue;}
           }
          }
        }
  } 

  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){
        QvsdEta[deta][sign][HF]->Fill( real_term[deta][sign][HF]/Npairs[deta][sign][HF], Npairs[deta][sign][HF]);
      }
    }
  }

  for(unsigned i = 0; i < towers->size(); ++i){

    const CaloTower & hit= (*towers)[i];

    double caloEta = hit.eta();
    double caloPhi = hit.phi();
    double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
    
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

  double QaQc = get2Real(Q3[1][0]/ETT[1], QcosTRK/QcountsTrk, Q3[1][1]/ETT[1], QsinTRK/QcountsTrk );
  double QcQb = get2Real(QcosTRK/QcountsTrk, Q3[0][0]/ETT[0], QsinTRK/QcountsTrk, Q3[0][1]/ETT[0]);
  double QaQb = get2Real(Q3[1][0]/ETT[1], Q3[0][0]/ETT[0], Q3[1][1]/ETT[1], -Q3[0][1]/ETT[0]);//an extra minus sign 

  c2_ac->Fill( QaQc, ETT[1]*QcountsTrk );
  c2_cb->Fill( QcQb, ETT[0]*QcountsTrk  );
  c2_ab->Fill( QaQb, ETT[1]*ETT[0] );

}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NdEtaBins = dEtaBins_.size() - 1;

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

  // edm::FileInPath fip1("CMEandCorrelation/ThreePointCorrelator/data/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
  // TFile f1(fip1.fullPath().c_str(),"READ");
  // effTable = (TH2D*)f1.Get("rTotalEff3D");

  edm::FileInPath fip1("CMEandCorrelation/ThreePointCorrelator/data/EPOS_eff.root");  
  TFile f1(fip1.fullPath().c_str(),"READ");
  effTable = (TH2D*)f1.Get("recoHist");

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", NetaBins, etaBinsArray);

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 20000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 20000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 20000,-1,1);

  for(int deta = 0; deta < NdEtaBins; deta++){
    for(int sign = 0; sign < 3; sign++){
      for(int HF = 0; HF < 2; HF++){       
        QvsdEta[deta][sign][HF] = fs->make<TH1D>(Form("QvsdEta_%d_%d_%d",deta,sign,HF), "", 20000,-1.0-0.00005, 1.0-0.00005);
      }
    }
  }
 
}
double 
ThreePointCorrelatorEtaGapNestedLoop::get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}

double ThreePointCorrelatorEtaGapNestedLoop::get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = N1*(N1-1)*N3;

      if( N == 0.0 ){return 0.0;}
      else{return (t1-t2)/N;}

}
double ThreePointCorrelatorEtaGapNestedLoop::get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double ThreePointCorrelatorEtaGapNestedLoop::get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = N1*(N1-1)*N3;

      if( N == 0.0 ){ return 0.0;}
      else{return (t1+t2)/N;}

}

double ThreePointCorrelatorEtaGapNestedLoop::get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double ThreePointCorrelatorEtaGapNestedLoop::get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double ThreePointCorrelatorEtaGapNestedLoop::get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double ThreePointCorrelatorEtaGapNestedLoop::get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
} 
// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorEtaGapNestedLoop::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorEtaGapNestedLoop::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorEtaGapNestedLoop);