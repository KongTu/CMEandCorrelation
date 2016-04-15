// -*- C++ -*-
//
// Package:    ThreePointCorrelatorNestedLoop
// Class:      ThreePointCorrelatorNestedLoop
// 
/**\class ThreePointCorrelatorNestedLoop ThreePointCorrelatorNestedLoop.cc CMEandCorrelation/ThreePointCorrelatorNestedLoop/src/ThreePointCorrelatorNestedLoop.cc

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

ThreePointCorrelatorNestedLoop::ThreePointCorrelatorNestedLoop(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");

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

  holesize_ = iConfig.getUntrackedParameter<double>("holesize");

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");

}


ThreePointCorrelatorNestedLoop::~ThreePointCorrelatorNestedLoop()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called for each event  ------------
void
ThreePointCorrelatorNestedLoop::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  double dEtaBinsArray[48];

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

  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByLabel(genParticleSrc_, genParticleCollection);

  double real_term[48];
  int Npairs[48];

  for(int eta = 0; eta < NdEtaBins; eta++){
    real_term[eta] = 0.;
    Npairs[eta] = 0;
  }

  for(unsigned it=0; it<genParticleCollection->size(); ++it) {

    const reco::GenParticle & genCand1 = (*genParticleCollection)[it];
    int status1 = genCand1.status();
    double genpt1 = genCand1.pt();
    double geneta1 = genCand1.eta();
    double genphi1 = genCand1.phi();
    int gencharge1 = genCand1.charge();

    if( status1 != 1  || gencharge1 != 1 ) continue;
    if( genpt1 < ptLow_ || genpt1 > ptHigh_ ) continue;

    cout << "loop 1 with " << it << endl;

      for(unsigned jt=0; jt<genParticleCollection->size(); ++jt) {

        const reco::GenParticle & genCand2 = (*genParticleCollection)[jt];
        int status2 = genCand2.status();
        double genpt2 = genCand2.pt();
        double geneta2 = genCand2.eta();
        double genphi2 = genCand2.phi();
        int gencharge2 = genCand2.charge();

        if( status2 != 1  || gencharge2 != 1 ) continue;//only plus sign
        if( genpt2 < ptLow_ || genpt2 > ptHigh_ ) continue;
        //if( it == jt ) continue;

        cout << "loop 2 with " << jt << endl;


          for(unsigned kt=0; kt<genParticleCollection->size(); ++kt) {

            const reco::GenParticle & genCand3 = (*genParticleCollection)[kt];
            int status3 = genCand3.status();
            double geneta3 = genCand3.eta();
            double genphi3 = genCand3.phi();
            int gencharge3 = genCand3.charge();

            if( status3 != 1  || gencharge3 == 0 ) continue;
            if( geneta3 < etaLowHF_ || geneta3 > etaHighHF_ ) continue;


            cout << "loop 3 with " << kt << endl;

            double deltaEta = fabs(geneta1 - geneta2);
              for(int deta = 0; deta < NdEtaBins; deta++){
                if( deltaEta > dEtaBinsArray[deta] && deltaEta < dEtaBinsArray[deta]  ){

                  real_term[deta] += cos( genphi1 + genphi2 - 2*genphi3 );
                  cout << "deta = " << deta << " , cos(p1+p2-2p3) = " << real_term[deta] << endl;
                  Npairs[deta]++;

                }
              }
          }
      }
  }


  for(int deta = 0; deta < NdEtaBins; deta++){

    QvsdEta[deta]->Fill( real_term[deta]/Npairs[deta], Npairs[deta]);

  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorNestedLoop::beginJob()
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

  const int NdEtaBins = dEtaBins_.size() - 1;


  for(int deta = 0; deta < NdEtaBins; deta++){
      
    QvsdEta[deta] = fs->make<TH1D>(Form("QvsdEta_%d",deta), "", 20000,-1.0-0.00005,1.0-0.00005);
  
  }

}
double 
ThreePointCorrelatorNestedLoop::get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}

double ThreePointCorrelatorNestedLoop::get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = N1*(N1-1)*N3;

      return (t1-t2)/N;
}
double ThreePointCorrelatorNestedLoop::get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double ThreePointCorrelatorNestedLoop::get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = N1*(N1-1)*N3;

      return (t1+t2)/N;
}

double ThreePointCorrelatorNestedLoop::get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double ThreePointCorrelatorNestedLoop::get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double ThreePointCorrelatorNestedLoop::get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double ThreePointCorrelatorNestedLoop::get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
} 
// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelatorNestedLoop::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelatorNestedLoop::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelatorNestedLoop::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelatorNestedLoop::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelatorNestedLoop::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelatorNestedLoop::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelatorNestedLoop);
