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

  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");

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

        }
        else if( trk.charge() == -1 ){

          Qcos[1] += cos( trk.phi() );
          Qsin[1] += sin( trk.phi() );

          Q2cos[1] += cos( 2*trk.phi() );
          Q2sin[1] += sin( 2*trk.phi() );

          Qcounts[1]++;

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

          HFqVcosMinus = HFqVcosMinus + w*cos( 2*caloPhi );
          HFqVsinMinus = HFqVsinMinus + w*sin( 2*caloPhi );
         
          HFminusCounts++; 
          ETTminus += w;

        }
        if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
          
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

  double QaQc = get2Real(HFqVcosMinus, QcosTRK, HFqVsinMinus, QsinTRK );
  double QaQb = get2Real(HFqVcosMinus, HFqVcosPlus, HFqVsinMinus, HFqVsinPlus);
  double QcQb = get2Real(QcosTRK, HFqVcosPlus, QsinTRK, HFqVsinPlus);

  c2_ac->Fill( QaQc );
  c2_cb->Fill( QcQb  );
  c2_ab->Fill( QaQb );





}
// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorNtrk::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);

  // double ntrkBinsFill[100];
  // const int nNtrkBins = ntrkBins_.size() - 1;
  // for(unsigned num = 0; num < ntrkBins_.size(); num++ ){

  //   ntrkBinsFill[num] = ntrkBins_[num];
  // }

//HF:

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 10000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 10000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 10000,-1,1);


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
