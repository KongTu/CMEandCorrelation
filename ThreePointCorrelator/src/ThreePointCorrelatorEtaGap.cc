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


//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* averageCosHF[2];//calculate the correction on v2
      TH1D* averageSinHF[2];
//end v2

      //culmulants: 
      TH2D* QvsdEtaPlusPlus[3];
      TH2D* QvsdEtaMinusMinus[3];
      TH2D* QvsdEtaPlusMinus[3];

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

      
      TH2D* EvsEta;
      TH2D* ETvsEta;

      TH1D* testDeta;
      TH1D* cbinHist;

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

  const int binSize_ = etaBins_.size() - 1 ;

// initialize Qcos and Qsin
  double Qcos[binSize_][2];
  double Qsin[binSize_][2];
  int Qcounts[binSize_][2];
  for(int eta = 0; eta < binSize_; eta++){
    for(int sign = 0; sign < 2; sign++){

      Qcounts[eta][sign] = 0;
      Qcos[eta][sign] = 0.;
      Qsin[eta][sign] = 0.;
    }
  }

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
        if(fabs(trk.eta()) > 2.4 || trk.pt() < 0.4) continue;
        if( messAcceptance_ ) {
          if( ( trk.phi() < (0.0 + holesize_) && trk.phi() > (0.0 - holesize_) )    ||
              ( trk.phi() < (2.09 + holesize_) && trk.phi() > (2.09 - holesize_) )  ||
              ( trk.phi() < (-2.09 + holesize_) && trk.phi() > (-2.09 - holesize_) )     ) continue;
        }
        nTracks++;  

        trkPhi->Fill( trk.phi() );//make sure if messAcceptance is on or off

        QcosTRK = QcosTRK + cos( 2*trk.phi() );
        QsinTRK = QsinTRK + sin( 2*trk.phi() );
        QcountsTrk++;

        for(int eta = 0; eta < binSize_; eta++){
          if( trk.eta() > etaBins_[eta] && trk.eta() < etaBins_[eta+1] ){

             if(trk.charge() == 1){

                Qcos[eta][0] = Qcos[eta][0] + cos( trk.phi() );
                Qsin[eta][0] = Qsin[eta][0] + sin( trk.phi() );
                Qcounts[eta][0]++;
             }
             if(trk.charge() == -1){
                
                Qcos[eta][1] = Qcos[eta][1] + cos( trk.phi() );
                Qsin[eta][1] = Qsin[eta][1] + sin( trk.phi() );
                Qcounts[eta][1]++;
             }
          }
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
        double energy = hit.emEnergy() + hit.hadEnergy();
        if( messAcceptance_ ){ 
          if( ( caloPhi < (0.0 + holesize_) && caloPhi > (0.0 - holesize_) )    ||
              ( caloPhi < (2.09 + holesize_) && caloPhi > (2.09 - holesize_) )  ||
              ( caloPhi < (-2.09 + holesize_) && caloPhi > (-2.09 - holesize_) )     ) continue;
        }

        hfPhi->Fill( caloPhi );//make sure if messAcceptance is on or off

        if( reverseBeam_ ) caloEta = -hit.eta();

        if( fabs(caloEta) < 5.0 && fabs(caloEta) > 3.0 ){

          EvsEta->Fill(caloEta, energy );
          ETvsEta->Fill(caloEta, w);
        }

        if( fabs(caloEta) > etaLowHF_ && fabs(caloEta) < etaHighHF_ ){

          HFqVcos = HFqVcos + w*cos( 2*caloPhi );
          HFqVsin = HFqVsin + w*sin( 2*caloPhi );
          
          HFcounts++;
          ETT += w;
        }
        if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

          HFqVcosMinus = HFqVcosMinus + w*cos( 2*caloPhi );
          HFqVsinMinus = HFqVsinMinus + w*sin( 2*caloPhi );
          
          HFcosSum[1]->Fill( w*cos( 2*caloPhi ) );
          HFsinSum[1]->Fill( w*sin( 2*caloPhi ) );
          weightSum[1]->Fill( w );

          HFminusCounts++; 
          ETTminus += w;

        }
        if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
          
          HFqVcosPlus = HFqVcosPlus + w*cos( 2*caloPhi );
          HFqVsinPlus = HFqVsinPlus + w*sin( 2*caloPhi );

          HFcosSum[0]->Fill( w*cos( 2*caloPhi ) );
          HFsinSum[0]->Fill( w*sin( 2*caloPhi ) );
          weightSum[0]->Fill( w );

          HFplusCounts++;
          ETTplus += w;
        }
  }

  for(int ieta = 0; ieta < binSize_; ieta++){
    for(int jeta = 0; jeta < binSize_; jeta++){
    
      if( ieta == jeta ) continue;

      double deltaEta = fabs(etaBins_[jeta] - etaBins_[ieta]);

      //QaQb (tracker)
      double pp_real = Qcos[ieta][0]*Qcos[jeta][0] - Qsin[ieta][0]*Qsin[jeta][0];
      double pp_imag = Qcos[ieta][0]*Qsin[jeta][0] + Qsin[ieta][0]*Qcos[jeta][0];

      double mm_real = Qcos[ieta][1]*Qcos[jeta][1] - Qsin[ieta][1]*Qsin[jeta][1];
      double mm_imag = Qcos[ieta][1]*Qsin[jeta][1] + Qsin[ieta][1]*Qcos[jeta][1];

      double pm_real = Qcos[ieta][0]*Qcos[jeta][1] - Qsin[ieta][0]*Qsin[jeta][1];
      double pm_imag = Qcos[ieta][0]*Qsin[jeta][1] + Qsin[ieta][0]*Qcos[jeta][1];

      int pp_count = Qcounts[ieta][0]*Qcounts[jeta][0];
      int mm_count = Qcounts[ieta][1]*Qcounts[jeta][1];
      int pm_count = Qcounts[ieta][0]*Qcounts[jeta][1];

      QaQbvsdEta[0][0]->Fill( deltaEta, pp_real );
      QaQbvsdEta[0][1]->Fill( deltaEta, pp_imag );

      QaQbvsdEta[1][0]->Fill( deltaEta, mm_real );
      QaQbvsdEta[1][1]->Fill( deltaEta, mm_imag );

      QaQbvsdEta[2][0]->Fill( deltaEta, pm_real );
      QaQbvsdEta[2][1]->Fill( deltaEta, pm_imag );

      NaNbvsdEta[0]->Fill( deltaEta, pp_count );
      NaNbvsdEta[1]->Fill( deltaEta, mm_count );
      NaNbvsdEta[2]->Fill( deltaEta, pm_count );


      //QaQc (tracker + HF+) and QbQc = QaQc
      for(int type = 0; type < 3; type++){

        double tempCosHF = 0.;
        double tempSinHF = 0.;
        double tempETT = 0.;

        if( type == 0 ){

          tempCosHF = HFqVcosPlus;
          tempSinHF = HFqVsinPlus;
          tempETT = ETTplus;
        }
        else if( type == 1 ){
          
          tempCosHF = HFqVcosMinus;
          tempSinHF = HFqVsinMinus;
          tempETT = ETTminus;

        }
        else if( type == 2 ){
          
          tempCosHF = HFqVcos;
          tempSinHF = HFqVsin;
          tempETT = ETT;

        }

        double pHF_real = Qcos[ieta][0]*tempCosHF + Qsin[ieta][0]*tempSinHF;
        double pHF_imag = -Qcos[ieta][0]*tempSinHF + Qsin[ieta][0]*tempCosHF;

        QaQcvsdEta[type][0][0]->Fill( deltaEta, pHF_real);
        QaQcvsdEta[type][0][1]->Fill( deltaEta, pHF_imag);

        QbvsdEta[type][0][0]->Fill( deltaEta, Qcos[jeta][0] );
        QbvsdEta[type][0][1]->Fill( deltaEta, Qsin[jeta][0] );

        NaNcvsdEta[type][0]->Fill( deltaEta, Qcounts[ieta][0]* tempETT );

        double mHF_real = Qcos[ieta][1]*tempCosHF + Qsin[ieta][1]*tempSinHF;
        double mHF_imag = -Qcos[ieta][1]*tempSinHF + Qsin[ieta][1]*tempCosHF;
        
        QaQcvsdEta[type][1][0]->Fill( deltaEta, mHF_real);
        QaQcvsdEta[type][1][1]->Fill( deltaEta, mHF_imag);

        QbvsdEta[type][1][0]->Fill( deltaEta, Qcos[jeta][1] );
        QbvsdEta[type][1][1]->Fill( deltaEta, Qsin[jeta][1] );

        NaNcvsdEta[type][1]->Fill( deltaEta, Qcounts[ieta][1]* tempETT );

        double pmHF_real = Qcos[ieta][0]*tempCosHF + Qsin[ieta][0]*tempSinHF;//this two terms are the same as pHF_real, just for format consistency
        double pmHF_imag = -Qcos[ieta][0]*tempSinHF + Qsin[ieta][0]*tempCosHF;
        
        QaQcvsdEta[type][2][0]->Fill( deltaEta, pmHF_real);
        QaQcvsdEta[type][2][1]->Fill( deltaEta, pmHF_imag);

        QbvsdEta[type][2][0]->Fill( deltaEta, Qcos[jeta][1] );
        QbvsdEta[type][2][1]->Fill( deltaEta, Qsin[jeta][1] );

        NaNcvsdEta[type][2]->Fill( deltaEta, Qcounts[ieta][0]* tempETT );// the same as the first one


      }

      //<Qa><Qb>:

      //plusplus:
      QaSinglevsdEta[0][0]->Fill( deltaEta, Qcos[ieta][0] );//plus
      QaSinglevsdEta[0][1]->Fill( deltaEta, Qsin[ieta][0] );

      QbSinglevsdEta[0][0]->Fill( deltaEta, Qcos[jeta][0] );//plus
      QbSinglevsdEta[0][1]->Fill( deltaEta, Qsin[jeta][0] );

      //minusminus
      QaSinglevsdEta[1][0]->Fill( deltaEta, Qcos[ieta][1] );//minus
      QaSinglevsdEta[1][1]->Fill( deltaEta, Qsin[ieta][1] );

      QbSinglevsdEta[1][0]->Fill( deltaEta, Qcos[jeta][1] );//minus
      QbSinglevsdEta[1][1]->Fill( deltaEta, Qsin[jeta][1] );

      //plusminus
      QaSinglevsdEta[2][0]->Fill( deltaEta, Qcos[ieta][0] );//plus
      QaSinglevsdEta[2][1]->Fill( deltaEta, Qsin[ieta][0] );

      QbSinglevsdEta[2][0]->Fill( deltaEta, Qcos[jeta][1] );//minus
      QbSinglevsdEta[2][1]->Fill( deltaEta, Qsin[jeta][1] );

    }
  }


//self normalize the Qvectors from HF:
  HFqVcosPlus = HFqVcosPlus/ETTplus;
  HFqVsinPlus = HFqVsinPlus/ETTplus;

  HFqVcosMinus = HFqVcosMinus/ETTminus;
  HFqVsinMinus = HFqVsinMinus/ETTminus;

  averageCosHF[0]->Fill( HFqVcosPlus );
  averageSinHF[0]->Fill( HFqVsinPlus );
  
  averageCosHF[1]->Fill( HFqVcosMinus );
  averageSinHF[1]->Fill( HFqVsinMinus );

  QcosTRK = QcosTRK/QcountsTrk;
  QsinTRK = QsinTRK/QcountsTrk;

  HFqVcos = HFqVcos/ETT;
  HFqVsin = HFqVsin/ETT;

  double QaQc = get2RealDup(HFqVcosMinus, QcosTRK, HFqVsinMinus, QsinTRK );
  double QaQb = get2RealDup(HFqVcosMinus, HFqVcosPlus, HFqVsinMinus, HFqVsinPlus);
  double QcQb = get2RealDup(QcosTRK, HFqVcosPlus, QsinTRK, HFqVsinPlus);

  c2_ac->Fill( QaQc );
  c2_cb->Fill( QcQb  );
  c2_ab->Fill( QaQb );

//3 point correlator
//self normalize the Qvectors;
 for(int eta = 0; eta < binSize_; eta++){
    for(int sign = 0; sign < 2; sign++){

      if( Qcounts[eta][sign] == 0 ) continue;

      Qcos[eta][sign] = Qcos[eta][sign]/Qcounts[eta][sign];
      Qsin[eta][sign] = Qsin[eta][sign]/Qcounts[eta][sign];
    }
  }

  for(int ieta = 0; ieta < binSize_; ieta++){
    for(int jeta = 0; jeta < binSize_; jeta++){

      if( ieta == jeta ) continue;
      
      double deltaEta = fabs(etaBins_[jeta] - etaBins_[ieta]);
      testDeta->Fill(deltaEta);
    
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
        else if ( type == 2 ){
          tempHFcos = HFqVcos;
          tempHFsin = HFqVsin;
        }
        else{
          return;
        }

        double totalQplusplus = get3RealDup(Qcos[ieta][0],Qcos[jeta][0], tempHFcos, Qsin[ieta][0], Qsin[jeta][0], tempHFsin );
        double totalQminusminus = get3RealDup(Qcos[ieta][1],Qcos[jeta][1], tempHFcos, Qsin[ieta][1], Qsin[jeta][1], tempHFsin );
        double totalQplusminus = get3RealDup(Qcos[ieta][0],Qcos[jeta][1], tempHFcos, Qsin[ieta][0], Qsin[jeta][1], tempHFsin );

        QvsdEtaPlusPlus[type]->Fill(deltaEta, totalQplusplus);
        QvsdEtaMinusMinus[type]->Fill(deltaEta, totalQminusminus);
        QvsdEtaPlusMinus[type]->Fill(deltaEta, totalQplusminus);

      }
    }
  }


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

  const int bins = dEtaBins_.size() - 1;
  const int temp = dEtaBins_.size();
  //const int NbinsEta = etaBins_.size() - 1;

  double dEtaBinsArray[48];
  for(int eta = 0; eta < temp; eta++){

    dEtaBinsArray[eta] = dEtaBins_[eta] - 0.0001;
  }
//HF:
  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 10000,-1,1);
  c2_ac = fs->make<TH1D>("c2_ac",";c2_ac", 10000,-1,1);
  c2_cb = fs->make<TH1D>("c2_cb",";c2_cb", 10000,-1,1);

  testDeta = fs->make<TH1D>("testDeta",";delta eta", bins, dEtaBinsArray);
  EvsEta = fs->make<TH2D>("EvsEta",";#eta;Energy(GeV)", 100, -5.0 , 5.0, 10000,0,500);
  ETvsEta = fs->make<TH2D>("ETvsEta",";#eta;E_{T}(GeV)", 100, -5.0 , 5.0, 10000,0,500);

  for(int sign = 0; sign < 2; sign++){

    HFsinSum[sign] = fs->make<TH1D>(Form("HFsinSum_%d", sign), ";HFsinSum", 2000, -1.0, 1.0 );
    HFcosSum[sign] = fs->make<TH1D>(Form("HFcosSum_%d", sign), ";HFcosSum", 2000, -1.0, 1.0 );
    weightSum[sign] = fs->make<TH1D>(Form("weightSum_%d", sign), ";weightSum", 3000, 0, 150 );

    averageCosHF[sign] = fs->make<TH1D>(Form("averageCosHF_%d", sign), ";averageCosHF", 2000, -1.0, 1.0);
    averageSinHF[sign] = fs->make<TH1D>(Form("averageSinHF_%d", sign), ";averageSinHF", 2000, -1.0, 1.0);
  }

  for(int type = 0; type < 3; type++){

    QvsdEtaPlusPlus[type] = fs->make<TH2D>(Form("QvsdEtaPlusPlus_%d", type),";#Delta#eta;Q_{#phi_{1,+}}Q_{#phi_{2,+}}Q^{*}_{2#phi_{3}}", bins, dEtaBinsArray, 40000,-2.0-0.00005,2.0-0.00005 );
    QvsdEtaMinusMinus[type] = fs->make<TH2D>(Form("QvsdEtaMinusMinus_%d", type),";#Delta#eta;Q_{#phi_{1,-}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", bins, dEtaBinsArray, 40000,-2.0-0.00005,2.0-0.00005 );
    QvsdEtaPlusMinus[type] = fs->make<TH2D>(Form("QvsdEtaPlusMinus_%d", type),";#Delta#eta;Q_{#phi_{1,+}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", bins, dEtaBinsArray, 40000,-2.0-0.00005,2.0-0.00005 );
    
    for(int sign = 0; sign < 3; sign++){

      NaNcvsdEta[type][sign] = fs->make<TH2D>(Form("NaNcvsdEta_%d_%d", type, sign), ";#Delta#eta;N_{a}N_{c}", bins, dEtaBinsArray, 3000,0,3000);
    }
  }

  for(int sign = 0; sign < 3; sign++){

      NaNbvsdEta[sign] = fs->make<TH2D>(Form("NaNbvsdEta_%d", sign), ";#Delta#eta;N_{a}N_{b}", bins, dEtaBinsArray, 500,0,500);

      for(int real = 0; real < 2; real++){

          QaQbvsdEta[sign][real] = fs->make<TH2D>(Form("QaQbvsdEta_%d_%d", sign, real), ";#Delta#eta;<Q_{a}Q_{b}>", bins, dEtaBinsArray, 40000,-200.0-0.005,200.0-0.005);
          QaSinglevsdEta[sign][real] = fs->make<TH2D>(Form("QaSinglevsdEta_%d_%d", sign, real), ";#Delta#eta;<Q_{a}>", bins, dEtaBinsArray, 20000,-10.0-0.0005,10.0-0.0005);
          QbSinglevsdEta[sign][real] = fs->make<TH2D>(Form("QbSinglevsdEta_%d_%d", sign, real), ";#Delta#eta;<Q_{b}>", bins, dEtaBinsArray, 20000,-10.0-0.0005,10.0-0.0005);

        for(int type = 0; type < 3; type++){
    
          QaQcvsdEta[type][sign][real] = fs->make<TH2D>(Form("QaQcvsdEta_%d_%d_%d", type, sign, real), ";#Delta#eta;<Q_{a}Q_{c}>", bins, dEtaBinsArray, 40000,-200.0-0.005,200.0-0.005);
          QbvsdEta[type][sign][real] = fs->make<TH2D>(Form("QbvsdEta_%d_%d_%d", type, sign, real), ";#Delta#eta;<Q_{b}>", bins, dEtaBinsArray, 20000,-10.0-0.0005,10.0-0.0005);

      }
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
