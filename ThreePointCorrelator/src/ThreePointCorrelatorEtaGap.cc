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
      TH1D* evtWeightedQp3;
      TH1D* evtWeight;
      
      TH2D* QvsdEtaPlusPlus;
      TH2D* QvsdEtaMinusMinus;
      TH2D* QvsdEtaPlusMinus;
      TH2D* QvsdEtaMinusPlus;

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
}


ThreePointCorrelatorEtaGap::~ThreePointCorrelatorEtaGap()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


double getReal(double cos1, double cos2, double cos3, double sin1, double sin2, double sin3){

  double t1 = cos1*cos2*cos3;
  double t2 = cos1*sin2*sin3;
  double t3 = cos2*sin1*sin3;
  double t4 = sin1*sin2*cos3;

  return t1+t2+t3-t4;

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
  if(bestvz < -15.0 || bestvz > 15.0) return;
  
  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

//loop over calo towers (HF)

  double HFqVcos = 0.;
  double HFqVsin = 0.;
  int HFcounts = 0;
  for(unsigned int i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();

        if( fabs(caloEta) < 3 || fabs(caloEta) > 5 ) continue;

        double CosTerm = cos( 2*caloPhi );
        double SinTerm = sin( 2*caloPhi );
        double w = hit.et();

        HFqVcos = HFqVcos + w*CosTerm;
        HFqVsin = HFqVsin + w*SinTerm;
        HFcounts++;
  }

// define eta bins:
  vector<double> etabins;
  double increment = 0.0;
  for(int eta = 0; eta < 49; eta++){

    double initial = -2.4;
    etabins.push_back( initial + increment );
    increment = increment + 0.1;
  }

// initialize Qcos and Qsin
  //q vector in tracker for particle 3, that correlate with HF. 
  double QcosP3 = 0.;
  double QsinP3 = 0.;
  int QcountsP3 = 0;

  double Qcos[48][2];
  double Qsin[48][2];
  int Qcounts[48][2];
  for(int eta = 0; eta < 48; eta++){
    for(int sign = 0; sign < 2; sign++){

      Qcounts[eta][sign] = 0;
      Qcos[eta][sign] = 0.;
      Qsin[eta][sign] = 0.;
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
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        if ( fabs(trk.eta()) > 2.4 || trk.pt() < 0.4  ) continue;
        nTracks++;   

        if( fabs(trk.eta() ) < 1.0 ){

          QcosP3 = QcosP3 + cos( 2*trk.phi() );
          QsinP3 = QsinP3 + sin( 2*trk.phi() );
          QcountsP3++;

        }
        for(unsigned eta = 0; eta < etabins.size()-1; eta++){
          if( trk.eta() > etabins[eta] && trk.eta() < etabins[eta+1] ){

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

  Ntrk->Fill(nTracks);



  double Q = QcosP3*HFqVcos + QsinP3*HFqVsin;
  double weight = QcountsP3*HFcounts;
  double weightedQ = Q/weight;
  double W2 = nTracks*(nTracks-1);
    
  evtWeight->Fill( W2 );
  evtWeightedQp3->Fill( W2*weightedQ );

  for(int ieta = 0; ieta < 48; ieta++){
    for(int jeta = 0; jeta < 48; jeta++){

      if( ieta == jeta ) continue;

      double totalQplusplus = getReal(Qcos[ieta][0],Qcos[jeta][0], HFqVcos, Qsin[ieta][0], Qsin[jeta][0], HFqVsin );
      double totalQminusminus = getReal(Qcos[ieta][1],Qcos[jeta][1], HFqVcos, Qsin[ieta][1], Qsin[jeta][1], HFqVsin );
      double totalQplusminus = getReal(Qcos[ieta][0],Qcos[jeta][1], HFqVcos, Qsin[ieta][0], Qsin[jeta][1], HFqVsin );
      double totalQminusplus = getReal(Qcos[ieta][1],Qcos[jeta][0], HFqVcos, Qsin[ieta][1], Qsin[jeta][0], HFqVsin );

      double Nplusplus = Qcounts[ieta][0]*Qcounts[jeta][0]*HFcounts;
      double Nminusminus = Qcounts[ieta][1]*Qcounts[jeta][1]*HFcounts;
      double Nplusminus = Qcounts[ieta][0]*Qcounts[jeta][1]*HFcounts;
      double Nminusplus = Qcounts[ieta][1]*Qcounts[jeta][0]*HFcounts;

      double deltaEta = fabs(etabins[jeta] - etabins[ieta]);
      if( deltaEta < 0.09999 ){
        cout << "deltaEta: " << deltaEta << endl;
        continue;
      }

      double normalizeQplusplus = totalQplusplus/Nplusplus;
      double normalizeQminusminus = totalQminusminus/Nminusminus;
      double normalizeQplusminus = totalQplusminus/Nplusminus;
      double normalizeQminusplus = totalQminusplus/Nminusplus;
      
      if( Nplusplus == 0 ){
        normalizeQplusplus = 0.0;
      }
      else if( Nminusminus == 0 ){
        normalizeQminusminus = 0.0;
      }
      else if( Nplusminus == 0 ){
        normalizeQplusminus = 0.0;
      }
      else if( Nminusplus == 0 ){
        normalizeQminusplus = 0.0;
      }

      QvsdEtaPlusPlus->Fill(deltaEta, normalizeQplusplus);
      QvsdEtaMinusMinus->Fill(deltaEta, normalizeQminusminus);
      QvsdEtaPlusMinus->Fill(deltaEta, normalizeQplusminus);
      QvsdEtaMinusPlus->Fill(deltaEta, normalizeQminusplus);

    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelatorEtaGap::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",200,0,200);
  evtWeight = fs->make<TH1D>("evtWeight",";evtWeight", 100000,0,100000);
  evtWeightedQp3 = fs->make<TH1D>("evtWeightedQp3",";evtWeightedQp3", 100000,0,500);
  QvsdEtaPlusPlus = fs->make<TH2D>("QvsdEtaPlusPlus",";#Delta#eta;Q_{#phi_{1,+}}Q_{#phi_{2,+}}Q^{*}_{2#phi_{3}}", 48,0,4.8,20000,-0.1,0.1 );
  QvsdEtaMinusMinus = fs->make<TH2D>("QvsdEtaMinusMinus",";#Delta#eta;Q_{#phi_{1,-}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", 48,0,4.8,20000,-0.1,0.1 );
  QvsdEtaPlusMinus = fs->make<TH2D>("QvsdEtaPlusMinus",";#Delta#eta;Q_{#phi_{1,+}}Q_{#phi_{2,-}}Q^{*}_{2#phi_{3}}", 48,0,4.8,20000,-0.1,0.1 );
  QvsdEtaMinusPlus = fs->make<TH2D>("QvsdEtaMinusPlus",";#Delta#eta;Q_{#phi_{1,-}}Q_{#phi_{2,+}}Q^{*}_{2#phi_{3}}", 48,0,4.8,20000,-0.1,0.1 );
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
