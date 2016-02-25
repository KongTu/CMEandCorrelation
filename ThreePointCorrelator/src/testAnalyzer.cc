// -*- C++ -*-
//
// Package:    testAnalyzer
// Class:      testAnalyzer
// 
/**\class testAnalyzer testAnalyzer.cc CMEandCorrelation/testAnalyzer/src/testAnalyzer.cc

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

class testAnalyzer : public edm::EDAnalyzer {
   public:
      explicit testAnalyzer(const edm::ParameterSet&);
      ~testAnalyzer();

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
      TH1D* Qp3;
      TH1D* evtWeightedQp3;
      TH1D* evtWeight;
      TH1D* averageCos;
      TH1D* averageSin;
      
      //culmulants:
      TH2D* QvsdEtaPlusPlus;
      TH2D* QvsdEtaMinusMinus;
      TH2D* QvsdEtaPlusMinus;
      TH2D* QvsdEtaMinusPlus;

      //two and single particle sum
      //HF:
      TH1D* HFcosSum;
      TH1D* HFsinSum;
      TH1D* weightSum; 

      TH1D* TRKcosPlusSum[48];
      TH1D* TRKsinPlusSum[48];
      TH1D* TRKcosMinusSum[48];
      TH1D* TRKsinMinusSum[48];

      TH2D* EvsEta;
      TH2D* ETvsEta;

      TH1D* testDeta;
      TH1D* testV2;

      int Nmin_;
      int Nmax_;

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
testAnalyzer::testAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  towerSrc_ = iConfig.getParameter<edm::InputTag>("towerSrc");
  
  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

}


testAnalyzer::~testAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

    
// ------------ method called for each event  ------------
void
testAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z();

  
  //first selection; vertices
  if(bestvz < -15.0 || bestvz > 15.0) return;
  
  Handle<CaloTowerCollection> towers;
  iEvent.getByLabel(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

//loop over calo towers (HF)

  double HFqVcosPlus = 0.;
  double HFqVsinPlus = 0.;
  int HFminusCounts = 0;
  double ETTplus = 0.;

  double HFqVcosMinus = 0.;
  double HFqVsinMinus = 0.;
  int HFplusCounts = 0;
  double ETTminus = 0.;

  double tV2 = 0.0;
  int norm = 0;

  for(unsigned i = 0; i < towers->size(); ++i){

        const CaloTower & hit= (*towers)[i];

        double caloEta = hit.eta();
        double caloPhi = hit.phi();
        double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );
        //double energy = hit.emEnergy() + hit.hadEnergy();

        for(unsigned j = 0; j < towers->size(); j++){

            const CaloTower & hit2 = (*towers)[j];

            if( i == j ) continue;
            if( (hit.eta() > 4.4 && hit.eta() < 5.0 && hit2.eta() < -4.4 && hit2.eta() > -5) || (hit2.eta() > 4.4 && hit2.eta() < 5.0 && hit.eta() < -4.4 && hit.eta() > -5) ){

              tV2 += cos( 2*( hit.phi() - hit2.phi() ));
              norm++;
            }
            
        }

        if( caloEta < -4.4 && caloEta > -5 ){

          HFqVcosMinus = HFqVcosMinus + w*cos( 2*caloPhi );
          HFqVsinMinus = HFqVsinMinus + w*sin( 2*caloPhi );
         
          HFminusCounts++; 
          ETTminus += w;

        }
        if( caloEta < 5 && caloEta > 4.4 ){
          
          HFqVcosPlus = HFqVcosPlus + w*cos( 2*caloPhi );
          HFqVsinPlus = HFqVsinPlus + w*sin( 2*caloPhi );
          
          HFplusCounts++;
          ETTplus += w;
        }


  }

  testV2->Fill( tV2/norm );

//weight by ET() and renormalize by ETT in order to have dimensionless 

  if( HFplusCounts == 0 || HFminusCounts == 0 ) return;

  HFqVcosPlus = HFqVcosPlus/ETTplus;
  HFqVsinPlus = HFqVsinPlus/ETTplus;

  HFqVcosMinus = HFqVcosMinus/ETTminus;
  HFqVsinMinus = HFqVsinMinus/ETTminus;

  double Q = HFqVcosMinus*HFqVcosPlus + HFqVsinMinus*HFqVsinPlus;
  double W2 = ETTminus*ETTplus;
  evtWeight->Fill( W2 );
  evtWeightedQp3->Fill( W2*Q );
  Qp3->Fill( Q );

}


// ------------ method called once each job just before starting event loop  ------------
void 
testAnalyzer::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  evtWeight = fs->make<TH1D>("evtWeight",";evtWeight", 10000000,0,500);
  evtWeightedQp3 = fs->make<TH1D>("evtWeightedQp3",";evtWeightedQp3", 1000000,0,50);
  Qp3 = fs->make<TH1D>("Qp3",";Qp3", 1000000,0,50);

  testV2 = fs->make<TH1D>("testV2",";cos", 100000,-1,1);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
testAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
testAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
testAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
testAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
testAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
testAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(testAnalyzer);
