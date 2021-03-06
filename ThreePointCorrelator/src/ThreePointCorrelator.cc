// -*- C++ -*-
//
// Package:    ThreePointCorrelator
// Class:      ThreePointCorrelator
// 
/**\class ThreePointCorrelator ThreePointCorrelator.cc CMEandCorrelation/ThreePointCorrelator/src/ThreePointCorrelator.cc

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

class ThreePointCorrelator : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelator(const edm::ParameterSet&);
      ~ThreePointCorrelator();

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
      std::string vertexSrc_;

      TH1D* Ntrk;
      TH1D* NtrkPlus;
      TH1D* NtrkMinus;
      
      TH2D* QvsNtrkPlusPlusPlus;
      TH2D* QvsNtrkMinusMinusMinus;
      TH2D* QvsNtrkPlusPlusMinus;
      TH2D* QvsNtrkMinusMinusPlus;

      TH2D* QvsNtrkPlusMinusPlus;
      TH2D* QvsNtrkPlusMinusMinus;
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
ThreePointCorrelator::ThreePointCorrelator(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");

}


ThreePointCorrelator::~ThreePointCorrelator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

double getQ3(double a1, double a2, double a3){

  double temp1 = cos(a1)*cos(a2)*cos(2*a3);
  double temp2 = cos(a1)*sin(a2)*sin(2*a3);
  double temp3 = sin(a1)*sin(2*a3)*cos(a2);
  double temp4 = sin(a1)*sin(a2)*cos(2*a3);

  return temp1 + temp2 + temp3 - temp4;
}

std::vector<double> combination;
std::vector< std::vector<double>> plusPlusCombination;
std::vector< std::vector<double>> minusMinusCombination;
std::vector< std::vector<double>> allCombination;

void go(unsigned offset, int k, std::vector<double> angle, std::string sign){

  if (k == 0) {
    if(sign == "plusplus"){
      plusPlusCombination.push_back(combination);
      return;
    }
    else if( sign == "minusminus"){
      minusMinusCombination.push_back(combination);
      return;
    }
    else if( sign == "all"){
      allCombination.push_back(combination);
      return;
    }
    else{

      std::cout << "wrong option! Can be either plusplus, minusminus, and all" << std::endl;
      return;
    }

  }
  for (unsigned i = offset; i <= angle.size() - k; ++i) {
    combination.push_back(angle[i]);
    go(i+1, k-1, angle, sign);
    combination.pop_back();
  }

} 

int choose(int n, int k){

    if (k == 0) return 1;
    return (n * choose(n - 1, k - 1)) / k;
}
    
// ------------ method called for each event  ------------
void
ThreePointCorrelator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

  vector<double> angle;
  vector<double> anglePlusPlus;
  vector<double> angleMinusMinus;

  int nTracks = 0;
  int nTracksPlus = 0;
  int nTracksMinus = 0;

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
        if( trk.charge() == 1 ) nTracksPlus++;
        if( trk.charge() == -1 ) nTracksMinus++;
      
        if( fabs( trk.eta() ) > 1.0 || trk.pt() < 0.4 ) continue;
        
        angle.push_back( trk.phi() );
        if( trk.charge() == 1 ) anglePlusPlus.push_back( trk.phi() );
        else if( trk.charge() == -1 ) angleMinusMinus.push_back( trk.phi() );
        
  } 
  //store all pairs combination in allCombination
  go(0,2,anglePlusPlus,"plusplus");
  combination.clear();
  go(0,2,angleMinusMinus,"minusminus");
  combination.clear();
  go(0,2,angle,"all");

  double total3QPlusPlusPlus = 0.;
  double total3QPlusPlusMinus = 0.;
  double total3QMinusMinusPlus = 0.;
  double total3QMinusMinusMinus = 0.;

  double total3QPlusMinusPlus = 0.;
  double total3QPlusMinusMinus = 0.;
  
  //start from plusplus pair
  for(unsigned i = 0; i < plusPlusCombination.size(); i++){
    for(unsigned j = 0; j < anglePlusPlus.size(); j++){

      if( anglePlusPlus[j] == plusPlusCombination[i][0] || anglePlusPlus[j] == plusPlusCombination[i][1] ) continue;
      total3QPlusPlusPlus = total3QPlusPlusPlus + getQ3(plusPlusCombination[i][0], plusPlusCombination[i][1], anglePlusPlus[j]);
    }

    for(unsigned j = 0; j < angleMinusMinus.size(); j++){

      //if( angleMinusMinus[j] == plusPlusCombination[i][0] || angleMinusMinus[j] == plusPlusCombination[i][1] ) continue;
      total3QPlusPlusMinus = total3QPlusPlusMinus + getQ3(plusPlusCombination[i][0], plusPlusCombination[i][1], angleMinusMinus[j]);
    }
  }

  double N3plus = (anglePlusPlus.size()-2)*choose(anglePlusPlus.size(), 2);
  double averageQ3plus = total3QPlusPlusPlus/N3plus;

  double N2plus1minus = (angleMinusMinus.size())*choose(anglePlusPlus.size(), 2);
  double averageQ2plus1minus = total3QPlusPlusMinus/N2plus1minus;

  QvsNtrkPlusPlusPlus->Fill(nTracks, averageQ3plus);
  QvsNtrkPlusPlusMinus->Fill(nTracks, averageQ2plus1minus);

  //start with minus minus pair
  for(unsigned i = 0; i < minusMinusCombination.size(); i++){
    for(unsigned j = 0; j < angleMinusMinus.size(); j++){

      if( angleMinusMinus[j] == minusMinusCombination[i][0] || angleMinusMinus[j] == minusMinusCombination[i][1] ) continue;
      total3QMinusMinusMinus = total3QMinusMinusMinus + getQ3(minusMinusCombination[i][0], minusMinusCombination[i][1], angleMinusMinus[j]);
    }

    for(unsigned j = 0; j < anglePlusPlus.size(); j++){

      //if( anglePlusPlus[j] == minusMinusCombination[i][0] || anglePlusPlus[j] == minusMinusCombination[i][1] ) continue;
      total3QMinusMinusPlus = total3QMinusMinusPlus + getQ3(minusMinusCombination[i][0], minusMinusCombination[i][1], anglePlusPlus[j]);
    }
  }

  double N3minus = (angleMinusMinus.size()-2)*choose(angleMinusMinus.size(), 2);
  double averageQ3minus = total3QMinusMinusMinus/N3minus;

  double N2minus1plus = (anglePlusPlus.size())*choose(angleMinusMinus.size(), 2);
  double averageQ2minus1plus = total3QMinusMinusPlus/N2minus1plus;

  QvsNtrkMinusMinusMinus->Fill(nTracks, averageQ3minus);
  QvsNtrkMinusMinusPlus->Fill(nTracks, averageQ2minus1plus);


  for(unsigned i = 0; i < anglePlusPlus.size(); i++){

    for(unsigned j = 0; j < angleMinusMinus.size(); j++){

      for( unsigned k = 0; k < anglePlusPlus.size(); k++){

        if( i == k ) continue;
        total3QPlusMinusPlus = total3QPlusMinusPlus + getQ3(anglePlusPlus[i], angleMinusMinus[j], anglePlusPlus[k]);
      }
      for( unsigned k = 0; k < angleMinusMinus.size(); k++){

        if( j == k ) continue;
        total3QPlusMinusMinus = total3QPlusMinusMinus + getQ3(anglePlusPlus[i], angleMinusMinus[j], angleMinusMinus[k]);
      }
    }
  }

  double Nplus = anglePlusPlus.size() * angleMinusMinus.size() * (anglePlusPlus.size() - 1);
  double Nminus = anglePlusPlus.size() * angleMinusMinus.size() * (angleMinusMinus.size() - 1);
 
  double averageQplus = total3QPlusMinusPlus/Nplus;
  double averageQminus = total3QPlusMinusMinus/Nminus;

  QvsNtrkPlusMinusPlus->Fill(nTracks, averageQplus);
  QvsNtrkPlusMinusMinus->Fill(nTracks, averageQminus);

  Ntrk->Fill(nTracks);
  NtrkPlus->Fill(nTracksPlus);
  NtrkMinus->Fill(nTracksMinus);
}


// ------------ method called once each job just before starting event loop  ------------
void 
ThreePointCorrelator::beginJob()
{

  edm::Service<TFileService> fs;
    
  TH3D::SetDefaultSumw2();

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",200,0,200);
  NtrkPlus = fs->make<TH1D>("NtrkPlus",";NtrkPlus",200,0,200);
  NtrkMinus = fs->make<TH1D>("NtrkMinus",";NtrkMinus",200,0,200);
  QvsNtrkPlusPlusPlus = fs->make<TH2D>("QvsNtrkPlusPlusPlus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);
  QvsNtrkPlusPlusMinus = fs->make<TH2D>("QvsNtrkPlusPlusMinus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);
  QvsNtrkMinusMinusPlus = fs->make<TH2D>("QvsNtrkMinusMinusPlus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);
  QvsNtrkMinusMinusMinus = fs->make<TH2D>("QvsNtrkMinusMinusMinus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);
  QvsNtrkPlusMinusPlus = fs->make<TH2D>("QvsNtrkPlusMinusPlus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);
  QvsNtrkPlusMinusMinus = fs->make<TH2D>("QvsNtrkPlusMinusMinus", ";Ntrk;<cos(#phi_{1} + #phi_{2} - 2#phi_{3})>", 300,0,300, 20000,-0.1,0.1);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ThreePointCorrelator::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ThreePointCorrelator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ThreePointCorrelator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ThreePointCorrelator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ThreePointCorrelator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ThreePointCorrelator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ThreePointCorrelator);
