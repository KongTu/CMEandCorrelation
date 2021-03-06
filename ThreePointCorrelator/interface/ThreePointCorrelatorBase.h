#ifndef ThreePointCorrelatorBase_
#define ThreePointCorrelatorBase_


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


using namespace std;
using namespace reco;
using namespace edm;


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
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);
      virtual double fRand(double fMin, double fMax);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* vtxZ;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* cbinHist;
      TH1D* q2_mag;
      TH1D* delEta3p[3][2];
      TH1D* delEta2p[3];

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* aveQ3[2][2];//calculate the correction on v2

//end v2

      TH1D* QvsdEta[48][3][2];
      TH1D* PvsdEta[48][3];

      TH1D* XY_real[48][3][2];TH1D* XY_imag[48][3][2];
      TH1D* XZ_real[48][3][2];TH1D* XZ_imag[48][3][2];
      TH1D* YZ_real[48][3][2];TH1D* YZ_imag[48][3][2];
      TH1D* X_real[48][3][2]; TH1D* X_imag[48][3][2];
      TH1D* Y_real[48][3][2]; TH1D* Y_imag[48][3][2];
      TH1D* Z_real[48][3][2]; TH1D* Z_imag[48][3][2];

      int Nmin_;
      int Nmax_;

      double etaTracker_;
      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;
      bool do3pTracker_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorEtaGapQ2 : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorEtaGapQ2(const edm::ParameterSet&);
      ~ThreePointCorrelatorEtaGapQ2();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);
      virtual double fRand(double fMin, double fMax);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* Ntrk_Q2;
      TH1D* vtxZ;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* cbinHist;
      TH1D* q2_mag;
      TH2D* q2_tracker_HF;
      TH1D* delEta3p[3][2];
      TH1D* delEta2p[3];

      TH1D* QnCQnC[2];
      TH1D* QnCQnC_ave[2];

      TH1D* QnQnA[2];
      TH1D* QnAQnB[2];
      TH1D* QnAQnC[2];
      TH1D* QnBQnC[2];

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* aveQ3[2][2];//calculate the correction on v2

//end v2

      TH2D* QvsV2[3][2];

      TH1D* QvsdEta[48][3][2];
      TH1D* PvsdEta[48][3];

      TH1D* XY_real[48][3][2];TH1D* XY_imag[48][3][2];
      TH1D* XZ_real[48][3][2];TH1D* XZ_imag[48][3][2];
      TH1D* YZ_real[48][3][2];TH1D* YZ_imag[48][3][2];
      TH1D* X_real[48][3][2]; TH1D* X_imag[48][3][2];
      TH1D* Y_real[48][3][2]; TH1D* Y_imag[48][3][2];
      TH1D* Z_real[48][3][2]; TH1D* Z_imag[48][3][2];

      int Nmin_;
      int Nmax_;

      double etaTracker_;
      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;
      double holeLeft_;
      double holeRight_;
      double holesize_;
      double q2max_;
      double q2min_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;
      bool do3pTracker_;
      bool doTrackerQ2_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorEtaGapTracker : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorEtaGapTracker(const edm::ParameterSet&);
      ~ThreePointCorrelatorEtaGapTracker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* vtxZ;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* cbinHist;
      TH1D* delEta3p[3];
      TH1D* delEta2p[3];

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* aveQ3[2][2];//calculate the correction on v2

//end v2

      TH1D* QvsdEta[48][3];
      TH1D* PvsdEta[48][3];

      TH1D* XY_real[48][3][2];TH1D* XY_imag[48][3][2];
      TH1D* XZ_real[48][3][2];TH1D* XZ_imag[48][3][2];
      TH1D* YZ_real[48][3][2];TH1D* YZ_imag[48][3][2];
      TH1D* X_real[48][3][2]; TH1D* X_imag[48][3][2];
      TH1D* Y_real[48][3][2]; TH1D* Y_imag[48][3][2];
      TH1D* Z_real[48][3][2]; TH1D* Z_imag[48][3][2];

      int Nmin_;
      int Nmax_;

      double etaTracker_;
      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;
      double holeLeft_;
      double holeRight_;
      double holesize_;
      double etaGap_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;
      bool do3pTracker_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorEtaTest : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorEtaTest(const edm::ParameterSet&);
      ~ThreePointCorrelatorEtaTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N2, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;

      TH1D* plusCount[1000];
      TH1D* minusCount[1000];
//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* aveQ3[2][2];//calculate the correction on v2

//end v2

      TH1D* QvsdEta[1000][3][2];

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;

};

class ThreePointCorrelatorEtaGapNestedLoop : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorEtaGapNestedLoop(const edm::ParameterSet&);
      ~ThreePointCorrelatorEtaGapNestedLoop();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;
      TH1D* trkPt;
      TH1D* trk_eta;

      //v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

      TH1D* QvsdEta[48][3][2];


      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;
      bool do3pTracker_;


      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorGen : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorGen(const edm::ParameterSet&);
      ~ThreePointCorrelatorGen();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;
      edm::InputTag genParticleSrc_;

      //correction table
      TH2D* effTable;
      TH2D* Ntrk2D;

      TH1D* Ntrk;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;
      TH1D* trkPt;
      TH1D* trk_eta;

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

//end v2
      TH1D* delEta3p[3][2];
      TH1D* delEta2p[3];

      TH1D* QvsdEta[48][3][2];
      TH1D* PvsdEta[48][3];

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorGenTracker : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorGenTracker(const edm::ParameterSet&);
      ~ThreePointCorrelatorGenTracker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;
      edm::InputTag genParticleSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* delEta3p[3];
      TH1D* delEta2p[3];

//v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;

//end v2

      TH1D* QvsdEta[48][3];
      TH1D* PvsdEta[48][3];

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

class ThreePointCorrelatorNestedLoop : public edm::EDAnalyzer {
   public:
      explicit ThreePointCorrelatorNestedLoop(const edm::ParameterSet&);
      ~ThreePointCorrelatorNestedLoop();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual double get3Real(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3);
      virtual double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3);
      virtual double get2Real( double R1, double R2, double I1, double I2);
      virtual double get2RealOverlap( double R1, double R2, double I1, double I2);
      virtual double get2Imag( double R1, double R2, double I1, double I2);
      virtual double get2ImagOverlap( double R1, double R2, double I1, double I2);

      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag towerSrc_;
      std::string vertexSrc_;
      edm::InputTag genParticleSrc_;

      //correction table
      TH2D* effTable;

      TH1D* Ntrk;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* cbinHist;
      TH1D* trkPt;
      TH1D* trk_eta;

      //v2
      TH1D* c2_ab;
      TH1D* c2_ac;
      TH1D* c2_cb;
      
      TH1D* QvsdEta[48][3][2];

      int Nmin_;
      int Nmax_;

      double etaLowHF_;
      double etaHighHF_;
      double vzLow_;
      double vzHigh_;
      double ptLow_;
      double ptHigh_;
      double offlineptErr_;
      double offlineDCA_;
      double holeLeft_;
      double holeRight_;
      double holesize_;

      bool useCentrality_;
      bool useBothSide_;
      bool reverseBeam_;
      bool messAcceptance_;
      bool doEffCorrection_;

      std::vector<double> etaBins_;
      std::vector<double> dEtaBins_;
      std::vector<double> ptBins_;

};

#endif