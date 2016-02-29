#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

void function(){

	cout << "loading function!" << endl;
}

//define eta bins:
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

//define all the 2D and 1D histograms:

TH2D* corrVsEtaPlusPlus = new TH2D("corrVsEtaPlusPlus", "corrVsEtaPlusPlus", Nbins, 0, 4.8, 100000, -0.00001, 0.00001);
TH2D* corrVsEtaMinusMinus = new TH2D("corrVsEtaMinusMinus", "corrVsEtaMinusMinus", Nbins, 0, 4.8, 100000, -0.00001, 0.00001);
TH2D* corrVsEtaPlusMinus = new TH2D("corrVsEtaPlusMinus", "corrVsEtaPlusMinus", Nbins, 0, 4.8, 100000, -0.00001, 0.00001);

TH1D* Qplusplus1D[48];
TH1D* Qminusminus1D[48];
TH1D* Qplusminus1D[48];
TH1D* Qminusplus1D[48];

TH1D* corrVsEtaPlusPlus1D[48];
TH1D* corrVsEtaMinusMinus1D[48];
TH1D* corrVsEtaPlusMinus1D[48];

TH1D* hist = makeHist("hist","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", 48,0,4.8, kBlack);
TH1D* hist1 = makeHist("hist1","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kBlack);
TH1D* hist2 = makeHist("hist2","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kRed);
TH1D* hist3 = makeHist("hist3","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kBlue);
TH1D* hist4 = makeHist("hist4","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kGreen-3);

TH1D* corrHistPlusPlus = makeHist("corrHistPlusPlus","","#Delta#eta", "correction factor", 48,0,4.8, kBlack);
TH1D* corrHistMinusMinus = makeHist("corrHistMinusMinus","","#Delta#eta", "correction factor", 48,0,4.8, kRed);
TH1D* corrHistPlusMinus = makeHist("corrHistPlusMinus","","#Delta#eta", "correction factor", 48,0,4.8, kGreen-3);


//****** member functions **********//

double getV2( TH1D* hist ){

	double integral2 = 0;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		if(temp != 0 ){
			integral2 = integral2 + temp*center;
		}
	}

	double v2 = sqrt( ( integral2/hist->GetEntries() ) );
	return v2;

}

double getNormalizedSum( TH1D* hist ){

	double integral2 = 0.;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		if(temp != 0 ){
			integral2 = integral2 + temp*center;
		}
	}

	double sum = integral2/( hist->GetEntries() );
	return sum;

}

double getWeightedSum( TH1D* hist1, TH1D* hist2 ){
	
	double integral1 = 0.;
	for(int i = 0; i < hist1->GetNbinsX(); i++ ){

		double temp = hist1->GetBinContent(i+1);
		double center = hist1->GetBinCenter(i+1);
		if(temp != 0 ){
			integral1 = integral1 + temp*center;
		}
	}

	double integral2 = 0.;
	for(int i = 0; i < hist2->GetNbinsX(); i++ ){

		double temp = hist2->GetBinContent(i+1);
		double center = hist2->GetBinCenter(i+1);
		if(temp != 0 ){
			integral2 = integral2 + temp*center;
		}
	}

	double sum = integral1/integral2;
	return sum;

}

double getReal(double cos1, double cos2, double cos3, double sin1, double sin2, double sin3){

  double t1 = cos1*cos2*cos3;
  double t2 = cos1*sin2*sin3;
  double t3 = cos2*sin1*sin3;
  double t4 = sin1*sin2*cos3;

  return t1+t2+t3-t4;

}

	