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

//eta bins and dEta bins
//double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
double etabins[] = {-2.4,0.0,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;

double ntrkbins[] = {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 150.0, 185.0, 220.0, 260.0, 300.0};
const int Nntrkbins = sizeof(ntrkbins) / sizeof(ntrkbins[0]) - 1;

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

double getSumOnly( TH1D* hist ){

	double integral2 = 0.;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		if(temp != 0 ){
			integral2 = integral2 + temp*center;
		}
	}

	double sum = integral2;
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



	