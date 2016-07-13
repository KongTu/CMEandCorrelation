#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260,300};//dummy bins
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double pPb_BinCenter[] = {131.3,162.1,196.7,231.5, 270.4};

const int Nmults = 5;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;


void makeTGraph_pPb_5TeV_Quan(){

	TFile* file[9];

	file[0] = new TFile("../dataPoints/pPb_eff_3part.root");

	TH1D* gr1 = (TH1D*) file[0]->Get("h3PartSS_raw");
	TH1D* gr2 = (TH1D*) file[0]->Get("h3PartOS_raw");

    double value1[22];
    double value1_error[22];
    double value2[22];
    double value2_error[22];
	
	double PbPb_centralityBinCenter_fill[22];

	double x1[22];
	double y1[22];
	double ey1[22];

	double x2[22];
	double y2[22];
	double ey2[22];

	double xbinwidth[22];

    for(int mult = 0; mult < gr1->GetNbinsX(); mult++){

        x1[mult] = gr1->GetBinCenter(mult+1);
        y1[mult] = gr1->GetBinContent(mult+1);
        ey1[mult] = gr1->GetBinError(mult+1);

        x2[mult] = gr2->GetBinCenter(mult+1);
        y2[mult] = gr2->GetBinContent(mult+1);
        ey2[mult] = gr2->GetBinError(mult+1);
    	
    	xbinwidth[mult] = 0.0;

    }

    TGraphErrors* gr5 = new TGraphErrors(gr1->GetNbinsX(), pPb_BinCenter, y1, xbinwidth, ey1);
    TGraphErrors* gr6 = new TGraphErrors(gr1->GetNbinsX(), pPb_BinCenter, y2, xbinwidth, ey2);
    
    TFile t1("../dataPoints/pPb_5TeV_data_Quan.root","RECREATE");

    gr5->Write();
    gr6->Write();



}