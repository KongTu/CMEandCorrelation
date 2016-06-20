#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260,300};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5, 270.4};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

//double PbPb_ntrkCentralityBinCenter[] = {151.6, 270.2, 441.9, 685.4,1024,1376,1721};
double PbPb_ntrkCentralityBinCenter[] = {404.1, 717.6, 1141, 1782, 2662.4, 3577.6, 4474};
double PbPb_centralityBinCenter[] = {65,55,45,35,25,15,7.5,2.5};

const int Nmults = 9;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;

void makeTGraph_ALICE_STAR(){


//ALICE:

    //VERZO:
	// double alice_centrality_oppo1[] = {1.08E-4, -2.6E-5, -6.4E-5, -2.5E-5, -1.4E-5, -1.0E-6, 0.0 };
	// double alice_centrality_oppo1_error[] = {8.0E-5, 4.1E-5, 2.3E-5, 1.4E-5, 9.0E-6, 6.0E-6, 4.0E-6};

	// double alice_centrality_same1[] = {-3.66E-4, -2.48E-4, -1.47E-4, -9.8E-5, -4.6E-5, -1.0E-5, 0.0};
	// double alice_centrality_same1_error[] = {8.1E-5, 4.1E-5, 2.4E-5, 1.4E-5, 9.0E-6, 6.0E-6, 4.0E-6};
	//don't use VERZO for now. 


	double alice_centrality_oppo1[9];
	double alice_centrality_oppo1_error[9];
	
	double alice_centrality_same1[9];
	double alice_centrality_same1_error[9];

//ZDC
	double p8370_d4x1y1_yval[] = { -2.4E-5, -1.6E-5, -5.0E-6, -2.3E-5, -7.8E-5, 6.0E-6, -4.65E-4, -1.63E-4 };
	double p8370_d4x1y1_yerrminus[] = { 2.5E-5, 1.0E-5, 1.2E-5, 1.8E-5, 3.4E-5, 7.6E-5, 2.25E-4, 8.95E-4 };

	double p8370_d4x1y2_yval[] = { -3.6E-5, -5.6E-5, -9.4E-5, -1.71E-4, -3.2E-4, -3.79E-4, -6.61E-4, -0.001213 };
	double p8370_d4x1y2_yerrminus[] = { 2.4E-5, 1.0E-5, 1.1E-5, 1.7E-5, 3.3E-5, 7.5E-5, 2.24E-4, 8.96E-4 };

//VEZERO
	double p8370_d3x1y1_yval[] = { 0.0, -1.0E-6, -1.4E-5, -2.5E-5, -6.4E-5, -2.6E-5, 1.08E-4, 3.15E-4, 9.89E-4 };
  	double p8370_d3x1y1_yerrminus[] = { 4.0E-6, 6.0E-6, 9.0E-6, 1.4E-5, 2.3E-5, 4.1E-5, 8.0E-5, 1.22E-4, 2.17E-4 };

	double p8370_d3x1y2_yval[] = { 0.0, -1.0E-5, -4.6E-5, -9.8E-5, -1.47E-4, -2.48E-4, -3.66E-4, -3.7E-4, -3.93E-4 };
  	double p8370_d3x1y2_yerrminus[] = { 4.0E-6, 6.0E-6, 9.0E-6, 1.4E-5, 2.4E-5, 4.1E-5, 8.1E-5, 1.77E-4, 3.77E-4 };

	double alice_centrality_oppo2[9];
	double alice_centrality_oppo2_error[9];
	
	double alice_centrality_same2[9];
	double alice_centrality_same2_error[9];

	//TPC cumulant
	double p8370_d2x1y1_yval[] = { -6.0E-6, 2.0E-6, -5.0E-6, -1.1E-5, -1.0E-5, -3.0E-5, -1.0E-5, 2.47E-4 };
	double p8370_d2x1y1_yerrminus[] = { 3.0E-6, 4.0E-6, 3.0E-6, 5.0E-6, 7.0E-6, 1.2E-5, 2.2E-5, 4.6E-5 };

	double p8370_d2x1y2_yval[] = { -1.0E-5, -2.3E-5, -5.4E-5, -9.7E-5, -1.6E-4, -2.69E-4, -3.72E-4, -4.45E-4 };
	double p8370_d2x1y2_yerrminus[] = { 3.0E-6, 4.0E-6, 3.0E-6, 5.0E-6, 7.0E-6, 1.2E-5, 2.2E-5, 4.6E-5};


	for(int mult = 0; mult < 8; mult++){

		alice_centrality_oppo1[mult] = p8370_d4x1y1_yval[7-mult];
		alice_centrality_oppo1_error[mult] = p8370_d4x1y1_yerrminus[7-mult];

		alice_centrality_same1[mult] = p8370_d4x1y2_yval[7-mult];
		alice_centrality_same1_error[mult] = p8370_d4x1y2_yerrminus[7-mult];

		alice_centrality_oppo2[mult] = p8370_d2x1y1_yval[7-mult];
		alice_centrality_oppo2_error[mult] = p8370_d2x1y1_yerrminus[7-mult];

		alice_centrality_same2[mult] = p8370_d2x1y2_yval[7-mult];
		alice_centrality_same2_error[mult] = p8370_d2x1y2_yerrminus[7-mult];

		
	}

//STAR 200GeV AuAu, CuCu

	double star_centrality_dump_same1[] = {-2.72475e-05, -4.84767e-05, -8.43622e-05, -0.000139391, -0.000212998, -0.000310464, -0.000449019, -0.000531625};
	double star_centrality_dump_same1_error[] = {2.51676e-06, 2.38477e-06, 1.96375e-06, 2.50148e-06, 3.87569e-06, 6.65112e-06, 1.36251e-05, 3.49563e-05};

	double star_centrality_dump_oppo1[] = {-1.16758e-05, -8.2939e-06, -7.9278e-06, -5.81744e-06, -3.32537e-06, 1.71515e-05, 5.88244e-05, 0.000205355};
	double star_centrality_dump_oppo1_error[] = {2.55046e-06, 2.41575e-06, 1.99013e-06, 2.52889e-06, 3.90957e-06, 6.65273e-06, 1.3518e-05, 3.40235e-05};

	double star_centrality_dump_same2[] = {-8.23628e-05, -0.000120783, -0.000170984, -0.000268721, -0.000360827, -0.000493512, -0.000666972, 100};
	double star_centrality_dump_same2_error[] = {7.62553e-06, 8.77766e-06, 7.82902e-06, 1.08416e-05, 1.69498e-05, 2.87805e-05, 5.33616e-05, 0.0};

	double star_centrality_dump_oppo2[] = {1.92707e-05, 2.76513e-05, 4.28017e-05, 9.16021e-05, 0.000133615, 0.00025283, 0.000513905, 100};
	double star_centrality_dump_oppo2_error[] = {7.68691e-06, 8.82822e-06, 7.84248e-06, 1.08158e-05, 1.676e-05, 2.82131e-05, 5.13577e-05, 0};


	double star_centrality_same1[9];
	double star_centrality_same1_error[9];
	double star_centrality_oppo1[9];
	double star_centrality_oppo1_error[9];

	double star_centrality_same2[9];
	double star_centrality_same2_error[9];
	double star_centrality_oppo2[9];
	double star_centrality_oppo2_error[9];

	double PbPb_centralityBinCenter_fill[9];

	for(int mult = 0; mult < 8; mult++){

		star_centrality_same1[mult] = star_centrality_dump_same1[7-mult];
		star_centrality_same1_error[mult] = star_centrality_dump_same1_error[7-mult];

		star_centrality_oppo1[mult] = star_centrality_dump_oppo1[7-mult];
		star_centrality_oppo1_error[mult] = star_centrality_dump_oppo1_error[7-mult];

		star_centrality_same2[mult] = star_centrality_dump_same2[7-mult];
		star_centrality_same2_error[mult] = star_centrality_dump_same2_error[7-mult];

		star_centrality_oppo2[mult] = star_centrality_dump_oppo2[7-mult];
		star_centrality_oppo2_error[mult] = star_centrality_dump_oppo2_error[7-mult];

		PbPb_centralityBinCenter_fill[mult] = 80-PbPb_centralityBinCenter[mult];

	}


	//use the flipped value to place on the axis, so that it can be later drawed reversely. 

	TGraphErrors* gr1 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, alice_centrality_same1, xbinwidth, alice_centrality_same1_error);
	TGraphErrors* gr2 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, alice_centrality_oppo1, xbinwidth, alice_centrality_oppo1_error);
	
	TGraphErrors* gr3 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, star_centrality_same1, xbinwidth, star_centrality_same1_error);
	TGraphErrors* gr4 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, star_centrality_oppo1, xbinwidth, star_centrality_oppo1_error);

	TGraphErrors* gr5 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, star_centrality_same2, xbinwidth, star_centrality_same2_error);
	TGraphErrors* gr6 = new TGraphErrors(8, PbPb_centralityBinCenter_fill, star_centrality_oppo2, xbinwidth, star_centrality_oppo2_error);

	TFile t1("../dataPoints/ALICE_STAR_backup_data.root","RECREATE");
    gr1->Write();
    gr2->Write();
    gr3->Write();
    gr4->Write();
    gr5->Write();
    gr6->Write();




}