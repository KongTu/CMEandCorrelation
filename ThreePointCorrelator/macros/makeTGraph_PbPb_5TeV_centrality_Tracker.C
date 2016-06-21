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

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_ntrkCentralityBinCenter[] = {624.8, 810.5, 1025, 1256};

double PbPb_centralityBinCenter[] = {75, 65, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5};

const int Nmults = 1;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;


void makeTGraph_PbPb_5TeV_centrality_Tracker(){

	TFile* file[9];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v12_4.root");

	TH1D* QvsdEta[30][48][3];

	TH1D* delEta3p[30][3];

	TH1D* Ntrk[10];

	for(int mult = 0; mult < Nmults; mult++){

		Ntrk[mult] = (TH1D*) file[mult]->Get("ana/Ntrk");
		cout << Ntrk[mult]->GetMean() << ", ";

		for(int sign = 0; sign < 3; sign++){

			delEta3p[mult][sign] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d",sign));
			
		}
	}

	TH1D* QaQb[30]; TH1D* QaQc[30]; TH1D* QcQb[30];
	TH1D* aveQ3[30][2][2];

	for(int mult = 0; mult < Nmults; mult++){

		QaQb[mult] = (TH1D*)file[mult]->Get("ana/c2_ab");
		QaQc[mult] = (TH1D*)file[mult]->Get("ana/c2_ac");
		QcQb[mult] = (TH1D*)file[mult]->Get("ana/c2_cb");

		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){

				aveQ3[mult][i][j] = (TH1D*)file[mult]->Get(Form("ana/aveQ3_%d_%d",i,j) );
			}
		}
	}

	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
			for(int sign = 0; sign < 3; sign++){
			  
				  QvsdEta[mult][deta][sign] = (TH1D*) file[mult]->Get( Form("ana/QvsdEta_%d_%d",deta,sign) );
				  
			}
		}
	}

	double v2[23];//get corrected v2_3

	for(int mult = 0; mult < Nmults; mult++){
		
		double meanQaQb = QaQb[mult]->GetMean();
		double meanQaQc = QaQc[mult]->GetMean();
		double meanQcQb = QcQb[mult]->GetMean();

		v2[mult] = sqrt( meanQaQc*meanQcQb/meanQaQb );

	}

	TH1D* hist1[3];
	for(int sign = 0; sign < 3; sign++){
		hist1[sign] = new TH1D(Form("hist1_%d_%d",sign), "", NntrkBins, ntrkBins);
	}

	double threeParticleNtrk[9][3];
	double threeParticleNtrkError[9][3];
	double totalWeight[9][3];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){

			for(int deta = 0; deta < 16; deta++){

				double Q_total_real_dEta = QvsdEta[mult][deta][sign]->GetMean();
				double Q_total_real_dEta_error = QvsdEta[mult][deta][sign]->GetMeanError();
				double deltaEtaWeight = delEta3p[mult][sign]->GetBinContent( deta+1 );

				threeParticleNtrk[mult][sign] += Q_total_real_dEta*deltaEtaWeight;
				threeParticleNtrkError[mult][sign] += (Q_total_real_dEta_error*Q_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
				totalWeight[mult][sign] += deltaEtaWeight;

			}	
		}
	}
	
	//pPb:
	for(int sign = 0; sign < 3; sign++){
		for(int mult = 0; mult < Nmults; mult++){

			//pPb(0,7)
			double value = threeParticleNtrk[mult][sign]/totalWeight[mult][sign];
			value = value/v2[mult];
			hist1[sign]->SetBinContent( mult+1, value);
			double error = threeParticleNtrkError[mult][sign]/(totalWeight[mult][sign]*totalWeight[mult][sign]);
			error = sqrt(error)/v2[mult];
			hist1[sign]->SetBinError( mult+1, error);

		}
	}

	TH1D* temp1 = (TH1D*)hist1[0]->Clone("temp1");
	TH1D* temp2 = (TH1D*)hist1[1]->Clone("temp1");

	temp1->Add(temp2, +1);
	temp1->Scale(0.5);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp5 = (TH1D*) hist1[2]->Clone("temp5");
	temp5->SetMarkerStyle(25);
	temp5->SetMarkerColor(kBlue);
	temp5->SetLineColor(kBlue);

    double value1[9];
    double value1_error[9];
    double value2[9];
    double value2_error[9];
	
	double PbPb_centralityBinCenter_fill[9];

    for(int mult = 0; mult < Nmults; mult++){

    	value1[mult] = temp1->GetBinContent(mult+1);
    	value1_error[mult] = temp1->GetBinError(mult+1);

    	value2[mult] = temp5->GetBinContent(mult+1);
    	value2_error[mult] = temp5->GetBinError(mult+1);

    }

    PbPb_centralityBinCenter_fill[0] = 15;
    
    TGraphErrors* gr5 = new TGraphErrors(1, PbPb_centralityBinCenter_fill, value1, xbinwidth, value1_error);
    TGraphErrors* gr6 = new TGraphErrors(1, PbPb_centralityBinCenter_fill, value2, xbinwidth, value2_error);
    
    TFile t1("../dataPoints/PbPb_5TeV_centrality_centrality_data_tracker.root","RECREATE");

    gr5->Write();
    gr6->Write();



}