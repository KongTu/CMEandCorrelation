#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_ntrkCentralityBinCenter[] = {151.6, 270.2, 441.9, 685.4,1024,1376,1721,40000};
const int Nmults = 8;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;


void makeTGraph_PbPb(){

	TFile* file[8];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_1.root");
	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_2.root");
	file[2] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_3.root");
	file[3] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_4.root");
	file[4] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_5.root");
	file[5] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_6.root");
	file[6] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_7.root");
	file[7] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_8.root");


	TH1D* QvsdEta[30][48][3][2];

	TH1D* delEta3p[30][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				delEta3p[mult][sign][HF] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d_%d",sign,HF));
			}
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
				for(int HF = 0; HF < 2; HF++){
			  
				  QvsdEta[mult][deta][sign][HF] = (TH1D*) file[mult]->Get( Form("ana/QvsdEta_%d_%d_%d",deta,sign,HF) );
				  
				}
			}
		}
	}

	double v2[23][3];//get corrected v2_3

	for(int mult = 0; mult < Nmults; mult++){
		
		double meanQaQb = QaQb[mult]->GetMean();
		double meanQaQc = QaQc[mult]->GetMean();
		double meanQcQb = QcQb[mult]->GetMean();

		double c2_a = meanQaQb*meanQaQc/meanQcQb;
		double c2_b = meanQaQb*meanQcQb/meanQaQc;
		double c2_ab = meanQaQb;

		double bCorr = (aveQ3[mult][0][0]->GetMean() * aveQ3[mult][0][0]->GetMean()) +  ( aveQ3[mult][0][1]->GetMean() * aveQ3[mult][0][1]->GetMean() );
		double aCorr = (aveQ3[mult][1][0]->GetMean() * aveQ3[mult][1][0]->GetMean()) +  ( aveQ3[mult][1][1]->GetMean() * aveQ3[mult][1][1]->GetMean() );
	
		double m1 = (aveQ3[mult][0][0]->GetMean() + aveQ3[mult][1][0]->GetMean())/2.0;
		double m2 = (aveQ3[mult][0][1]->GetMean() + aveQ3[mult][1][1]->GetMean())/2.0;

		double abCorr = m1*m1 + m2*m2;

		v2[mult][0] = sqrt(c2_b - bCorr);
		v2[mult][1] = sqrt(c2_a - aCorr );
		v2[mult][2] = sqrt(c2_ab - abCorr );
	}

	TH1D* hist1[3][2];
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF), "", NntrkBins, ntrkBins);
		}
	}

	double threeParticleNtrk[8][3][2];
	double threeParticleNtrkError[8][3][2];
	double totalWeight[8][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				for(int deta = 0; deta < 15; deta++){

					double Q_total_real_dEta = QvsdEta[mult][deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error = QvsdEta[mult][deta][sign][HF]->GetMeanError();
					double deltaEtaWeight = delEta3p[mult][sign][HF]->GetBinContent( deta+1 );

					threeParticleNtrk[mult][sign][HF] += Q_total_real_dEta*deltaEtaWeight;
					threeParticleNtrkError[mult][sign][HF] += (Q_total_real_dEta_error*Q_total_real_dEta_error)*(deltaEtaWeight*deltaEtaWeight);
					totalWeight[mult][sign][HF] += deltaEtaWeight;

				}
			}
		}
	}
	
	//pPb:
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			for(int mult = 0; mult < 8; mult++){

				//pPb(0,7)
				double value = threeParticleNtrk[mult][sign][HF]/totalWeight[mult][sign][HF];
				value = value/v2[mult][HF];
				hist1[sign][HF]->SetBinContent( mult+1, value);
				double error = threeParticleNtrkError[mult][sign][HF]/(totalWeight[mult][sign][HF]*totalWeight[mult][sign][HF]);
				hist1[sign][HF]->SetBinError( mult+1, sqrt(error));

			}
		}
	}

	TH1D* temp1 = (TH1D*)hist1[0][0]->Clone("temp1");
	TH1D* temp2 = (TH1D*)hist1[0][1]->Clone("temp1");
	TH1D* temp3 = (TH1D*)hist1[1][0]->Clone("temp1");
	TH1D* temp4 = (TH1D*)hist1[1][1]->Clone("temp1");
	
	temp1->Add(temp2, +1);
	temp1->Add(temp3, +1);
	temp1->Add(temp4, +1);

	temp1->Scale(0.25);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp5 = (TH1D*) hist1[2][0]->Clone("temp5");
	TH1D* temp6 = (TH1D*) hist1[2][1]->Clone("temp6");
	
	temp5->Add(temp6, +1);
	
	temp5->Scale(0.5);
	temp5->SetMarkerStyle(25);
	temp5->SetMarkerColor(kBlue);
	temp5->SetLineColor(kBlue);


    double value1[8];
    double value1_error[8];
    double value2[8];
    double value2_error[8];


    for(int mult = 0; mult < 8; mult++){

    	value1[mult] = temp1->GetBinContent(mult+1);
    	value1_error[mult] = temp1->GetBinError(mult+1);

    	value2[mult] = temp5->GetBinContent(mult+1);
    	value2_error[mult] = temp5->GetBinError(mult+1);

    }

	TGraphErrors* gr5 = new TGraphErrors(8, PbPb_ntrkBinCenter, value1, xbinwidth, value1_error);
    TGraphErrors* gr6 = new TGraphErrors(8, PbPb_ntrkBinCenter, value2, xbinwidth, value2_error);

    TFile t1("../dataPoints/PbPb_data.root","RECREATE");
    gr5->Write();
    gr6->Write();



}