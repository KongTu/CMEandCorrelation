#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
//rebin option1:
double dEtaReBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,4.8};
const int NdEtaReBins = sizeof(dEtaReBins) / sizeof(dEtaReBins[0]) - 1;

double dEtaReBinCenter[] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.6,3.1,3.8,4.5};

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
const int Nmults = 3;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;

double weightedAverage(double a1, double a2, double a3, double eta1, double eta2, double eta3){

	double temp1 = a1*eta1 + a2*eta2 + a3*eta3;
	double temp2 = (a1+a2+a3);

	return temp1/temp2;
}

double weightedAverageError(double a1, double a2, double a3, double etaError1, double etaError2, double etaError3){

	double temp1 = (a1/(a1+a2+a3))*(a1/(a1+a2+a3));
	double temp2 = etaError1*etaError1;
	double temp3 = (a2/(a1+a2+a3))*(a2/(a1+a2+a3));
	double temp4 = etaError2*etaError2;
	double temp5 = (a3/(a1+a2+a3))*(a3/(a1+a2+a3));
	double temp6 = etaError3*etaError3;

	double total = temp1*temp2 + temp3*temp4 + temp5*temp6;

	return sqrt(total);

}

void plotDeltaEtaResultsTracker(){

	gStyle->SetErrorX(0);

	TFile* file[2];

	//file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");
	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_Systematics_v9.root");

	TH1D* QvsdEta[48][3];
	TH1D* delEta3p[3];

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){

			delEta3p[sign] = (TH1D*) file[0]->Get(Form("ana/delEta3p_%d",sign));
		}
	}
	
	TH1D* QaQb;
	TH1D* QaQc; 
	TH1D* QcQb;
	TH1D* aveQ3[2][2];

	QaQb = (TH1D*)file[0]->Get("ana/c2_ab");
	QaQc = (TH1D*)file[0]->Get("ana/c2_ac");
	QcQb = (TH1D*)file[0]->Get("ana/c2_cb");

	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){

			aveQ3[i][j] = (TH1D*)file[0]->Get(Form("ana/aveQ3_%d_%d",i,j) );
		}
	}
	
	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){
		  
		  QvsdEta[deta][sign] = (TH1D*) file[0]->Get( Form("ana/QvsdEta_%d_%d",deta,sign) );
			  
		}
	}
	
	double v2[4];//get corrected v2_3

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;
	double c2_c = meanQaQc*meanQcQb/meanQaQb;

	double bCorr = (aveQ3[0][0]->GetMean() * aveQ3[0][0]->GetMean()) +  ( aveQ3[0][1]->GetMean() * aveQ3[0][1]->GetMean() );
	double aCorr = (aveQ3[1][0]->GetMean() * aveQ3[1][0]->GetMean()) +  ( aveQ3[1][1]->GetMean() * aveQ3[1][1]->GetMean() );

	double m1 = (aveQ3[0][0]->GetMean() + aveQ3[1][0]->GetMean())/2.0;
	double m2 = (aveQ3[0][1]->GetMean() + aveQ3[1][1]->GetMean())/2.0;

	double abCorr = m1*m1 + m2*m2;

	v2[0] = sqrt(c2_b - bCorr);
	v2[1] = sqrt(c2_a - aCorr );
	v2[2] = sqrt(c2_ab - abCorr );
	v2[3] = sqrt(c2_c);

	cout << "v2[0]: " << v2[0] << endl;
	cout << "v2[1]: " << v2[1] << endl;
	cout << "v2[2]: " << v2[2] << endl;
	cout << "v2[3]: " << v2[3] << endl;

	TH1D* hist1[3];
	for(int sign = 0; sign < 3; sign++){
		hist1[sign] = new TH1D(Form("hist1_%d",sign),"test", NdEtaReBins2, dEtaReBins2);
	}	
	
	for(int deta = 0; deta < NdEtaReBins2; deta++){
		for(int sign = 0; sign < 3; sign++){

			if(deta < 8){

				double Q_total_real_dEta1 = QvsdEta[3*deta][sign]->GetMean();
				double Q_total_real_dEta_error1 = QvsdEta[3*deta][sign]->GetMeanError();

				double Q_total_real_dEta2 = QvsdEta[3*deta+1][sign]->GetMean();
				double Q_total_real_dEta_error2 = QvsdEta[3*deta+1][sign]->GetMeanError();

				double Q_total_real_dEta3 = QvsdEta[3*deta+2][sign]->GetMean();
				double Q_total_real_dEta_error3 = QvsdEta[3*deta+2][sign]->GetMeanError();

				double weight1 = delEta3p[sign]->GetBinContent( 3*deta+1 );
				double weight2 = delEta3p[sign]->GetBinContent( 3*deta+2 );
				double weight3 = delEta3p[sign]->GetBinContent( 3*deta+3 );
				
				double value = weightedAverage(weight1, weight2, weight3, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3);
				double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );				

				hist1[sign]->SetBinContent(deta+1, value );
				hist1[sign]->SetBinError(deta+1, error );

			}
			if(deta == 8){

				hist1[sign]->SetBinContent(deta+1, 10000.0 );
				hist1[sign]->SetBinError(deta+1, 0.0 );

			}
			if(deta == 9){

				double Q_total_real_dEta1 = QvsdEta[27][sign]->GetMean();
				double Q_total_real_dEta_error1 = QvsdEta[27][sign]->GetMeanError();
				
				double Q_total_real_dEta2 = QvsdEta[28][sign]->GetMean();
				double Q_total_real_dEta_error2 = QvsdEta[28][sign]->GetMeanError();

				double weight1 = delEta3p[sign]->GetBinContent( 28 );
				double weight2 = delEta3p[sign]->GetBinContent( 29 );

				double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
				double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
				cout << "value: " << value << endl;
				cout << "error: " << error << endl;
				hist1[sign]->SetBinContent(deta+1, value );
				hist1[sign]->SetBinError(deta+1,  error);
			}
		}
	}

	TH1D* base1 = makeHist("base1", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012, 0.0013);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	TH1D* temp1 = (TH1D*)hist1[0]->Clone("temp1");
	temp1->Add(hist1[1], +1);
	temp1->Scale(0.5);
	temp1->Scale(1.0/v2[3]);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[2]->Clone("temp2");
	temp2->SetMarkerStyle(25);
	temp2->Scale(1.0/v2[3]);
	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);

	TCanvas* c1 = makeCanvas("c1", "c1");
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

	temp1->Draw("Psame");
	temp2->Draw("Psame");

	TCanvas* c2 = makeCanvas("c2","c2");
	temp2->Add(temp1, -1);
	temp2->GetYaxis()->SetRangeUser(-0.0007,0.0018);

	temp2->Draw();

















}