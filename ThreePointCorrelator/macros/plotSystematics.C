#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
//rebin option:
double dEtaReBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,4.8};
const int NdEtaReBins = sizeof(dEtaReBins) / sizeof(dEtaReBins[0]) - 1;
//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

//short range:
double dEtaBinsShortRange[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
const int NdEtaBinsShortRange = sizeof(dEtaBinsShortRange) / sizeof(dEtaBinsShortRange[0]) - 1;

double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
const int Nmults = 4;

bool doShortRange_ = false;

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

void plotSystematics(){

	gStyle->SetErrorX(0);
	TGaxis::SetMaxDigits(3);

	string Label1; 
	string Label2;
	string Label3;
	string Label4; 

	Label1 = "pPb default";
	Label2 = "pPb |Vz| < 3 ";
	Label3 = "PbPb default";
	Label4 = "PbPb |Vz| < 3 ";


	TFile* file[16];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");
	file[1] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_Systematics_v1.root");
	file[2] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v2.root");
	file[3] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_Systematics_v1.root");

	file[4] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_Systematics_v5.root");

	TH1D* QvsdEta[16][48][3][2];

	TH1D* delEta3p[16][3][2];

	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				delEta3p[mult][sign][HF] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d_%d",sign,HF));
			}
		}
	}

	TH1D* QaQb[16]; TH1D* QaQc[16]; TH1D* QcQb[16];
	TH1D* aveQ3[16][2][2];



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

	double v2[16][3];//get corrected v2_3

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
	
//tracker 	
	for(int sign = 0; sign < 3; sign++){
		delEta3p[4][sign][0] = (TH1D*) file[4]->Get(Form("ana/delEta3p_%d",sign));
		delEta3p[4][sign][1] = (TH1D*) file[4]->Get(Form("ana/delEta3p_%d",sign));

	}

	QaQb[4] = (TH1D*)file[4]->Get("ana/c2_ab");
	QaQc[4] = (TH1D*)file[4]->Get("ana/c2_ac");
	QcQb[4] = (TH1D*)file[4]->Get("ana/c2_cb");

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){

			QvsdEta[4][deta][sign][0] = (TH1D*) file[4]->Get( Form("ana/QvsdEta_%d_%d",deta,sign) );
			QvsdEta[4][deta][sign][1] = (TH1D*) file[4]->Get( Form("ana/QvsdEta_%d_%d",deta,sign) );

		}
	}

	double meanQaQb = QaQb[4]->GetMean();
	double meanQaQc = QaQc[4]->GetMean();
	double meanQcQb = QcQb[4]->GetMean();

	double c2_c = meanQcQb*meanQaQc/meanQaQb;
	v2[4][0] = sqrt(c2_c);
	v2[4][1] = sqrt(c2_c);

//


	if( !doShortRange_ ){

		TH1D* hist1[8][3][2];
		TH1D* hist2[8][3][2];
		for(int mult = 0; mult < 5; mult++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
					hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaReBins2, dEtaReBins2);
					hist2[mult][sign][HF] = new TH1D(Form("hist2_%d_%d_%d",mult,sign,HF),"test", NdEtaReBins2, dEtaReBins2);
				}
			}
		}

		for(int mult = 0; mult < 5; mult++){
		for(int deta = 0; deta < NdEtaReBins2; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){

					if(deta < 8){

						double Q_total_real_dEta1 = QvsdEta[mult][3*deta][sign][HF]->GetMean();
						double Q_total_real_dEta_error1 = QvsdEta[mult][3*deta][sign][HF]->GetMeanError();

						double Q_total_real_dEta2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMean();
						double Q_total_real_dEta_error2 = QvsdEta[mult][3*deta+1][sign][HF]->GetMeanError();

						double Q_total_real_dEta3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMean();
						double Q_total_real_dEta_error3 = QvsdEta[mult][3*deta+2][sign][HF]->GetMeanError();

						double weight1 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+1 );
						double weight2 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+2 );
						double weight3 = delEta3p[mult][sign][HF]->GetBinContent( 3*deta+3 );
						
						double value = weightedAverage(weight1, weight2, weight3, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3);
						double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );
						
						hist1[mult][sign][HF]->SetBinContent(deta+1, value );
						hist1[mult][sign][HF]->SetBinError(deta+1, error );


					}

					// if(deta >= 8){

					// hist1[mult][sign][HF]->SetBinContent(deta+1, 0.0 );
					// hist1[mult][sign][HF]->SetBinError(deta+1, 0.0 );

					// }
					else{

						double Q_total_real_dEta1 = QvsdEta[mult][27][sign][HF]->GetMean();
						double Q_total_real_dEta_error1 = QvsdEta[mult][27][sign][HF]->GetMeanError();
						
						double Q_total_real_dEta2 = QvsdEta[mult][28][sign][HF]->GetMean();
						double Q_total_real_dEta_error2 = QvsdEta[mult][28][sign][HF]->GetMeanError();

						double weight1 = delEta3p[mult][sign][HF]->GetBinContent( 28 );
						double weight2 = delEta3p[mult][sign][HF]->GetBinContent( 29 );

						double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
						double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
						
						hist1[mult][sign][HF]->SetBinContent(deta+1, value );
						hist1[mult][sign][HF]->SetBinError(deta+1,  error);
					}


				}
			}
		}
	}
	
	}
	else{
	 
		TH1D* hist1[8][3][2];
		TH1D* hist2[8][3][2];
		for(int mult = 0; mult < Nmults; mult++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){
					hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaBinsShortRange, dEtaBinsShortRange);
					hist2[mult][sign][HF] = new TH1D(Form("hist2_%d_%d_%d",mult,sign,HF),"test", NdEtaBinsShortRange, dEtaBinsShortRange);
				}
			}
		}

		for(int mult = 0; mult < Nmults; mult++){
			for(int deta = 0; deta < NdEtaBinsShortRange; deta++){
				for(int sign = 0; sign < 3; sign++){
					for(int HF = 0; HF < 2; HF++){

					double Q_total_real_dEta = QvsdEta[mult][deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error = QvsdEta[mult][deta][sign][HF]->GetMeanError();
					
					hist1[mult][sign][HF]->SetBinContent(deta+1, Q_total_real_dEta );
					hist1[mult][sign][HF]->SetBinError(deta+1,  Q_total_real_dEta_error);
						
					}
				}
			}
		}
	}

	TH1D* base1 = makeHist("base1","same-sign, Pb-going","#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8);
    base1->GetXaxis()->SetTitleColor(kBlack);
    base1->GetYaxis()->SetRangeUser(-0.002, 0.002);
    base1->GetYaxis()->SetTitleOffset(1.0);

//tracker    

 //    TH1D* tracker1 = (TH1D*) hist1[4][0][0]->Clone("tracker1");
 //    TH1D* tracker2 = (TH1D*) hist1[4][1][0]->Clone("tracker2");
 //    tracker1->Add(tracker2, +1);
 //    tracker1->Scale(0.5);
 //    tracker1->Scale(1.0/v2[4][0]);
 //    tracker1->SetMarkerStyle(20);
 //    tracker1->SetMarkerColor(kBlack);
	// tracker1->SetLineColor(kBlack);

 //    TH1D* tracker3 = (TH1D*) hist1[4][2][0]->Clone("tracker3");
 //    tracker3->Scale(1.0/v2[4][0]);
 //    tracker3->SetMarkerColor(kGreen+3);
	// tracker3->SetLineColor(kGreen+3);
	// tracker3->SetMarkerStyle(21);

//default:
	//like sign Pb-going
	TH1D* temp1 = (TH1D*) hist1[0][0][0]->Clone("temp1");
	TH1D* temp2 = (TH1D*) hist1[0][1][0]->Clone("temp2");

	temp1->Add(temp2, +1);
	temp1->Scale(0.5);
	temp1->Scale(1.0/v2[0][0]);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);
	temp1->SetMarkerStyle(20);

	//like sign p-going 
	TH1D* temp3 = (TH1D*) hist1[0][0][1]->Clone("temp3");
	TH1D* temp4 = (TH1D*) hist1[0][1][1]->Clone("temp4");

	temp3->Add(temp4, +1);
	temp3->Scale(0.5);
	temp3->Scale(1.0/v2[0][1]);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);
	temp3->SetMarkerStyle(20);

	//unlike sign Pb-going
	TH1D* temp5 = (TH1D*) hist1[0][2][0]->Clone("temp5");
	temp5->Scale(1.0/v2[0][0]);
	temp5->SetMarkerColor(kBlue);
	temp5->SetLineColor(kBlue);
	temp5->SetMarkerStyle(21);

	//unlike sign p-going
	TH1D* temp6 = (TH1D*) hist1[0][2][1]->Clone("temp6");
	temp6->Scale(1.0/v2[0][1]);
	temp6->SetMarkerColor(kBlue);
	temp6->SetLineColor(kBlue);
	temp6->SetMarkerStyle(21);

	//like sign PbPb
	TH1D* temp7 = (TH1D*) hist1[2][0][0]->Clone("temp7");
	TH1D* temp8 = (TH1D*) hist1[2][0][1]->Clone("temp8");
	TH1D* temp9 = (TH1D*) hist1[2][1][0]->Clone("temp9");
	TH1D* temp10 = (TH1D*) hist1[2][1][1]->Clone("temp10");
	temp7->Add(temp8, +1);
	temp7->Add(temp9, +1);
	temp7->Add(temp10, +1);
	temp7->Scale(0.25);
	temp7->Scale(1.0/v2[2][2]);
	temp7->SetMarkerColor(kRed);
	temp7->SetLineColor(kRed);
	temp7->SetMarkerStyle(20);

	//unlike sign PbPb
	TH1D* temp10_1 = (TH1D*) hist1[2][2][0]->Clone("temp10_1");
	TH1D* temp10_2 = (TH1D*) hist1[2][2][1]->Clone("temp10_2");
	temp10_1->Add(temp10_2, +1);
	temp10_1->Scale(0.5);
	temp10_1->Scale(1.0/v2[2][2]);
	temp10_1->SetMarkerColor(kBlue);
	temp10_1->SetLineColor(kBlue);
	temp10_1->SetMarkerStyle(21);

//systematics:

	//like sign Pb-going
	TH1D* temp11 = (TH1D*) hist1[1][0][0]->Clone("temp11");
	TH1D* temp22 = (TH1D*) hist1[1][1][0]->Clone("temp22");

	temp11->Add(temp22, +1);
	temp11->Scale(0.5);
	temp11->Scale(1.0/v2[1][0]);
	temp11->SetMarkerColor(kBlack);
	temp11->SetLineColor(kBlack);
	temp11->SetMarkerStyle(20);

	//like sign p-going 
	TH1D* temp33 = (TH1D*) hist1[1][0][1]->Clone("temp33");
	TH1D* temp44 = (TH1D*) hist1[1][1][1]->Clone("temp44");

	temp33->Add(temp44, +1);
	temp33->Scale(0.5);
	temp33->Scale(1.0/v2[1][1]);
	temp33->SetMarkerColor(kBlack);
	temp33->SetLineColor(kBlack);
	temp33->SetMarkerStyle(20);

	//unlike sign Pb-going
	TH1D* temp55 = (TH1D*) hist1[1][2][0]->Clone("temp55");
	temp55->Scale(1.0/v2[1][0]);
	temp55->SetMarkerColor(kGreen+3);
	temp55->SetLineColor(kGreen+3);
	temp55->SetMarkerStyle(21);

	//unlike sign p-going
	TH1D* temp66 = (TH1D*) hist1[1][2][1]->Clone("temp66");
	temp66->Scale(1.0/v2[1][1]);
	temp66->SetMarkerColor(kGreen+3);
	temp66->SetLineColor(kGreen+3);
	temp66->SetMarkerStyle(21);

	//like sign PbPb
	TH1D* temp77 = (TH1D*) hist1[3][0][0]->Clone("temp77");
	TH1D* temp88 = (TH1D*) hist1[3][0][1]->Clone("temp88");
	TH1D* temp99 = (TH1D*) hist1[3][1][0]->Clone("temp99");
	TH1D* temp100 = (TH1D*) hist1[3][1][1]->Clone("temp100");
	temp77->Add(temp88, +1);
	temp77->Add(temp99, +1);
	temp77->Add(temp100, +1);
	temp77->Scale(0.25);
	temp77->Scale(1.0/v2[3][2]);
	temp77->SetMarkerColor(kBlack);
	temp77->SetLineColor(kBlack);
	temp77->SetMarkerStyle(20);

	//unlike sign PbPb
	TH1D* temp10_11 = (TH1D*) hist1[3][2][0]->Clone("temp10_11");
	TH1D* temp10_22 = (TH1D*) hist1[3][2][1]->Clone("temp10_22");
	temp10_11->Add(temp10_22, +1);
	temp10_11->Scale(0.5);
	temp10_11->Scale(1.0/v2[3][2]);
	temp10_11->SetMarkerColor(kGreen+3);
	temp10_11->SetLineColor(kGreen+3);
	temp10_11->SetMarkerStyle(21);


//Canvas:
	TCanvas* c1 = makeMultiCanvas("c1","c1",3,2);
	c1->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	base1->Draw();
	temp1->Draw("Psame");
	temp11->Draw("Psame");

	TLegend *w1 = new TLegend(0.35,0.2,0.5,0.4);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(temp1, Label1.c_str());
    w1->AddEntry(temp11, Label2.c_str());
    w1->Draw("same");

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base2 = (TH1D*)base1->Clone("base2");
	base2->SetTitle("same-sign, p-going");
	base2->Draw();
	temp3->Draw("Psame");
	temp33->Draw("Psame");

	TLegend *w2 = new TLegend(0.35,0.2,0.5,0.4);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp3, Label1.c_str());
    w2->AddEntry(temp33, Label2.c_str());
    w2->Draw("same");

    c1->cd(3);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base_PbPb = (TH1D*)base1->Clone("base_PbPb");
	base_PbPb->SetTitle("same-sign, PbPb");
	base_PbPb->Draw();
	temp7->Draw("Psame");
	temp77->Draw("Psame");


	TLatex* PbPb_sign = new TLatex(0.2, 0.6, "No HM trigger in PbPb");
    PbPb_sign->SetNDC();
    PbPb_sign->SetTextSize(23);
    PbPb_sign->SetTextFont(43);
    PbPb_sign->SetTextColor(kBlack);
    //PbPb_sign->Draw("same");

	TLegend *w3 = new TLegend(0.35,0.2,0.5,0.4);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(temp7, Label3.c_str());
    w3->AddEntry(temp77, Label4.c_str());
    w3->Draw("same");

	c1->cd(4);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base3 = (TH1D*)base1->Clone("base3");
	base3->SetTitle("opposite-sign, Pb-going");
	base3->Draw();
	temp5->Draw("Psame");
	temp55->Draw("Psame");

	TLegend *w4 = new TLegend(0.35,0.2,0.5,0.4);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(20);
    w4->SetTextFont(43);
    w4->AddEntry(temp5, Label1.c_str());
    w4->AddEntry(temp55, Label2.c_str());
    w4->Draw("same");

	c1->cd(5);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base4 = (TH1D*)base1->Clone("base4");
	base4->SetTitle("opposite-sign, p-going");
	base4->Draw();
	temp6->Draw("Psame");
	temp66->Draw("Psame");

	TLegend *w5 = new TLegend(0.35,0.2,0.5,0.4);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(20);
    w5->SetTextFont(43);
    w5->AddEntry(temp6, Label1.c_str());
    w5->AddEntry(temp66, Label2.c_str());
    w5->Draw("same");

	c1->cd(6);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base_PbPb_unlike = (TH1D*)base1->Clone("base_PbPb_unlike");
	base_PbPb_unlike->SetTitle("opposite-sign, PbPb");
	base_PbPb_unlike->Draw();
	temp10_1->Draw("Psame");
	temp10_11->Draw("Psame");
    //PbPb_sign->Draw("same");

	TLegend *w6 = new TLegend(0.35,0.2,0.5,0.4);
    w6->SetLineColor(kWhite);
    w6->SetFillColor(0);
    w6->SetTextSize(20);
    w6->SetTextFont(43);
    w6->AddEntry(temp10_1, Label3.c_str());
    w6->AddEntry(temp10_11, Label4.c_str());
    w6->Draw("same");
	
	TH1D* base5 = makeHist("base5","same-sign, Pb-going","#Delta#eta", "default - systematic checks", 48,0,4.8);
    base5->GetXaxis()->SetTitleColor(kBlack);
    base5->GetYaxis()->SetRangeUser(-0.0005, 0.0005);
    base5->GetYaxis()->SetTitleOffset(1.5);

    TH1D* ratio1 = (TH1D*)temp1->Clone("ratio1");
    TH1D* ratio11 = (TH1D*)temp11->Clone("ratio11");

    TH1D* ratio3 = (TH1D*)temp3->Clone("ratio3");
    TH1D* ratio33 = (TH1D*)temp33->Clone("ratio33");

    TH1D* ratio5 = (TH1D*)temp5->Clone("ratio5");
    TH1D* ratio55 = (TH1D*)temp55->Clone("ratio55");

    TH1D* ratio6 = (TH1D*)temp6->Clone("ratio6");
    TH1D* ratio66 = (TH1D*)temp66->Clone("ratio66");

    TH1D* ratio7 = (TH1D*)temp7->Clone("ratio7");
    TH1D* ratio77 = (TH1D*)temp77->Clone("ratio77");

    TH1D* ratio8 = (TH1D*)temp10_1->Clone("ratio8");
    TH1D* ratio88 = (TH1D*)temp10_11->Clone("ratio88");

   	TCanvas* c2 = makeMultiCanvas("c2","c2",3,2);
	c2->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	base5->Draw();
	ratio1->Add( ratio11, -1 );
	ratio1->Draw("Psame");

	c2->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base6 = (TH1D*)base5->Clone("base6");
	base6->SetTitle("same-sign, p-going");
	base6->Draw();
	ratio3->Add( ratio33, -1 );
	ratio3->Draw("Psame");

	c2->cd(3);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base6_PbPb = (TH1D*)base5->Clone("base6_PbPb");
	base6_PbPb->SetTitle("same-sign, PbPb");
	base6_PbPb->Draw();
	ratio7->Add( ratio77, -1 );
	ratio7->Draw("Psame");
    //PbPb_sign->Draw("same");

	c2->cd(4);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base7 = (TH1D*)base5->Clone("base7");
	base7->SetTitle("opposite-sign, Pb-going");
	base7->Draw();
	ratio5->Add( ratio55, -1 );	
	ratio5->Draw("Psame");

	c2->cd(5);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base8 = (TH1D*)base5->Clone("base8");
	base8->SetTitle("opposite-sign, p-going");
	base8->Draw();
	ratio6->Add( ratio66, -1 );
	ratio6->Draw("Psame");

	c2->cd(6);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.14);
	TH1D* base8_PbPb = (TH1D*)base5->Clone("base8_PbPb");
	base8_PbPb->SetTitle("opposite-sign, PbPb");
	base8_PbPb->Draw();
	ratio8->Add( ratio88, -1 );
	ratio8->Draw("Psame");
    //PbPb_sign->Draw("same");
    
 //    TCanvas* c3 = makeCanvas("c3","c3");
 //    gPad->SetTicks();
	// gPad->SetLeftMargin(0.15);
	// gPad->SetBottomMargin(0.14);
	// TH1D* base9_tracker = (TH1D*)base1->Clone("base9_tracker");
	// base9_tracker->SetTitle("PbPb tracker");
	// base9_tracker->GetXaxis()->SetRangeUser(0,3.0);
	// base9_tracker->Draw();

	// tracker1->Draw("Psame");
	// tracker3->Draw("Psame");
	// temp7->Draw("Psame");
	// temp10_1->Draw("Psame");
	
	// TLegend *w7 = new TLegend(0.55,0.2,0.7,0.4);
 //    w7->SetLineColor(kWhite);
 //    w7->SetFillColor(0);
 //    w7->SetTextSize(20);
 //    w7->SetTextFont(43);
 //    w7->AddEntry(tracker1, "same, tracker");
 //    w7->AddEntry(temp7, "same, HF");
 //    w7->AddEntry(tracker3, "opposite, tracker");
 //    w7->AddEntry(temp10_1, "opposite, HF");
 //    w7->Draw("same");

	// TLatex* r3 = new TLatex(0.2, 0.84, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
 //    r3->SetNDC();
 //    r3->SetTextSize(23);
 //    r3->SetTextFont(43);
 //    r3->SetTextColor(kBlack);
 //    r3->Draw("same");

 //    TLatex* lmult = new TLatex(0.2, 0.76, "185 #leq N^{offline}_{trk} < 220");
 //    lmult->SetNDC();
 //    lmult->SetTextSize(23);
 //    lmult->SetTextFont(43);
 //    lmult->SetTextColor(kBlack);
 //    lmult->Draw("same");

	c1->Print("../systematics/systematics_v1.pdf");
	c2->Print("../systematics/systematics_diff_v1.pdf");



}