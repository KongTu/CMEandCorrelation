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

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

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

void plotDeltaEtaResultChecks(){

	gStyle->SetErrorX(0);

	TFile* file[16];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_v46.root");//no eff
	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v5_1.root");//default
	file[2] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v7.root");//right eff

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

	TH1D* hist1[8][3][2];
	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
				hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaReBins2, dEtaReBins2);
			}
		}
	}

	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaReBins2; deta++){
			for(int sign = 0; sign < 3; sign++){
				for(int HF = 0; HF < 2; HF++){

					if(deta < 9){

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

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012,0.001);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);


	TH1D* temp1 = (TH1D*)hist1[0][0][0]->Clone("temp1");
	TH1D* temp2 = (TH1D*) hist1[0][0][1]->Clone("temp2");
	TH1D* temp3 = (TH1D*)hist1[0][1][0]->Clone("temp3");
	TH1D* temp4 = (TH1D*) hist1[0][1][1]->Clone("temp4");

	temp1->Add( temp2, +1);
	temp1->Add( temp3, +1);
	temp1->Add( temp4, +1);

	temp1->Scale(0.25);
	temp1->Scale(1.0/v2[0][2]);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp11 = (TH1D*) hist1[0][2][0]->Clone("temp11");
	temp11->Add( hist1[0][2][1], +1);
	temp11->Scale(0.5);
	temp11->Scale(1.0/v2[0][2]);
	temp11->SetMarkerStyle(25);
	temp11->SetMarkerSize(1.4);
	temp11->SetMarkerColor(kBlue);
	temp11->SetLineColor(kBlue);

	TH1D* temp5 = (TH1D*)hist1[1][0][0]->Clone("temp5");
	temp5->Add(hist1[1][0][1], +1);
	temp5->Add(hist1[1][1][0], +1);
	temp5->Add(hist1[1][1][1], +1);

	temp5->Scale(0.25);
	temp5->Scale(1.0/v2[1][2]);
	temp5->SetMarkerStyle(20);
	temp5->SetMarkerSize(1.4);
	temp5->SetMarkerColor(kRed);
	temp5->SetLineColor(kRed);

	TH1D* temp9 = (TH1D*) hist1[1][2][0]->Clone("temp9");
	TH1D* temp10 = (TH1D*) hist1[1][2][1]->Clone("temp10");
	temp9->Add(temp10, +1);
	temp9->Scale(0.5);
	temp9->Scale(1.0/v2[1][2]);
	temp9->SetMarkerStyle(21);
	temp9->SetMarkerColor(kBlue);
	temp9->SetMarkerSize(1.4);
	temp9->SetLineColor(kBlue);

	TH1D* right1 = (TH1D*)hist1[2][0][0]->Clone("right1");
	TH1D* right2 = (TH1D*) hist1[2][0][1]->Clone("right2");
	TH1D* right3 = (TH1D*)hist1[2][1][0]->Clone("right3");
	TH1D* right4 = (TH1D*) hist1[2][1][1]->Clone("right4");

	right1->Add( right2, +1);
	right1->Add( right3, +1);
	right1->Add( right4, +1);

	right1->Scale(0.25);
	right1->Scale(1.0/v2[2][2]);
	right1->SetMarkerStyle(20);
	right1->SetMarkerSize(1.4);
	right1->SetMarkerColor(kRed);
	right1->SetLineColor(kRed);

	TH1D* right11 = (TH1D*) hist1[2][2][0]->Clone("right11");
	right11->Add( hist1[2][2][1], +1);
	right11->Scale(0.5);
	right11->Scale(1.0/v2[2][2]);
	right11->SetMarkerStyle(21);
	right11->SetMarkerSize(1.4);
	right11->SetMarkerColor(kBlue);
	right11->SetLineColor(kBlue);

    TLatex* r41 = new TLatex(0.24, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.05, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.55,0.95, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.052);
    
    TLatex* r44 = new TLatex(0.69,0.95, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(24);
    r44->SetTextFont(53);
    
    TLatex* r45 = new TLatex(0.05, 0.80, "30-35%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r46 = new TLatex(0.24, 0.80, "30-35%");
    r46->SetNDC();
    r46->SetTextSize(23);
    r46->SetTextFont(43);
    r46->SetTextColor(kBlack);

    TLegend *w4 = new TLegend(0.3,0.2,0.95,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->SetNColumns(2);
    w4->AddEntry(temp1, "  ");
    w4->AddEntry(temp11, "  No tracking correction");
    
    w4->AddEntry(temp5, "  ");
    w4->AddEntry(temp9, "  default correction");

	TLatex* latex1 = new TLatex(0.28, 0.37, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.39, 0.37, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    
    TCanvas* c4 = new TCanvas("c4","c4",1000,600);
	c4->Divide(2,1,0,0);
	c4->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

    double value1[50];
    double value1_error[50];
    double value2[50];
    double value2_error[50];
    double value3[50];
    double value3_error[50];
    double value4[50];
    double value4_error[50];
    double value5[50];
    double value5_error[50];
    double value6[50];
    double value6_error[50];

    TBox* box1[50];
    TBox* box2[50];
    TBox* box3[50];
    TBox* box4[50];
    TBox* box5[50];
    TBox* box6[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value1[deta] = temp1->GetBinContent(deta+1);
    	value1_error[deta] = temp1->GetBinError(deta+1);

    	value2[deta] = temp11->GetBinContent(deta+1);
    	value2_error[deta] = temp11->GetBinError(deta+1);

    	value3[deta] = temp5->GetBinContent(deta+1);
    	value3_error[deta] = temp5->GetBinError(deta+1);

    	value4[deta] = temp9->GetBinContent(deta+1);
    	value4_error[deta] = temp9->GetBinError(deta+1);

    	value5[deta] = right1->GetBinContent(deta+1);
    	value5_error[deta] = right1->GetBinError(deta+1);

    	value6[deta] = right11->GetBinContent(deta+1);
    	value6_error[deta] = right11->GetBinError(deta+1);

    }

	for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.065;
    	double ye = total_systematics_pPb;

    	box1[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value1[deta]-ye,dEtaReBinCenter2[deta]+xe,value1[deta]+ye);
		box1[deta]->SetFillColor(kRed);
        box1[deta]->SetFillColorAlpha(kGray+1,0.4);
        box1[deta]->SetFillStyle(1001);
    	box1[deta]->SetLineWidth(0);
    	box1[deta]->SetLineColor(kRed);
        box1[deta]->Draw("SAME");

		box2[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value2[deta]-ye,dEtaReBinCenter2[deta]+xe,value2[deta]+ye);
		box2[deta]->SetFillColor(kBlue);
        box2[deta]->SetFillColorAlpha(kGray+1,0.4);
        box2[deta]->SetFillStyle(1001);
    	box2[deta]->SetLineWidth(0);
    	box2[deta]->SetLineColor(kBlue);
        box2[deta]->Draw("SAME");

        box3[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value3[deta]-ye,dEtaReBinCenter2[deta]+xe,value3[deta]+ye);
		box3[deta]->SetFillColor(kRed);
        box3[deta]->SetFillColorAlpha(kGray+1,0.4);
        box3[deta]->SetFillStyle(1001);
    	box3[deta]->SetLineWidth(0);
    	box3[deta]->SetLineColor(kRed);
        box3[deta]->Draw("SAME");

		box4[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value4[deta]-ye,dEtaReBinCenter2[deta]+xe,value4[deta]+ye);
		box4[deta]->SetFillColor(kBlue);
        box4[deta]->SetFillColorAlpha(kGray+1,0.4);
        box4[deta]->SetFillStyle(1001);
    	box4[deta]->SetLineWidth(0);
    	box4[deta]->SetLineColor(kBlue);
        box4[deta]->Draw("SAME");
    }

    temp1->Draw("Psame");
	temp11->Draw("Psame");
	temp5->Draw("Psame");
	temp9->Draw("Psame");
	r41->Draw("same");
    w4->Draw("same");
    latex1->Draw("same");
    latex2->Draw("same");
    r46->Draw("same");

    c4->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base1->Draw();
	r45->Draw("same");

    double value5[50];
    double value5_error[50];
    double value6[50];
    double value6_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value5[deta] = temp5->GetBinContent(deta+1);
    	value5_error[deta] = temp5->GetBinError(deta+1);

    	value6[deta] = temp9->GetBinContent(deta+1);
    	value6_error[deta] = temp9->GetBinError(deta+1);

    }	

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.065;
    	double ye = total_systematics_PbPb;

    	box1[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value1[deta]-ye,dEtaReBinCenter2[deta]+xe,value1[deta]+ye);
		box1[deta]->SetFillColor(kRed);
        box1[deta]->SetFillColorAlpha(kGray+1,0.4);
        box1[deta]->SetFillStyle(1001);
    	box1[deta]->SetLineWidth(0);
    	box1[deta]->SetLineColor(kRed);
        box1[deta]->Draw("SAME");

		box2[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value2[deta]-ye,dEtaReBinCenter2[deta]+xe,value2[deta]+ye);
		box2[deta]->SetFillColor(kBlue);
        box2[deta]->SetFillColorAlpha(kGray+1,0.4);
        box2[deta]->SetFillStyle(1001);
    	box2[deta]->SetLineWidth(0);
    	box2[deta]->SetLineColor(kBlue);
        box2[deta]->Draw("SAME");

    	box5[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value5[deta]-ye,dEtaReBinCenter2[deta]+xe,value5[deta]+ye);
		box5[deta]->SetFillColor(kRed);
        box5[deta]->SetFillColorAlpha(kGray+1,0.4);
        box5[deta]->SetFillStyle(1001);
    	box5[deta]->SetLineWidth(0);
    	box5[deta]->SetLineColor(kRed);
        box5[deta]->Draw("SAME");

		box6[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value6[deta]-ye,dEtaReBinCenter2[deta]+xe,value6[deta]+ye);
		box6[deta]->SetFillColor(kBlue);
        box6[deta]->SetFillColorAlpha(kGray+1,0.4);
        box6[deta]->SetFillStyle(1001);
    	box6[deta]->SetLineWidth(0);
    	box6[deta]->SetLineColor(kBlue);
        box6[deta]->Draw("SAME");
    }

    TLegend *w40 = new TLegend(0.3,0.22,0.55,0.4);
    w40->SetLineColor(kWhite);
    w40->SetFillColor(0);
    w40->SetTextSize(23);
    w40->SetTextFont(45);
    w40->SetNColumns(2);
    w40->SetTextAlign();
    w40->AddEntry(temp5, " ");
    w40->AddEntry(temp9, " ");

    TLatex* r434 = new TLatex(0.50,0.29, "HM correction");
    r434->SetNDC();
    r434->SetTextSize(21);
    r434->SetTextFont(43);

    temp1->Draw("Psame");
	temp11->Draw("Psame");
	right1->Draw("Psame");
	right11->Draw("Psame");
	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");
	w40->Draw("same");
	latex1->Draw("same");
    latex2->Draw("same");
    r434->Draw("same");







	//c2->Print("../results/deltaEtaResults.pdf");



}