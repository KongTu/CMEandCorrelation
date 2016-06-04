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

void plotDeltaEtaResult(){

	gStyle->SetErrorX(0);

	TFile* file[16];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");
	//file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_MB_v4_4.root");

	//file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_.root");
	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v3_4.root");
	
	file[2] = new TFile("../rootfiles/CME_QvsdEta_pPb_MB_v4_1.root");


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
	TH1D* base2 = makeHist("base2", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012,0.001);
	//base1->GetYaxis()->SetRangeUser(-0.01, 0.03);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0012,0.001);
	//base2->GetYaxis()->SetRangeUser(-0.01, 0.03);

	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT (opposite - same)");
	base3->GetYaxis()->SetRangeUser(-0.0007,0.0018);
	base3->GetYaxis()->SetTitleOffset(1.2);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.0007,0.0018);
	base4->GetYaxis()->SetTitleOffset(1.3);
	base4->GetXaxis()->SetTitleOffset(1.3);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);

	TCanvas* c2 = new TCanvas("c2","c2",1200,500);
	c2->Divide(3,1,0,0);
	c2->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

	TH1D* temp1_plot = (TH1D*)hist1[0][0][0]->Clone("temp1_plot");
	temp1_plot->Add(hist1[0][1][0], +1);
	temp1_plot->Scale(0.5);
	temp1_plot->Scale(1.0/v2[0][0]);
	temp1_plot->SetMarkerStyle(20);
	temp1_plot->SetMarkerSize(1.4);
	temp1_plot->SetMarkerColor(kRed);
	temp1_plot->SetLineColor(kRed);


	TH1D* temp1 = (TH1D*)hist1[0][0][0]->Clone("temp1");
	TH1D* neg1  = (TH1D*)hist1[0][1][0]->Clone("neg1");
	neg1->Scale(1.0/v2[0][0]);
	neg1->SetMarkerStyle(24);
	neg1->SetMarkerSize(1.4);
	neg1->SetMarkerColor(kBlack);
	neg1->SetLineColor(kBlack);
	
	temp1->Scale(1.0/v2[0][0]);
	temp1->SetMarkerStyle(20);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	TH1D* temp2 = (TH1D*) hist1[0][2][0]->Clone("temp2");
	temp2->SetMarkerStyle(21);
	temp2->Scale(1.0/v2[0][0]);
	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);


	TH1D* temp3_plot = (TH1D*)hist1[0][0][1]->Clone("temp3_plot");
	temp3_plot->Add(hist1[0][1][1], +1);
	temp3_plot->Scale(0.5);
	temp3_plot->Scale(1.0/v2[0][1]);
	temp3_plot->SetMarkerStyle(20);
	temp3_plot->SetMarkerSize(1.4);
	temp3_plot->SetMarkerColor(kRed);
	temp3_plot->SetLineColor(kRed);
	
	TH1D* temp3 = (TH1D*)hist1[0][0][1]->Clone("temp3");
	TH1D* neg2 = (TH1D*) hist1[0][1][1]->Clone("neg2");
	neg2->Scale(1.0/v2[0][1]);
	neg2->SetMarkerStyle(24);
	neg2->SetMarkerSize(1.4);
	neg2->SetMarkerColor(kBlack);
	neg2->SetLineColor(kBlack);
	//temp3->Add(hist1[0][1][1], +1);
	//temp3->Scale(0.5);
	temp3->Scale(1.0/v2[0][1]);
	temp3->SetMarkerStyle(20);
	temp3->SetMarkerSize(1.4);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);

	TH1D* temp4 = (TH1D*) hist1[0][2][1]->Clone("temp4");
	temp4->SetMarkerStyle(21);
	temp4->Scale(1.0/v2[0][1]);	
	temp4->SetMarkerSize(1.4);
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);

    TLatex* r3 = new TLatex(0.25, 0.84, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r3->SetNDC();
    r3->SetTextSize(23);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* lmult = new TLatex(0.25, 0.76, "260 #leq N^{offline}_{trk} < 300");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

    double value1[50];
    double value1_error[50];
    double value2[50];
    double value2_error[50];
    double value3[50];
    double value3_error[50];
    double value4[50];
    double value4_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value1[deta] = temp1->GetBinContent(deta+1);
    	value1_error[deta] = temp1->GetBinError(deta+1);

    	value2[deta] = temp2->GetBinContent(deta+1);
    	value2_error[deta] = temp2->GetBinError(deta+1);

    	value3[deta] = temp3->GetBinContent(deta+1);
    	value3_error[deta] = temp3->GetBinError(deta+1);

    	value4[deta] = temp4->GetBinContent(deta+1);
    	value4_error[deta] = temp4->GetBinError(deta+1);
    }

    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];

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
    }

	temp1->Draw("Psame");
	temp2->Draw("Psame");
	neg1->Draw("Psame");

    TLatex* r6 = new TLatex(0.73, 0.2, "Pb-going");
    r6->SetNDC();
    r6->SetTextSize(23);
    r6->SetTextFont(43);
    r6->SetTextColor(kBlack);
    r6->Draw("same");

//PAD1:
	c2->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base2->Draw();

    TBox *box3[50];
    TBox *box4[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.065;
    	double ye = total_systematics_pPb;

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

    temp3->Draw("Psame");
	temp4->Draw("Psame");
	neg2->Draw("Psame");

    TLatex* r32 = new TLatex(0.11, 0.84, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r32->SetNDC();
    r32->SetTextSize(23);
    r32->SetTextFont(43);
    r32->SetTextColor(kBlack);
    r32->Draw("same");

    TLatex* r62 = new TLatex(0.73, 0.2, "p-going");
    r62->SetNDC();
    r62->SetTextSize(23);
    r62->SetTextFont(43);
    r62->SetTextColor(kBlack);
    r62->Draw("same");

	c2->cd(3);
	gPad->SetBottomMargin(0.13);
	gPad->SetTicks();
	gPad->SetTopMargin(0.06);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	base2->Draw();

	TH1D* temp5_plot = (TH1D*)hist1[1][0][0]->Clone("temp5_plot");
	temp5_plot->Add(hist1[1][0][1], +1);
	temp5_plot->Add(hist1[1][1][0], +1);
	temp5_plot->Add(hist1[1][1][1], +1);

	temp5_plot->Scale(0.25);
	temp5_plot->Scale(1.0/v2[1][2]);
	temp5_plot->SetMarkerStyle(20);
	temp5_plot->SetMarkerSize(1.4);
	temp5_plot->SetMarkerColor(kRed);
	temp5_plot->SetLineColor(kRed);

	TH1D* temp5 = (TH1D*)hist1[1][0][0]->Clone("temp5");
	TH1D* temp6 = (TH1D*)hist1[1][1][0]->Clone("temp6");
	TH1D* temp7 = (TH1D*)hist1[1][0][1]->Clone("temp6");
	TH1D* temp8 = (TH1D*)hist1[1][1][1]->Clone("temp6");

	temp5->Add(temp7, +1);
	temp5->Scale(0.5);
	temp5->Scale(1.0/v2[1][2]);
	temp5->SetMarkerStyle(20);
	temp5->SetMarkerSize(1.4);
	temp5->SetMarkerColor(kRed);
	temp5->SetLineColor(kRed);
	
	temp6->Add(temp8, +1);
	temp6->Scale(0.5);
	temp6->Scale(1.0/v2[1][2]);
	temp6->SetMarkerStyle(24);
	temp6->SetMarkerSize(1.4);
	temp6->SetMarkerColor(kBlack);
	temp6->SetLineColor(kBlack);

	TH1D* temp9 = (TH1D*) hist1[1][2][0]->Clone("temp9");
	TH1D* temp10 = (TH1D*) hist1[1][2][1]->Clone("temp10");
	temp9->Add(temp10, +1);
	temp9->Scale(0.5);
	temp9->Scale(1.0/v2[1][2]);
	temp9->SetMarkerStyle(21);
	temp9->SetMarkerColor(kBlue);
	temp9->SetMarkerSize(1.4);
	temp9->SetLineColor(kBlue);


	TLegend *w2 = new TLegend(0.60,0.15,0.75,0.3);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(45);
    w2->AddEntry(temp5, " pos-pos ");
    w2->AddEntry(temp6, " neg-neg ");
    w2->AddEntry(temp9, " pos-neg ");
    w2->Draw("same");

   	TLatex* r1 = new TLatex(0.58,0.95, "CMS");
    r1->SetNDC();
    r1->SetTextSize(0.06);
    r1->Draw("same");

    TLatex* r2 = new TLatex(0.72,0.95, "Preliminary");
    r2->SetNDC();
    r2->SetTextSize(22);
    r2->SetTextFont(53);
    r2->Draw("same");

    TLatex* r4 = new TLatex(0.11, 0.84, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");

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

    TBox *box5[50];
    TBox *box6[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.065;
    	double ye = total_systematics_PbPb;

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

    temp5->Draw("Psame");
	temp6->Draw("Psame");
	temp9->Draw("Psame");

//plotting difference between unlike and like sign:

    TCanvas* c3 = new TCanvas("c3","c3", 600, 600);
 	gPad->SetLeftMargin(0.16);
 	gPad->SetBottomMargin(0.13);
 	gPad->SetTopMargin(0.1);
	gPad->SetTicks();
	base3->Draw();

	TH1D* diff1 = (TH1D*) temp2->Clone("diff1");
	diff1->SetMarkerStyle(20);
	diff1->SetMarkerSize(1.3);
	diff1->SetMarkerColor(kRed);
	diff1->SetLineColor(kRed);
	diff1->Add(temp1, -1);
	
	TH1D* diff2 = (TH1D*) temp4->Clone("diff2");
	diff2->SetMarkerStyle(21);
	diff2->SetMarkerSize(1.3);
	diff2->SetMarkerColor(kBlue);
	diff2->SetLineColor(kBlue);
	diff2->Add(temp3, -1);

	diff1->Draw("Psame");
	diff2->Draw("Psame");

    TLatex* r33 = new TLatex(0.25, 0.82, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    //r33->Draw("same");

    TLatex* lmult = new TLatex(0.20, 0.82, "185 #leq N^{offline}_{trk} < 220");
    lmult->SetNDC();
    lmult->SetTextSize(26);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

	TH1D* diff3 = (TH1D*) temp9->Clone("diff3");
	diff3->Add(temp5, -1);
	diff3->SetMarkerStyle(28);
	diff3->SetMarkerSize(1.4);
	diff3->SetMarkerColor(kBlack);
	diff3->SetLineColor(kBlack);
	diff3->Draw("Psame");

   	TLatex* r11 = new TLatex(0.62,0.91, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.71,0.91, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

	TLegend *w3 = new TLegend(0.40,0.5,0.6,0.7);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(diff1, "pPb #sqrt{s_{NN}} = 5.02 TeV, Pb-going");
    w3->AddEntry(diff2, "pPb #sqrt{s_{NN}} = 5.02 TeV, p-going");
    w3->AddEntry(diff3, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    w3->Draw("same");

	TH1D* base5 = (TH1D*) base1->Clone("base3");
	base5->GetYaxis()->SetTitleOffset(1.3);
	base5->GetXaxis()->SetTitleOffset(0.95);
	base5->GetYaxis()->SetTitleSize(base5->GetYaxis()->GetTitleSize()*1.3);
	base5->GetXaxis()->SetTitleSize(base5->GetXaxis()->GetTitleSize()*1.3);
	base5->GetYaxis()->SetLabelSize(base5->GetYaxis()->GetLabelSize()*1.3);
	base5->GetXaxis()->SetLabelSize(base5->GetXaxis()->GetLabelSize()*1.3);
	base5->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base6 = (TH1D*) base5->Clone("base6");
	base6->GetYaxis()->SetRangeUser(-0.01, 0.04);
	base6->GetXaxis()->SetRangeUser(-0.2, 4.8);

    TLatex* r41 = new TLatex(0.24, 0.85, "pPb #sqrt{s_{NN}} = 5.02 TeV, Pb-going");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.05, 0.85, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
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
    
    TLatex* r45 = new TLatex(0.24, 0.78, "185 #leq N^{offline}_{trk} < 220");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);
    
    
    TCanvas* c4 = new TCanvas("c4","c4",1000,600);
	c4->Divide(2,1,0,0);
	c4->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base5->Draw();

    double value1[50];
    double value1_error[50];
    double value2[50];
    double value2_error[50];
    double value3[50];
    double value3_error[50];
    double value4[50];
    double value4_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value1[deta] = temp1_plot->GetBinContent(deta+1);
    	value1_error[deta] = temp1_plot->GetBinError(deta+1);

    	value2[deta] = temp2->GetBinContent(deta+1);
    	value2_error[deta] = temp2->GetBinError(deta+1);

    	value3[deta] = temp3_plot->GetBinContent(deta+1);
    	value3_error[deta] = temp3_plot->GetBinError(deta+1);

    	value4[deta] = temp4->GetBinContent(deta+1);
    	value4_error[deta] = temp4->GetBinError(deta+1);
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
    }

    temp1_plot->Draw("Psame");
	temp2->Draw("Psame");
	r41->Draw("same");
	r45->Draw("same");

    c4->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base5->Draw();

    double value5[50];
    double value5_error[50];
    double value6[50];
    double value6_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value5[deta] = temp5_plot->GetBinContent(deta+1);
    	value5_error[deta] = temp5_plot->GetBinError(deta+1);

    	value6[deta] = temp9->GetBinContent(deta+1);
    	value6_error[deta] = temp9->GetBinError(deta+1);

    }	

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.065;
    	double ye = total_systematics_PbPb;

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

    temp5_plot->Draw("Psame");
	temp9->Draw("Psame");
	r42->Draw("same");
	r43->Draw("same");
	r44->Draw("same");

	TLegend *w4 = new TLegend(0.60,0.25,0.75,0.4);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->AddEntry(temp5_plot, "  same");
    w4->AddEntry(temp9, "  opposite");
    w4->Draw("same");

    TLatex* r433 = new TLatex(0.47,0.94, "CMS");
    r433->SetNDC();
    r433->SetTextSize(0.06);
    
    TLatex* r444 = new TLatex(0.62,0.94, "Preliminary");
    r444->SetNDC();
    r444->SetTextSize(21);
    r444->SetTextFont(53);
    
    TLatex* r333 = new TLatex(0.25, 0.84, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r333->SetNDC();
    r333->SetTextSize(21);
    r333->SetTextFont(43);
    r333->SetTextColor(kBlack);

    TLatex* r455 = new TLatex(0.25, 0.74, "0 #leq N^{offline}_{trk} < 35");
    r455->SetNDC();
    r455->SetTextSize(21);
    r455->SetTextFont(43);
    r455->SetTextColor(kBlack);

	TLatex* r755 = new TLatex(0.25, 0.82, "185 #leq N^{offline}_{trk} < 220");
    r755->SetNDC();
    r755->SetTextSize(22);
    r755->SetTextFont(43);
    r755->SetTextColor(kBlack);

//low-multiplicity

    TH1D* base7 = (TH1D*) base6->Clone("base7");
    base7->GetYaxis()->SetTitleOffset(1.9);
	base7->GetXaxis()->SetTitleOffset(1.3);
	base7->GetYaxis()->SetTitleSize(base7->GetYaxis()->GetTitleSize()*0.7);
	base7->GetXaxis()->SetTitleSize(base7->GetXaxis()->GetTitleSize()*0.7);
	base7->GetYaxis()->SetLabelSize(base7->GetYaxis()->GetLabelSize()*0.7);
	base7->GetXaxis()->SetLabelSize(base7->GetXaxis()->GetLabelSize()*0.7);
	base7->GetXaxis()->SetNdivisions(5,6,0);
	base7->GetYaxis()->SetNdivisions(5,6,0);

    TH1D* base8 = (TH1D*) base5->Clone("base8");
    base8->GetYaxis()->SetRangeUser(-0.0011,0.0012);
    base8->GetYaxis()->SetTitleOffset(1.9);
	base8->GetXaxis()->SetTitleOffset(1.3);
	base8->GetYaxis()->SetTitleSize(base8->GetYaxis()->GetTitleSize()*0.7);
	base8->GetXaxis()->SetTitleSize(base8->GetXaxis()->GetTitleSize()*0.7);
	base8->GetYaxis()->SetLabelSize(base8->GetYaxis()->GetLabelSize()*0.75);
	base8->GetXaxis()->SetLabelSize(base8->GetXaxis()->GetLabelSize()*0.75);
	base8->GetXaxis()->SetNdivisions(5,6,0);

	TLatex* r7 = new TLatex(0.72, 0.2, "Pb-going");
    r7->SetNDC();
    r7->SetTextSize(23);
    r7->SetTextFont(43);
    r7->SetTextColor(kBlack);

	TLatex* r72 = new TLatex(0.67, 0.2, "p-going");
    r72->SetNDC();
    r72->SetTextSize(23);
    r72->SetTextFont(43);
    r72->SetTextColor(kBlack);


	TH1D* temp11 = (TH1D*)hist1[2][0][0]->Clone("temp11");
	temp11->Add(hist1[2][1][0], +1);
	temp11->Scale(0.5);
	temp11->Scale(1.0/v2[2][0]);
	temp11->SetMarkerStyle(20);
	temp11->SetMarkerSize(1.4);
	temp11->SetMarkerColor(kRed);
	temp11->SetLineColor(kRed);

	TH1D* temp12 = (TH1D*) hist1[2][2][0]->Clone("temp12");
	temp12->SetMarkerStyle(21);
	temp12->Scale(1.0/v2[2][0]);
	temp12->SetMarkerColor(kBlue);
	temp12->SetMarkerSize(1.4);
	temp12->SetLineColor(kBlue);

	TH1D* temp13 = (TH1D*)hist1[2][0][1]->Clone("temp13");
	temp13->Add(hist1[2][1][1], +1);
	temp13->Scale(0.5);
	temp13->Scale(1.0/v2[2][1]);
	temp13->SetMarkerStyle(20);
	temp13->SetMarkerSize(1.4);
	temp13->SetMarkerColor(kRed);
	temp13->SetLineColor(kRed);

	TH1D* temp14 = (TH1D*) hist1[2][2][1]->Clone("temp14");
	temp14->SetMarkerStyle(21);
	temp14->Scale(1.0/v2[2][1]);	
	temp14->SetMarkerSize(1.4);	
	temp14->SetMarkerColor(kBlue);
	temp14->SetLineColor(kBlue);


    TCanvas* c5 = new TCanvas("c5","c5",700,700);
	c5->Divide(2,2,0,0);
	c5->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.115);
	gPad->SetTopMargin(0.07);
	gPad->SetTicks();
	base7->Draw();
	temp11->Draw("Psame");
	temp12->Draw("Psame");
	r455->Draw("same");
	r333->Draw("same");
    r7->Draw("same");

	c5->cd(2);
	gPad->SetBottomMargin(0.115);
	gPad->SetTopMargin(0.07);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base7->Draw();
	temp13->Draw("Psame");
	temp14->Draw("Psame");
	r433->Draw("same");
	r444->Draw("same");
    r72->Draw("same");
	

	TLegend *w5 = new TLegend(0.55,0.60,0.7,0.80);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(23);
    w5->SetTextFont(45);
    w5->AddEntry(temp11, "  same");
    w5->AddEntry(temp12, "  opposite");
    w5->Draw("same");

    c5->cd(3);
    gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.09);
	gPad->SetTicks();
	base8->Draw();
    r7->Draw("same");
    r755->Draw("same");

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

    }

	temp1_plot->Draw("Psame");
	temp2->Draw("Psame");

    c5->cd(4);
	gPad->SetBottomMargin(0.15);
	gPad->SetTopMargin(0.09);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base8->Draw();
    r72->Draw("same");

    for(int deta = 0; deta < NdEtaReBins2; deta++){

		double xe = 0.065;
		double ye = total_systematics_pPb;

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


	temp3_plot->Draw("Psame");
	temp4->Draw("Psame");

	


	TFile t1("../dataPoints/diff.root","RECREATE");
	diff1->Write();
	diff2->Write();
	diff3->Write();
	temp1->Write();
	temp2->Write();
	temp3->Write();
	temp4->Write();
	temp5->Write();
	temp9->Write();
//


	//c2->Print("../results/deltaEtaResults.pdf");

	return;

}