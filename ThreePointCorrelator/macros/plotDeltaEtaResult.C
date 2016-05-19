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
const int Nmults = 2;

void plotDeltaEtaResult(){

	TFile* file[16];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3and4.root");
	//file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");

	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_7and8.root");
	// file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_50_100_v3_7.root");

	TH1D* QvsdEta[16][48][3][2];

	TH1D* delEta3p[16][3][2];

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
	TH1D* hist2[8][3][2];
	for(int mult = 0; mult < Nmults; mult++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
				hist1[mult][sign][HF] = new TH1D(Form("hist1_%d_%d_%d",mult,sign,HF),"test", NdEtaBins, dEtaBins);
				hist2[mult][sign][HF] = new TH1D(Form("hist2_%d_%d_%d",mult,sign,HF),"test", NdEtaBins, dEtaBins);
			}
		}
	}

	for(int mult = 0; mult < Nmults; mult++){
		for(int deta = 0; deta < NdEtaBins; deta++){
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

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8,kBlack);
	TH1D* base2 = makeHist("base2", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0014, 0.0018);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0014, 0.0018);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.9);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.9);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#Delta#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT");
	base3->GetYaxis()->SetRangeUser(-0.001,0.002);
	base3->GetYaxis()->SetTitleOffset(1.9);
	base3->GetXaxis()->SetTitleOffset(3.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.0);
	base3->GetYaxis()->SetNdivisions(6);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.001,0.002);
	base4->GetYaxis()->SetTitleOffset(1.9);
	base4->GetXaxis()->SetTitleOffset(3.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);

	TCanvas *c2 = new TCanvas("c2","c2",1,1,1000,800);
	c2->Range(0,0,1,1);
	TPad* pad2[8];

	pad2[0] = new TPad("pad20", "pad20",0.0,   0.33, 0.53,   1);
	pad2[1] = new TPad("pad21", "pad21",0.53,   0.33, 1.0, 1);

	pad2[2] = new TPad("pad28", "pad28",0.0,      0.01, 0.53,   0.32);
	pad2[3] = new TPad("pad29", "pad29",0.53,      0.01, 1.0, 0.32);


	for(int i = 0; i < 4; i++){

	pad2[i]->SetLeftMargin(0.0);
	pad2[i]->SetRightMargin(0);
	pad2[i]->SetTopMargin(0.0);
	pad2[i]->SetBottomMargin(0);
	pad2[i]->Draw();

	}

	pad2[0]->SetLeftMargin(0.22);
	pad2[2]->SetLeftMargin(0.22);

	pad2[1]->SetRightMargin(0.1);
	pad2[3]->SetRightMargin(0.1);

	pad2[0]->SetTopMargin(0.1);
	pad2[1]->SetTopMargin(0.1);

	pad2[2]->SetBottomMargin(0.28);
	pad2[3]->SetBottomMargin(0.28);

	pad2[2]->SetTopMargin(0.1);
	pad2[3]->SetTopMargin(0.1);

	pad2[0]->cd();
	pad2[0]->SetTicks();
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	base1->Draw();

	TH1D* temp1 = (TH1D*)hist1[0][0][0]->Clone("temp1");
	temp1->Add(hist1[0][1][0], +1);
	temp1->Scale(0.5);
	temp1->Scale(1.0/v2[0][0]);
	temp1->SetMarkerStyle(24);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	temp1->Draw("Psame");

	TH1D* temp2 = (TH1D*) hist1[0][2][0]->Clone("temp2");
	temp2->SetMarkerStyle(25);
	temp2->Scale(1.0/v2[0][0]);
	temp2->SetMarkerColor(kBlue);
	temp2->SetLineColor(kBlue);

	temp2->Draw("Psame");

	TH1D* temp3 = (TH1D*)hist1[0][0][1]->Clone("temp3");
	temp3->Add(hist1[0][1][1], +1);
	temp3->Scale(0.5);
	temp3->Scale(1.0/v2[0][1]);
	temp3->SetMarkerStyle(20);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);

	temp3->Draw("Psame");

	TH1D* temp4 = (TH1D*) hist1[0][2][1]->Clone("temp4");
	temp4->SetMarkerStyle(21);
	temp4->Scale(1.0/v2[0][1]);	
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);

	temp4->Draw("Psame");

    TLatex* r3 = new TLatex(0.6, 0.82, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r3->SetNDC();
    r3->SetTextSize(23);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* lmult = new TLatex(0.6, 0.74, "185 #leq N^{offline}_{trk} < 260");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

	TLegend *w1 = new TLegend(0.55,0.1,0.7,0.3);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(temp1, "Pb-going, like sign");
    w1->AddEntry(temp2, "Pb-going, unlike sign");
    w1->AddEntry(temp3, "p-going, like sign");
    w1->AddEntry(temp4, "p-going, unlike sign");
    w1->Draw("same");

	pad2[1]->cd();
	pad2[1]->SetTicks();
	base2->Draw();

	TH1D* temp5 = (TH1D*)hist1[1][0][0]->Clone("temp5");
	TH1D* temp6 = (TH1D*)hist1[1][1][0]->Clone("temp6");
	TH1D* temp7 = (TH1D*)hist1[1][0][1]->Clone("temp6");
	TH1D* temp8 = (TH1D*)hist1[1][1][1]->Clone("temp6");

	temp5->Add(temp6, +1);
	temp5->Add(temp7, +1);
	temp5->Add(temp8, +1);

	temp5->Scale(0.25);
	temp5->Scale(1.0/v2[1][2]);
	temp5->SetMarkerStyle(24);
	temp5->SetMarkerColor(kRed);
	temp5->SetLineColor(kRed);

	temp5->Draw("Psame");

	TH1D* temp9 = (TH1D*) hist1[1][2][0]->Clone("temp9");
	TH1D* temp10 = (TH1D*) hist1[1][2][1]->Clone("temp10");
	temp9->Add(temp10, +1);
	temp9->Scale(0.5);
	temp9->Scale(1.0/v2[1][2]);
	temp9->SetMarkerStyle(25);
	temp9->SetMarkerColor(kBlue);
	temp9->SetLineColor(kBlue);

	temp9->Draw("Psame");

	TLegend *w2 = new TLegend(0.55,0.1,0.7,0.3);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp5, "like sign");
    w2->AddEntry(temp9, "unlike sign");
    w2->Draw("same");

   	TLatex* r1 = new TLatex(0.56,0.92, "CMS");
    r1->SetNDC();
    r1->SetTextSize(0.05);
    r1->Draw("same");

    TLatex* r2 = new TLatex(0.69,0.92, "Preliminary");
    r2->SetNDC();
    r2->SetTextSize(22);
    r2->SetTextFont(53);
    r2->Draw("same");

    TLatex* r4 = new TLatex(0.42, 0.82, "PbPb #sqrt{s_{NN}} = 2.76 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");	

	pad2[2]->cd();
	pad2[2]->SetTicks();
	base3->Draw();

	TH1D* diff1 = (TH1D*) temp2->Clone("diff1");
	diff1->Add(temp1, -1);
	
	TH1D* diff2 = (TH1D*) temp4->Clone("diff2");
	diff2->Add(temp3, -1);

	TLegend *w3 = new TLegend(0.55,0.65,0.7,0.8);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(diff1, "unlike-like sign, Pb-going");
    w3->AddEntry(diff2, "unlike-like sign, p-going");
    w3->Draw("same");


	diff1->Draw("Psame");
	diff2->Draw("Psame");

	pad2[3]->cd();
	pad2[3]->SetTicks();
	base4->Draw();

	TH1D* diff3 = (TH1D*) temp9->Clone("diff3");
	diff3->Add(temp5, -1);

	diff3->Draw("Psame");


	return;

}