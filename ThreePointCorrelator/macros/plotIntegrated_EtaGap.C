#include "RiceStyle.h"

using namespace std;

void plotIntegrated_EtaGap(){

	gStyle->SetErrorX(0);

	TFile* file[9];
	file[0] = new TFile("../dataPoints/PbPb5TeV_60_70_etaGap_data.root");
	file[1] = new TFile("../dataPoints/PbPb5TeV_60_70_etaGap_2p4_data.root");

	TGraphErrors* gr1[2];
	for(int i = 0; i < 2; i++){

		gr1[i] = (TGraphErrors*) file[0]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr2[2];
	for(int i = 0; i < 2; i++){

		gr2[i] = (TGraphErrors*) file[1]->Get(Form("Graph;%d", i+1));
	}
	

//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "#eta gap", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 1000,-1,5.0,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "#eta gap", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (oppo-same)", 1000,-1,5.0,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0009, 0.0007);
	base1->GetXaxis()->SetRangeUser(-1, 5.0);

	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0009, 0.001);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.23);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,5,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base2->GetYaxis()->SetTitleOffset(1.23);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	//c1->Divide(2,1,0.01,0.01);
	//c1->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.12);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base1->Draw();


	gr1[0]->SetMarkerStyle(20);
	gr1[0]->SetMarkerSize(1.4);
	gr1[0]->SetMarkerColor(kRed);
	gr1[0]->SetLineColor(kRed);
	gr1[0]->Draw("Psame");

	gr1[1]->SetMarkerStyle(21);
	gr1[1]->SetMarkerSize(1.4);
	gr1[1]->SetMarkerColor(kBlue);
	gr1[1]->SetLineColor(kBlue);
	gr1[1]->Draw("Psame");


	TLatex* r11 = new TLatex(0.18,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);

    TLatex* r4 = new TLatex(0.18, 0.78, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
   
    TLatex* r5 = new TLatex(0.18, 0.71., "60-70%");
    r5->SetNDC();
    r5->SetTextSize(23);
    r5->SetTextFont(43);
    r5->SetTextColor(kBlack);

    TLatex* r6 = new TLatex(0.18, 0.67., "0.3 < p_{T} < 3.0 GeV/c");
    r6->SetNDC();
    r6->SetTextSize(23);
    r6->SetTextFont(43);
    r6->SetTextColor(kBlack);
	

    r11->Draw("same");
    r22->Draw("same");
    r4->Draw("same");
    r5->Draw("same");
    r6->Draw("same");

	// c1->cd(2);
	// gPad->SetTicks();
	// gPad->SetLeftMargin(0.12);
	// gPad->SetBottomMargin(0.13);
	// gPad->SetRightMargin(0.05);
	// gStyle->SetPadBorderMode(0.1);
	// gStyle->SetOptTitle(0);

	// base2->Draw();

	gr2[0]->SetMarkerStyle(24);
	gr2[0]->SetMarkerSize(1.4);
	gr2[0]->SetMarkerColor(kRed);
	gr2[0]->SetLineColor(kRed);
	gr2[0]->Draw("Psame");

	gr2[1]->SetMarkerStyle(25);
	gr2[1]->SetMarkerSize(1.4);
	gr2[1]->SetMarkerColor(kBlue);
	gr2[1]->SetLineColor(kBlue);
	gr2[1]->Draw("Psame");

	TLegend *w1 = new TLegend(0.53,0.7,0.92,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);

    w1->SetNColumns(2);

    w1->AddEntry(gr1[0], "  ", "P");
    w1->AddEntry(gr1[1], "  |#eta_{#alpha,#beta}|<0.8", "P");

    w1->AddEntry(gr2[0], "  ", "P");
    w1->AddEntry(gr2[1], "  |#eta_{#alpha,#beta}|<2.4", "P");
    w1->Draw("same");
	
	TLatex* latex1 = new TLatex(0.53, 0.82, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.61, 0.82, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");

 //    TLatex* r7 = new TLatex(0.38, 0.20., "0.3 < p_{T} < 3.0 GeV/c, 4.4 < |#eta_{c}| < 5.0");
 //    r7->SetNDC();
 //    r7->SetTextSize(23);
 //    r7->SetTextFont(43);
 //    r7->SetTextColor(kBlack);
 //    r7->Draw("same");

    TGraphErrors* gr1_diff = new TGraphErrors(7);
    TGraphErrors* gr2_diff = new TGraphErrors(4);

    for(int i = 0; i < 7; i++){

    	double x1,y1,ey1;
    	gr1[0]->GetPoint(i, x1, y1);
    	ey1 = gr1[0]->GetErrorY(i);
		
		double x2,y2,ey2;
    	gr1[1]->GetPoint(i, x2, y2);
    	ey2 = gr1[1]->GetErrorY(i);

    	gr1_diff->SetPoint(i, x1, y2-y1);
    	gr1_diff->SetPointError(i, 0,sqrt(ey1*ey1+ey2*ey2));

    }

    for(int i = 0; i < 4; i++){

    	double x1,y1,ey1;
    	gr2[0]->GetPoint(i, x1, y1);
    	ey1 = gr2[0]->GetErrorY(i);
		
		double x2,y2,ey2;
    	gr2[1]->GetPoint(i, x2, y2);
    	ey2 = gr2[1]->GetErrorY(i);

    	gr2_diff->SetPoint(i, x1, y2-y1);
    	gr2_diff->SetPointError(i, 0,sqrt(ey1*ey1+ey2*ey2));

    }


    TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base2->Draw();

	gr1_diff->SetMarkerStyle(20);
	gr1_diff->SetMarkerSize(1.4);
	gr1_diff->SetMarkerColor(kRed);
	gr1_diff->SetLineColor(kRed);
	gr1_diff->Draw("LPsame");

	gr2_diff->SetMarkerStyle(21);
	gr2_diff->SetMarkerSize(1.4);
	gr2_diff->SetMarkerColor(kBlue);
	gr2_diff->SetLineColor(kBlue);
	gr2_diff->Draw("LPsame");

	TLegend *w1 = new TLegend(0.70,0.2,0.88,0.35);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);
    w1->AddEntry(gr1_diff, "|#eta_{#alpha,#beta}|<0.8", "P");
    w1->AddEntry(gr2_diff, "|#eta_{#alpha,#beta}|<2.4", "P");

    w1->Draw("same");
	r11->Draw("same");
    r22->Draw("same");
    r4->Draw("same");
    r5->Draw("same");
    r6->Draw("same");
}