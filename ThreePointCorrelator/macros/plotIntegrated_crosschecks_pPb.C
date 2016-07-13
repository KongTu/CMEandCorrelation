#include "RiceStyle.h"

using namespace std;

//double PbPb_centralityBinCenter[] = {75, 65, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5};
double PbPb_centralityBinCenter[] = {75,67.5,62.5,57.5,52.5,47.5,42.5,37.5,32.5};
double PbPb_centralityBinCenter_tracker[] = {65, 55, 45, 35};

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

void plotIntegrated_crosschecks_pPb(){

	gStyle->SetErrorX(0);

	TFile* file[10];
	file[0] = new TFile("../dataPoints/pPb_5TeV_data_tracker.root");
	file[1] = new TFile("../dataPoints/pPb_5TeV_data_Quan.root");

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

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 5000,0.1,10000,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0003, 0.0006);
	base1->GetXaxis()->SetRangeUser(81, 1700);

	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0003, 0.013);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.23);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,18,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base2->GetYaxis()->SetTitleOffset(1.23);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (OS - SS)");
	base3->GetYaxis()->SetRangeUser(-0.0007,0.0016);
	base3->GetYaxis()->SetTitleOffset(1.23);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.0006,0.0016);
	base4->GetYaxis()->SetTitleOffset(1.9);
	base4->GetXaxis()->SetTitleOffset(3.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);
	

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	base1->GetXaxis()->SetLabelOffset(999);
	base1->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis2 = new TGaxis(81,
	                            -0.0003,
	                            1700,
	                            -0.0003,
	                            81,
	                            1700,
	                            510,"G");
	newaxis2->SetLabelOffset(0.01);
	newaxis2->SetLabelFont(42);
	
	base1->Draw();
	newaxis2->Draw("SS");
	gPad->Update();
	gPad->SetLogx(1);

    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];
    TBox *box5[50];
    TBox *box6[50];

    double xe[9];

    for(int mult = 0; mult < 6; mult++){

    	xe[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe[mult] = 6;
    	double ye = total_systematics_pPb*0.05;

    	double x1;
    	double value1;
    	gr1[0]->GetPoint(mult+3, x1, value1);

    	double x2;
    	double value2;
    	gr1[1]->GetPoint(mult+3, x2, value2);

    	box1[mult] = new TBox(x1-xe[mult],value1-ye,x1+xe[mult],value1+ye);
		box1[mult]->SetFillColor(kRed);
	 	box1[mult]->SetFillColorAlpha(kGray+1,0.3);
        box1[mult]->SetFillStyle(1001);
    	box1[mult]->SetLineWidth(0);
    	box1[mult]->SetLineColor(kRed);
        box1[mult]->Draw("SAME");

		box2[mult] = new TBox(x2-xe[mult],value2-ye,x2+xe[mult],value2+ye);
		box2[mult]->SetFillColor(kBlue);
	 	box2[mult]->SetFillColorAlpha(kGray+1,0.3);
        box2[mult]->SetFillStyle(1001);
    	box2[mult]->SetLineWidth(0);
    	box2[mult]->SetLineColor(kBlue);
        box2[mult]->Draw("SAME");

    }

    double xe2[11]; 

    for(int mult = 0; mult < 9; mult++){

    	xe2[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe2[mult] = 6;

    	double ye = total_systematics_pPb*0.05;

    	double x1;
    	double value1;
    	gr2[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr2[1]->GetPoint(mult, x2, value2);

    	box3[mult] = new TBox(x1-xe2[mult],value1[mult]-ye,x1+xe2[mult],value1[mult]+ye);
		box3[mult]->SetFillColor(kRed);
        box3[mult]->SetFillColorAlpha(kGray+1,0.3);
        box3[mult]->SetFillStyle(1001);
    	box3[mult]->SetLineWidth(0);
    	box3[mult]->SetLineColor(kRed);
        box3[mult]->Draw("SAME");

		box4[mult] = new TBox(x2-xe2[mult],value2[mult]-ye,x2+xe2[mult],value2[mult]+ye);
		box4[mult]->SetFillColor(kBlue);
        box4[mult]->SetFillColorAlpha(kGray+1,0.3);
        box4[mult]->SetFillStyle(1001);
    	box4[mult]->SetLineWidth(0);
    	box4[mult]->SetLineColor(kBlue);
        box4[mult]->Draw("SAME");
    }


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

    TLatex* r4 = new TLatex(0.19, 0.84, "#sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");
	

	TLatex* latex1 = new TLatex(0.48, 0.83, "SS");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    //latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.55, 0.83, "OS");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    //latex2->Draw("same");

	TLegend *w2 = new TLegend(0.45,0.25,0.8,0.5);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(gr1[0], "pPb Kong, SS", "P");
    w2->AddEntry(gr1[1], "pPb Kong, OS", "P");
    w2->AddEntry(gr2[0], "pPb Quan, SS", "P");
    w2->AddEntry(gr2[1], "pPb Quan, OS", "P");
    w2->Draw("same");

    TLatex* r3 = new TLatex(0.69, 0.94, "PbPb centrality(%)");
    r3->SetNDC();
    r3->SetTextSize(21);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* cent1[7];
    cent1[0] = new TLatex(0.55, 0.91, "55");
    cent1[1] = new TLatex(0.70, 0.91, "45");
    cent1[2] = new TLatex(0.83, 0.91, "35");
    cent1[3] = new TLatex(0.35.5, 0.91, "65");
    cent1[4] = new TLatex(0.75, 0.91, "15");
    cent1[5] = new TLatex(0.80, 0.91, "7.5");
    cent1[6] = new TLatex(0.84, 0.91, "2.5");

    for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(13);
	    cent1[i]->SetTextFont(63);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

    TLine* l1[7];
    l1[0] = new TLine(404.1,0.00057, 404.1, 0.0006);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.00057, 717.6, 0.0006);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.00057, 1141, 0.0006);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(81.412,0.00057, 81.412, 0.0006);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(197,0.00057, 197, 0.0006);
    l1[4]->SetLineWidth(2);
    l1[4]->Draw("Lsame");

    l1[5] = new TLine(3577.6,0.00057, 3577.6, 0.0006);
    l1[5]->SetLineWidth(2);
    //l1[5]->Draw("Lsame");

    l1[6] = new TLine(4474,0.00057, 4474, 0.0006);
    l1[6]->SetLineWidth(2);
    //l1[6]->Draw("Lsame");

   	TLatex* r11 = new TLatex(0.8,0.85, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
//    r22->Draw("same");

}