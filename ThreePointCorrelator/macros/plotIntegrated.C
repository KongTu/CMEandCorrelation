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

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};
double PbPb_5TeV_ntrkBinCenter[] = {103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};

double PbPb_ntrkCentralityBinCenter[] = {625.275, 811.118, 1025.79, 1257.64};

double PbPb_centralityBinCenter[] = {57.5, 52.5, 47.5, 42.5, 37.5, 32.5};

const int Nmults = 23;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

void plotIntegrated(){

	gStyle->SetErrorX(0);

	TFile* file[6];
	file[0] = new TFile("../dataPoints/pPb_data.root");
	file[1] = new TFile("../dataPoints/PbPb5TeV_data.root");
	file[2] = new TFile("../dataPoints/PbPb_5TeV_centrality_data.root");
	file[3] = new TFile("../dataPoints/EPOS_data.root");
	file[4] = new TFile("../dataPoints/ALICE_STAR_data.root");
	file[5] = new TFile("../dataPoints/PbPb_5TeV_centrality_centrality_data.root");

	TGraphErrors* gr1[4];
	for(int i = 0; i < 4; i++){

		gr1[i] = (TGraphErrors*) file[0]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr2[3];
	for(int i = 0; i < 3; i++){

		gr2[i] = (TGraphErrors*) file[1]->Get(Form("Graph;%d", i+1));
	}
	
	TGraphErrors* gr3[2];
	for(int i = 0; i < 2; i++){

		gr3[i] = (TGraphErrors*) file[2]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr4[4];
	for(int i = 0; i < 4; i++){

		gr4[i] = (TGraphErrors*) file[3]->Get(Form("Graph;%d", i+1));
	}
	
	TGraphErrors* gr5[6];
	for(int i = 0; i < 6; i++){

		gr5[i] = (TGraphErrors*) file[4]->Get(Form("Graph;%d", i+1));
	}
	TGraphErrors* gr6[2];
	for(int i = 0; i < 2; i++){

		gr6[i] = (TGraphErrors*) file[5]->Get(Form("Graph;%d", i+1));
	}


	for(int i = 0; i < 4; i++){
		for(int mult = 0; mult < 3; mult++){
			double x, y;
			gr1[i]->GetPoint(mult, x, y);
			gr1[i]->SetPoint(mult, x, 100);
		}
	}
	// for(int i = 0; i < 2; i++){
	// 	for(int mult = 0; mult < 3; mult++){
	// 		double x, y;
	// 		gr2[i]->GetPoint(mult, x, y);
	// 		gr2[i]->SetPoint(mult, x, 100);
	// 	}
	// }
	// 
	

//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT", 5000,0.1,10000,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0009, 0.0008);
	base1->GetXaxis()->SetRangeUser(70, 2000);

	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.001, 0.015);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,18,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	base2->GetYaxis()->SetTitleOffset(1.9);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetRangeUser(-0.0006,0.0012);
	base3->GetYaxis()->SetTitleOffset(1.9);
	base3->GetXaxis()->SetTitleOffset(3.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.0);
	base3->GetYaxis()->SetNdivisions(6);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.0006,0.0012);
	base4->GetYaxis()->SetTitleOffset(1.9);
	base4->GetXaxis()->SetTitleOffset(3.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);

	TCanvas* c1 = new TCanvas("c1","c1",1,1,1300,650);
	c1->Divide(2,1,0.01,0.01);
	c1->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	base1->Draw();

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
    	double ye = total_systematics_pPb;

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

    double xe2[9]; 

    for(int mult = 0; mult < 9; mult++){

    	xe2[mult] = 7*log(1.9*(mult+1));
    	if(mult == 0) xe2[mult] = 6;

    	double ye = total_systematics_PbPb;

    	double x1;
    	double value1;
    	gr2[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr2[1]->GetPoint(mult, x2, value2);


    	box3[mult] = new TBox(PbPb_5TeV_ntrkBinCenter[mult]-xe2[mult],value1[mult]-ye,PbPb_5TeV_ntrkBinCenter[mult]+xe2[mult],value1[mult]+ye);
		box3[mult]->SetFillColor(kRed);
        box3[mult]->SetFillColorAlpha(kGray+1,0.3);
        box3[mult]->SetFillStyle(1001);
    	box3[mult]->SetLineWidth(0);
    	box3[mult]->SetLineColor(kRed);
        box3[mult]->Draw("SAME");

		box4[mult] = new TBox(PbPb_5TeV_ntrkBinCenter[mult]-xe2[mult],value2[mult]-ye,PbPb_5TeV_ntrkBinCenter[mult]+xe2[mult],value2[mult]+ye);
		box4[mult]->SetFillColor(kBlue);
        box4[mult]->SetFillColorAlpha(kGray+1,0.3);
        box4[mult]->SetFillStyle(1001);
    	box4[mult]->SetLineWidth(0);
    	box4[mult]->SetLineColor(kBlue);
        box4[mult]->Draw("SAME");
    }

    double xe3[4];
    for(int mult = 0; mult < 4; mult++){

    	xe3[mult] = 30*log(1.9*(mult+1));
    	if(mult == 0) xe3[mult] = 35;
    	double ye = total_systematics_PbPb;

    	double x1;
    	double value1;
    	gr3[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr3[1]->GetPoint(mult, x2, value2);


    	box5[mult] = new TBox(PbPb_ntrkCentralityBinCenter[mult]-xe3[mult],value1[mult]-ye,PbPb_ntrkCentralityBinCenter[mult]+xe3[mult],value1[mult]+ye);
		box5[mult]->SetFillColor(kRed);
        box5[mult]->SetFillColorAlpha(kGray+1,0.3);
        box5[mult]->SetFillStyle(1001);
    	box5[mult]->SetLineWidth(0);
    	box5[mult]->SetLineColor(kRed);
        box5[mult]->Draw("SAME");

		box6[mult] = new TBox(PbPb_ntrkCentralityBinCenter[mult]-xe3[mult],value2[mult]-ye,PbPb_ntrkCentralityBinCenter[mult]+xe3[mult],value2[mult]+ye);
		box6[mult]->SetFillColor(kBlue);
        box6[mult]->SetFillColorAlpha(kGray+1,0.3);
        box6[mult]->SetFillStyle(1001);
    	box6[mult]->SetLineWidth(0);
    	box6[mult]->SetLineColor(kBlue);
        box6[mult]->Draw("SAME");
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

	// gr1[2]->SetMarkerStyle(20);
	// gr1[2]->SetMarkerColor(kRed);
	// gr1[2]->SetLineColor(kRed);
	// gr1[2]->Draw("Psame");

	// gr1[3]->SetMarkerStyle(21);
	// gr1[3]->SetMarkerColor(kBlue);
	// gr1[3]->SetLineColor(kBlue);
	// gr1[3]->Draw("Psame");

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

	gr2[2]->SetMarkerStyle(25);
	gr2[2]->SetMarkerSize(1.4);
	gr2[2]->SetMarkerColor(kBlue);
	gr2[2]->SetLineColor(kBlue);
	gr2[2]->SetFillColorAlpha(kBlue,0.2);
   	gr2[2]->SetFillStyle(1001);
	//gr2[2]->DrawClone("E3same");

	gr3[0]->SetMarkerStyle(24);
	gr3[0]->SetMarkerSize(1.4);
	gr3[0]->SetMarkerColor(kRed);
	gr3[0]->SetLineColor(kRed);
	gr3[0]->Draw("Psame");

	gr3[1]->SetMarkerStyle(25);
	gr3[1]->SetMarkerSize(1.4);
	gr3[1]->SetMarkerColor(kBlue);
	gr3[1]->SetLineColor(kBlue);
	gr3[1]->Draw("Psame");

	gr4[0]->SetMarkerStyle(24);
	gr4[0]->SetMarkerSize(1.4);
	gr4[0]->SetMarkerColor(kBlack);
	gr4[0]->SetLineColor(kBlack);
	gr4[0]->SetLineStyle(5);
	gr4[0]->SetLineWidth(2.0);
	gr4[0]->Draw("same");

	gr4[1]->SetMarkerStyle(25);
	gr4[1]->SetMarkerSize(1.4);
	gr4[1]->SetMarkerColor(kGreen+2);
	gr4[1]->SetLineColor(kGreen+2);
	gr4[1]->SetLineStyle(5);
	gr4[1]->SetLineWidth(2.0);
	gr4[1]->Draw("same");

	gr5[0]->SetMarkerStyle(26);
	gr5[0]->SetMarkerSize(1.4);
	gr5[0]->SetMarkerColor(kBlack);
	gr5[0]->SetLineColor(kBlack);
	//gr5[0]->Draw("Psame");

	gr5[1]->SetMarkerStyle(27);
	gr5[1]->SetMarkerSize(1.4);
	gr5[1]->SetMarkerColor(kBlack);
	gr5[1]->SetLineColor(kBlack);
	//gr5[1]->Draw("Psame");

	gr5[2]->SetMarkerStyle(28);
	gr5[2]->SetMarkerSize(1.4);
	gr5[2]->SetMarkerColor(kBlack);
	gr5[2]->SetLineColor(kBlack);
	//gr5[2]->Draw("Psame");

	gr5[3]->SetMarkerStyle(30);
	gr5[3]->SetMarkerSize(1.4);
	gr5[3]->SetMarkerColor(kBlack);
	gr5[3]->SetLineColor(kBlack);
	//gr5[3]->Draw("Psame");
	
	gr6[0]->SetMarkerStyle(24);
	gr6[0]->SetMarkerSize(1.4);
	gr6[0]->SetMarkerColor(kRed);
	gr6[0]->SetLineColor(kRed);
	//gr6[0]->Draw("Psame");

	gr6[1]->SetMarkerStyle(25);
	gr6[1]->SetMarkerSize(1.4);
	gr6[1]->SetMarkerColor(kBlue);
	gr6[1]->SetLineColor(kBlue);
	//gr6[1]->Draw("Psame");

	TLegend *w1 = new TLegend(0.48,0.64,0.88,0.77);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);
    //w1->SetHeader("CMS");

    w1->SetNColumns(2);

    w1->AddEntry(gr1[0], "  ", "P");
    w1->AddEntry(gr1[1], "  pPb, Pb-going", "P");

    w1->AddEntry(gr2[0], "  ", "P");
    w1->AddEntry(gr2[1], "  PbPb", "P");

    //w1->AddEntry(gr5[0], "  ", "P");
    //w1->AddEntry(gr5[1], "ALICE PbPb #sqrt{s_{NN}} = 2.76 TeV", "P");

	//w1->AddEntry(gr5[2], "  ", "P");
    //w1->AddEntry(gr5[3], "STAR AuAu #sqrt{s_{NN}} = 200 GeV", "P");

    w1->Draw("same");

    TLatex* r4 = new TLatex(0.18, 0.78, "#sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");
	

	TLatex* latex1 = new TLatex(0.47, 0.78, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.56, 0.78, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");

	TLegend *w2 = new TLegend(0.45,0.15,0.85,0.25);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(gr4[1], "EPOS pPb, oppo", "L");
    w2->AddEntry(gr4[0], "EPOS pPb, same", "L");
    w2->Draw("same");

    TLatex* r3 = new TLatex(0.63, 0.84, "PbPb centrality(%)");
    r3->SetNDC();
    r3->SetTextSize(21);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* cent1[7];
    cent1[0] = new TLatex(0.54, 0.91, "55");
    cent1[1] = new TLatex(0.66, 0.91, "45");
    cent1[2] = new TLatex(0.76, 0.91, "35");
    cent1[3] = new TLatex(0.69, 0.91, "25");
    cent1[4] = new TLatex(0.75, 0.91, "15");
    cent1[5] = new TLatex(0.80, 0.91, "7.5");
    cent1[6] = new TLatex(0.84, 0.91, "2.5");

    for(int i = 0; i < 3; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(13);
	    cent1[i]->SetTextFont(63);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

    TLine* l1[7];
    l1[0] = new TLine(404.1,0.00077, 404.1, 0.0008);
    l1[0]->SetLineWidth(2);
    l1[0]->Draw("Lsame");

    l1[1] = new TLine(717.6,0.00077, 717.6, 0.0008);
    l1[1]->SetLineWidth(2);
    l1[1]->Draw("Lsame");

    l1[2] = new TLine(1141,0.00077, 1141, 0.0008);
    l1[2]->SetLineWidth(2);
    l1[2]->Draw("Lsame");

    l1[3] = new TLine(1782,0.00077, 1782, 0.0008);
    l1[3]->SetLineWidth(2);
    //l1[3]->Draw("Lsame");

    l1[4] = new TLine(2662,0.00077, 2662, 0.0008);
    l1[4]->SetLineWidth(2);
    //l1[4]->Draw("Lsame");

    l1[5] = new TLine(3577.6,0.00077, 3577.6, 0.0008);
    l1[5]->SetLineWidth(2);
    //l1[5]->Draw("Lsame");

    l1[6] = new TLine(4474,0.00077, 4474, 0.0008);
    l1[6]->SetLineWidth(2);
    //l1[6]->Draw("Lsame");

   	TLatex* r11 = new TLatex(0.18,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base10 = makeHist("base10", "", "Centrality (%)", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{EP})#GT", 100,0,100,kBlack);

	// Remove the current axis
	base10->GetXaxis()->SetLabelOffset(999);
	base10->GetXaxis()->SetTickLength(0);

	// Redraw the new axis
	gPad->Update();
	TGaxis *newaxis1 = new TGaxis(70,
	                            -0.0009,
	                            0,
	                            -0.0009,
	                            0,
	                            70,
	                            510,"-");
	newaxis1->SetLabelOffset(-0.03);
	newaxis1->SetLabelFont(42);

	base10->GetYaxis()->SetRangeUser(-0.0009, 0.0008);
	base10->GetXaxis()->SetRangeUser(0, 70);
	base10->GetXaxis()->SetTitleColor(kBlack);

	
	fixedFontHist1D(base10,1.1,1.25);

	base10->GetYaxis()->SetTitleOffset(1.3);
	base10->GetYaxis()->SetTitleSize(base10->GetYaxis()->GetTitleSize()*1.4);
	base10->GetXaxis()->SetTitleSize(base10->GetXaxis()->GetTitleSize()*1.4);
	base10->GetYaxis()->SetLabelSize(base10->GetYaxis()->GetLabelSize()*1.4);
	base10->GetXaxis()->SetLabelSize(base10->GetXaxis()->GetLabelSize()*1.4);
	base10->GetXaxis()->SetNdivisions(8,18,0);
	base10->GetYaxis()->SetNdivisions(4,6,0);

	base10->Draw("");
	newaxis1->Draw("Asame");

	double PbPb_centralityBinCenter_fill[6];
	for(int mult = 0; mult < 6; mult++){
    	PbPb_centralityBinCenter_fill[mult] = 70 - PbPb_centralityBinCenter[mult];

	}
	
	TBox* box7[6];
	TBox* box8[6];

    double xe4;
    for(int mult = 0; mult < 6; mult++){

    	xe4 = 1;
    	double ye = total_systematics_PbPb;

    	double x1;
    	double value1;
    	gr6[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr6[1]->GetPoint(mult, x2, value2);


    	box7[mult] = new TBox(PbPb_centralityBinCenter_fill[mult]-xe4,value1[mult]-ye,PbPb_centralityBinCenter_fill[mult]+xe4,value1[mult]+ye);
		box7[mult]->SetFillColor(kRed);
        box7[mult]->SetFillColorAlpha(kGray+1,0.3);
        box7[mult]->SetFillStyle(1001);
    	box7[mult]->SetLineWidth(0);
    	box7[mult]->SetLineColor(kRed);
        box7[mult]->Draw("SAME");

		box8[mult] = new TBox(PbPb_centralityBinCenter_fill[mult]-xe4,value2[mult]-ye,PbPb_centralityBinCenter_fill[mult]+xe4,value2[mult]+ye);
		box8[mult]->SetFillColor(kBlue);
        box8[mult]->SetFillColorAlpha(kGray+1,0.3);
        box8[mult]->SetFillStyle(1001);
    	box8[mult]->SetLineWidth(0);
    	box8[mult]->SetLineColor(kBlue);
        box8[mult]->Draw("SAME");
    }


	gr5[0]->Draw("Psame");
	gr5[1]->Draw("Psame");
	gr5[2]->Draw("Psame");
	gr5[3]->Draw("Psame");

	gr6[0]->Draw("Psame");
	gr6[1]->Draw("Psame");

	TLegend *w2 = new TLegend(0.3,0.64,0.88,0.77);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(18);
    w2->SetTextFont(43);
    //w2->SetHeader("CMS");

    w2->SetNColumns(2);

    w2->AddEntry(gr6[0], "  ", "P");
    w2->AddEntry(gr6[1], "CMS PbPb #sqrt{s_{NN}} = 5.02 TeV", "P");

    w2->AddEntry(gr5[0], "  ", "P");
    w2->AddEntry(gr5[1], "ALICE PbPb #sqrt{s_{NN}} = 2.76 TeV", "P");

	w2->AddEntry(gr5[2], "  ", "P");
    w2->AddEntry(gr5[3], "STAR AuAu #sqrt{s_{NN}} = 200 GeV", "P");

    w2->Draw("same");
	
	TLatex* latex1 = new TLatex(0.29, 0.78, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.39, 0.78, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    latex2->Draw("same");

   	TLatex* r33 = new TLatex(0.18,0.84, "CMS");
    r33->SetNDC();
    r33->SetTextSize(0.04);
    r33->Draw("same");

    TLatex* r44 = new TLatex(0.27,0.84, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(24);
    r44->SetTextFont(53);
    r44->Draw("same");






}