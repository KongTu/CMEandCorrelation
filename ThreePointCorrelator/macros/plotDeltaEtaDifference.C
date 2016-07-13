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

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_5TeV_ntrkBinCenter[] = {103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};

double PbPb_ntrkCentralityBinCenter[] = {625.275, 811.118, 1025.79, 1257.64};

double PbPb_centralityBinCenter[] = {57.5, 52.5, 47.5, 42.5, 37.5, 32.5};

const int Nmults = 23;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

void plotDeltaEtaDifference(){

	gStyle->SetErrorX(0);

	TFile* file[7];
	file[0] = new TFile("../dataPoints/pPb_data.root");
	file[1] = new TFile("../dataPoints/PbPb5TeV_data.root");
	file[2] = new TFile("../dataPoints/PbPb_5TeV_centrality_data.root");
	file[3] = new TFile("../dataPoints/EPOS_data.root");
	file[4] = new TFile("../dataPoints/diff.root");
	file[5] = new TFile("../dataPoints/EPOS_PbPb_data.root");


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

	TH1D* gr5[3];
	for(int i = 0; i < 3; i++){

		gr5[i] = (TH1D*) file[4]->Get(Form("diff%d", i+1));
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

		gr1[i]->SetPoint(9,x,100);
	}

//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "[#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}]", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "[#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}]", 5000,0.1,10000,kBlack);

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
	base3->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (OS-SS)");
	base3->GetYaxis()->SetRangeUser(-0.00025,0.0016);
	base3->GetYaxis()->SetTitleOffset(1.1);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);

	TH1D* base4 = makeHist("base4", "", "#Delta#eta", "[#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}]", 48,0,4.8,kBlack);
	base4->GetYaxis()->SetRangeUser(-0.0012,0.001);
	//base4->GetYaxis()->SetRangeUser(-0.01, 0.03);
	base4->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base4->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base4,1.1,1.25);
	base4->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c} (OS-SS)");
	base4->GetYaxis()->SetRangeUser(-0.00025,0.0016);
	base4->GetYaxis()->SetTitleOffset(1.1);
	base4->GetXaxis()->SetTitleOffset(1.1);
	base4->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.2);
	base4->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.2);
	base4->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.3);
	base4->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.3);
	base4->GetXaxis()->SetNdivisions(5,9,0);
	base4->GetYaxis()->SetNdivisions(5,9,0);

    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];
    TBox *box5[50];
    TBox *box6[50];

    double xe[9];
    double xe2[9]; 
    double xe3[4];


    TLatex* r4 = new TLatex(0.18, 0.78, "#sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(24);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
	
	TLatex* latex1 = new TLatex(0.47, 0.78, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.56, 0.78, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    TLatex* latex3 = new TLatex(0.56, 0.22, "oppo - same");
    latex3->SetNDC();
    latex3->SetTextSize(25);
    latex3->SetTextFont(43);
    latex3->SetTextColor(kBlack);

    TLatex* r3 = new TLatex(0.69, 0.94, "PbPb centrality(%)");
    r3->SetNDC();
    r3->SetTextSize(21);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);

    TLatex* cent1[7];
    cent1[0] = new TLatex(0.55, 0.91, "55");
    cent1[1] = new TLatex(0.69, 0.91, "45");
    cent1[2] = new TLatex(0.80, 0.91, "35");
    cent1[3] = new TLatex(0.38, 0.91, "65");
    cent1[4] = new TLatex(0.75, 0.91, "15");
    cent1[5] = new TLatex(0.80, 0.91, "7.5");
    cent1[6] = new TLatex(0.84, 0.91, "2.5");

    for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(13);
	    cent1[i]->SetTextFont(63);
	    cent1[i]->SetTextColor(kBlack);
	    //cent1[i]->Draw("same");
    }

   	TLatex* r11 = new TLatex(0.18,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);

    TLatex* r22 = new TLatex(0.27,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);

    TCanvas* c3 = new TCanvas("c3","c3",1,1,1300,650);
    c3->Divide(2,1,0.01,0.01);
    c3->cd(2);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gPad->SetLogx(1);
	gStyle->SetOptTitle(0);

	base3->GetXaxis()->SetLabelOffset(999);
	base3->GetXaxis()->SetTickLength(0);
	
	TGaxis *newaxis2 = new TGaxis(81,
	                            -0.00025,
	                            2000,
	                            -0.00025,
	                            81,
	                            2000,
	                            510,"G");
	newaxis2->SetLabelOffset(0.0);
	newaxis2->SetLabelFont(42);
	newaxis2->SetLabelSize(0.045);

	base3->Draw();
	newaxis2->Draw("same");

	TGraphErrors* new_gr1_pPb = new TGraphErrors(9);
	TGraphErrors* new_gr2_pPb = new TGraphErrors(9);
	
	TGraphErrors* new_gr1_PbPb = new TGraphErrors(9);
	TGraphErrors* new_gr2_PbPb = new TGraphErrors(9);

	TGraphErrors* new_gr1_EPOS = new TGraphErrors(3);
	TGraphErrors* new_gr2_EPOS = new TGraphErrors(3);

	TGraphErrors* new_gr1_EPOS_PbPb = new TGraphErrors(6);

	for(int mult = 0; mult < 6; mult++){

		xe[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe[mult] = 6;
    	double ye = total_systematics_pPb;

//pPb Pb-going
		double x1, y1, y1_error;
		gr1[0]->GetPoint(mult+3, x1, y1);
		y1_error = gr1[0]->GetErrorY(mult+3);

		double x2, y2, y2_error;
		gr1[1]->GetPoint(mult+3, x2, y2);
		y2_error = gr1[1]->GetErrorY(mult+3);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_pPb->SetPoint(mult+3, x1, y2-y1);
		new_gr1_pPb->SetPointError(mult+3, 0, total_error);

		box1[mult] = new TBox(x1-xe[mult],y2-y1-ye,x1+xe[mult],y2-y1+ye);
		box1[mult]->SetFillColor(kRed);
	 	box1[mult]->SetFillColorAlpha(kGray+1,0.3);
        box1[mult]->SetFillStyle(1001);
    	box1[mult]->SetLineWidth(0);
    	box1[mult]->SetLineColor(kRed);
        box1[mult]->Draw("SAME");

//pPb p-going
		double x1, y1, y1_error;
		gr1[2]->GetPoint(mult+3, x1, y1);
		y1_error = gr1[2]->GetErrorY(mult+3);

		double x2, y2, y2_error;
		gr1[3]->GetPoint(mult+3, x2, y2);
		y2_error = gr1[3]->GetErrorY(mult+3);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		if(mult == 5) y2 = 100;
		new_gr2_pPb->SetPoint(mult+3, x1, y2-y1);
		new_gr2_pPb->SetPointError(mult+3, 0, total_error);

		box2[mult] = new TBox(x1-xe[mult],y2-y1-ye,x1+xe[mult],y2-y1+ye);
		box2[mult]->SetFillColor(kRed);
	 	box2[mult]->SetFillColorAlpha(kGray+1,0.3);
        box2[mult]->SetFillStyle(1001);
    	box2[mult]->SetLineWidth(0);
    	box2[mult]->SetLineColor(kRed);
        box2[mult]->Draw("SAME");

	}


	for(int mult = 0; mult < 9; mult++){

		xe2[mult] = 7*log(1.1*(mult+1));
    	if(mult == 0) xe2[mult] = 6;

    	double ye = total_systematics_PbPb;

//PbPb Ntrk 
		double x1, y1, y1_error;
		gr2[0]->GetPoint(mult, x1, y1);
		y1_error = gr2[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr2[1]->GetPoint(mult, x2, y2);
		y2_error = gr2[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_PbPb->SetPoint(mult, x1, y2-y1);
		new_gr1_PbPb->SetPointError(mult, 0, total_error);

		box1[mult] = new TBox(x1-xe2[mult],y2-y1-ye,x1+xe2[mult],y2-y1+ye);
		box1[mult]->SetFillColor(kRed);
	 	box1[mult]->SetFillColorAlpha(kGray+1,0.3);
        box1[mult]->SetFillStyle(1001);
    	box1[mult]->SetLineWidth(0);
    	box1[mult]->SetLineColor(kRed);
        box1[mult]->Draw("SAME");

	}

	for(int mult = 0; mult < 4; mult++){

		xe3[mult] = 30*log(1.9*(mult+1));
    	if(mult == 0) xe3[mult] = 35;
    	double ye = total_systematics_PbPb;
	
//PbPb Centrality		
		double x1, y1, y1_error;
		gr3[0]->GetPoint(mult, x1, y1);
		y1_error = gr3[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr3[1]->GetPoint(mult, x2, y2);
		y2_error = gr3[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr2_PbPb->SetPoint(mult, x1, y2-y1);
		new_gr2_PbPb->SetPointError(mult, 0, total_error);

		box2[mult] = new TBox(x1-xe3[mult],y2-y1-ye,x1+xe3[mult],y2-y1+ye);
		box2[mult]->SetFillColor(kRed);
	 	box2[mult]->SetFillColorAlpha(kGray+1,0.3);
        box2[mult]->SetFillStyle(1001);
    	box2[mult]->SetLineWidth(0);
    	box2[mult]->SetLineColor(kRed);
        box2[mult]->Draw("SAME");

	}

	for(int mult = 0; mult < 3; mult++){

//EPOS
		double x1, y1, y1_error;
		gr4[0]->GetPoint(mult, x1, y1);
		y1_error = gr4[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr4[1]->GetPoint(mult, x2, y2);
		y2_error = gr4[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);
		new_gr1_EPOS->SetPoint(mult, x1, y2-y1);
		new_gr1_EPOS->SetPointError(mult, 0, total_error);

		double x1, y1, y1_error;
		gr4[2]->GetPoint(mult, x1, y1);
		y1_error = gr4[2]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr4[3]->GetPoint(mult, x2, y2);
		y2_error = gr4[3]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr2_EPOS->SetPoint(mult, x1, y2-y1);
		new_gr2_EPOS->SetPointError(mult, 0, total_error);
	}
	
	for(int mult = 0; mult < 6; mult++){

//EPOS PbPb
		double x1, y1, y1_error;
		gr6[0]->GetPoint(mult, x1, y1);
		y1_error = gr6[0]->GetErrorY(mult);

		double x2, y2, y2_error;
		gr6[1]->GetPoint(mult, x2, y2);
		y2_error = gr6[1]->GetErrorY(mult);

		double total_error = sqrt(y1_error*y1_error + y2_error*y2_error);

		new_gr1_EPOS_PbPb->SetPoint(mult, x1, y2-y1);
		new_gr1_EPOS_PbPb->SetPointError(mult, 0, total_error);


	}
	new_gr1_pPb->SetMarkerStyle(20);
	new_gr1_pPb->SetMarkerSize(1.4);
	new_gr1_pPb->SetMarkerColor(kRed);
	new_gr1_pPb->SetLineColor(kRed);
	new_gr1_pPb->Draw("Psame");

	new_gr2_pPb->SetMarkerStyle(21);
	new_gr2_pPb->SetMarkerSize(1.4);
	new_gr2_pPb->SetMarkerColor(kBlue);
	new_gr2_pPb->SetLineColor(kBlue);
	new_gr2_pPb->Draw("Psame");

	new_gr1_PbPb->SetMarkerStyle(28);
	new_gr1_PbPb->SetMarkerSize(1.4);
	new_gr1_PbPb->SetMarkerColor(kBlack);
	new_gr1_PbPb->SetLineColor(kBlack);
	new_gr1_PbPb->Draw("Psame");

	new_gr2_PbPb->SetMarkerStyle(28);
	new_gr2_PbPb->SetMarkerSize(1.4);
	new_gr2_PbPb->SetMarkerColor(kBlack);
	new_gr2_PbPb->SetLineColor(kBlack);
	new_gr2_PbPb->Draw("Psame");

	new_gr1_EPOS->SetMarkerStyle(28);
	new_gr1_EPOS->SetMarkerSize(1.4);
	new_gr1_EPOS->SetMarkerColor(kGreen+2);
	new_gr1_EPOS->SetLineColor(kGreen+2);
	new_gr1_EPOS->SetLineWidth(3.0);
	new_gr1_EPOS->SetLineStyle(5.0);
	new_gr1_EPOS->Draw("same");

	new_gr2_EPOS->SetMarkerStyle(28);
	new_gr2_EPOS->SetMarkerSize(1.4);
	new_gr2_EPOS->SetMarkerColor(kBlack);
	new_gr2_EPOS->SetLineColor(kBlack);
	new_gr2_EPOS->SetLineWidth(3.0);
	new_gr2_EPOS->SetLineStyle(5.0);
	new_gr2_EPOS->Draw("same");
	
	new_gr1_EPOS_PbPb->SetMarkerSize(1.4);
	new_gr1_EPOS_PbPb->SetMarkerColor(kBlue);
	new_gr1_EPOS_PbPb->SetLineColor(kBlue);
	new_gr1_EPOS_PbPb->SetLineWidth(3.0);
	new_gr1_EPOS_PbPb->SetLineStyle(5.0);
	new_gr1_EPOS_PbPb->Draw("same");

	TLegend *w3 = new TLegend(0.50,0.63,0.7,0.8);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    //w3->AddEntry(new_gr1_pPb, "pPb, #Psi_{EP}(Pb-going)", "P");
    //w3->AddEntry(new_gr2_pPb, "pPb, #Psi_{EP}(p-going)", "P");
    w3->AddEntry(new_gr1_EPOS, "pPb EPOS, #Psi_{EP}(Pb-going)", "L" );    
    w3->AddEntry(new_gr2_EPOS, "pPb EPOS, #Psi_{EP}(p-going)", "L" );    
    //w3->AddEntry(new_gr1_PbPb, "PbPb", "P");
    w3->AddEntry(new_gr1_EPOS_PbPb, "PbPb EPOS", "L" );
    w3->Draw("same");

    r4->Draw("same");
    r11->Draw("same");
    r22->Draw("same");
    r3->Draw("same");
    for(int i = 0; i < 4; i++){
    	cent1[i]->SetNDC();
    	cent1[i]->SetTextSize(13);
	    cent1[i]->SetTextFont(63);
	    cent1[i]->SetTextColor(kBlack);
	    cent1[i]->Draw("same");
    }

    TLine* l2[7];
    l2[0] = new TLine(404.1,0.00157, 404.1, 0.0016);
    l2[0]->SetLineWidth(2);
    l2[0]->Draw("Lsame");

    l2[1] = new TLine(717.6,0.00157, 717.6, 0.0016);
    l2[1]->SetLineWidth(2);
    l2[1]->Draw("Lsame");

    l2[2] = new TLine(1141,0.00157, 1141, 0.0016);
    l2[2]->SetLineWidth(2);
    l2[2]->Draw("Lsame");

    l2[3] = new TLine(81.412,0.00157, 81.412, 0.0016);
    l2[3]->SetLineWidth(2);
    //l2[3]->Draw("Lsame");

    l2[4] = new TLine(197,0.00157, 197, 0.0016);
    l2[4]->SetLineWidth(2);
    l2[4]->Draw("Lsame");

//	latex3->Draw("same");

    c3->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.14);
	gPad->SetBottomMargin(0.13);
	gPad->SetRightMargin(0.05);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
    base4->Draw();


    TLatex* r33 = new TLatex(0.25, 0.82, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    //r33->Draw("same");

    TLatex* lmult = new TLatex(0.20, 0.84, "185 #leq N^{offline}_{trk} < 220");
    lmult->SetNDC();
    lmult->SetTextSize(24);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

    double value5[50];
    double value5_error[50];
    double value6[50];
    double value6_error[50];
    double value7[50];
    double value7_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value5[deta] = gr5[0]->GetBinContent(deta+1);
    	value5_error[deta] = gr5[0]->GetBinError(deta+1);

    	value6[deta] = gr5[1]->GetBinContent(deta+1);
    	value6_error[deta] = gr5[1]->GetBinError(deta+1);

    	value7[deta] = gr5[2]->GetBinContent(deta+1);
    	value7_error[deta] = gr5[2]->GetBinError(deta+1);

    }

    TBox *box5[50];
    TBox *box6[50];
    TBox *box7[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe4 = 0.065;
    	double ye = total_systematics_PbPb;

    	box5[deta] = new TBox(dEtaReBinCenter2[deta]-xe4,value5[deta]-ye,dEtaReBinCenter2[deta]+xe4,value5[deta]+ye);
		box5[deta]->SetFillColor(kRed);
        box5[deta]->SetFillColorAlpha(kGray+1,0.4);
        box5[deta]->SetFillStyle(1001);
    	box5[deta]->SetLineWidth(0);
    	box5[deta]->SetLineColor(kRed);
        box5[deta]->Draw("SAME");

		box6[deta] = new TBox(dEtaReBinCenter2[deta]-xe4,value6[deta]-ye,dEtaReBinCenter2[deta]+xe4,value6[deta]+ye);
		box6[deta]->SetFillColor(kBlue);
        box6[deta]->SetFillColorAlpha(kGray+1,0.4);
        box6[deta]->SetFillStyle(1001);
    	box6[deta]->SetLineWidth(0);
    	box6[deta]->SetLineColor(kBlue);
        box6[deta]->Draw("SAME");

        box7[deta] = new TBox(dEtaReBinCenter2[deta]-xe4,value7[deta]-ye,dEtaReBinCenter2[deta]+xe4,value7[deta]+ye);
		box7[deta]->SetFillColor(kBlue);
        box7[deta]->SetFillColorAlpha(kGray+1,0.4);
        box7[deta]->SetFillStyle(1001);
    	box7[deta]->SetLineWidth(0);
    	box7[deta]->SetLineColor(kBlue);
        box7[deta]->Draw("SAME");
    }

    gr5[0]->Draw("Psame");
    gr5[1]->Draw("Psame");
    gr5[2]->Draw("Psame");

   	TLatex* r11 = new TLatex(0.67,0.91, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.77,0.91, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

    TLatex* r44 = new TLatex(0.20, 0.78, "#sqrt{s_{NN}} = 5.02 TeV");
    r44->SetNDC();
    r44->SetTextSize(23);
    r44->SetTextFont(43);
    r44->SetTextColor(kBlack);
    r44->Draw("same");

	TLegend *w4 = new TLegend(0.5,0.63,0.7,0.8);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(20);
    w4->SetTextFont(43);
    w4->AddEntry(gr5[0], "pPb, #Psi_{EP}(Pb-going)", "P");
    w4->AddEntry(gr5[1], "pPb, #Psi_{EP}(p-going)", "P");
    w4->AddEntry(gr5[2], "PbPb", "P");
    w4->Draw("same");


//	latex3->Draw("same");

}