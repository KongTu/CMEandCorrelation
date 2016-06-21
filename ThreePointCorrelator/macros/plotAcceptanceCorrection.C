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


void plotAcceptanceCorrection(){

	TFile* file[2];
	file[0] = new TFile("../dataPoints/referencePoint.root");
	file[1] = new TFile("../dataPoints/acceptancePoint.root");

	TH1D* temp1[2];
	TH1D* temp_11[2];
	TH1D* temp3[2];
	TH1D* temp_33[2];

	TH1D* temp2[2];
	TH1D* temp_22[2];
	TH1D* temp4[2];
	TH1D* temp_44[2];

	for(int i = 0; i < 2; i++){

		temp1[i] = (TH1D*) file[i]->Get("temp1");
		temp_11[i] = (TH1D*) file[i]->Get("temp_11");
		temp3[i] = (TH1D*) file[i]->Get("temp3");
		temp_33[i] = (TH1D*) file[i]->Get("temp_33");

		temp2[i] = (TH1D*) file[i]->Get("temp2");
		temp_22[i] = (TH1D*) file[i]->Get("temp_22");
		temp4[i] = (TH1D*) file[i]->Get("temp4");
		temp_44[i] = (TH1D*) file[i]->Get("temp_44");

	}

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "same", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,4.8,kBlack);
	TH1D* base2 = makeHist("base2", "oppo", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT/v_{2,c}", 48,0,4.8,kBlack);

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
	base1->GetXaxis()->SetTitleOffset(0.9);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.9);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);

	TCanvas* c1 = makeMultiCanvas("c1","c1",2,1);
	c1->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

	temp1[0]->SetMarkerColor(kRed);
	temp1[1]->SetMarkerColor(kBlack);
	temp2[1]->SetLineColor(kBlue);
	temp2[1]->SetMarkerColor(kBlue);

	temp1[0]->Draw("Psame");
	temp1[1]->Draw("Psame");
	temp2[1]->Draw("Psame");

	c1->cd(2);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base2->Draw();

	temp_11[0]->SetMarkerColor(kRed);
	temp_11[0]->SetMarkerStyle(20);
	temp_11[1]->SetMarkerColor(kBlack);
	temp_11[1]->SetMarkerStyle(20);
	temp_22[1]->SetMarkerColor(kBlue);
	temp_22[1]->SetMarkerStyle(24);

	temp_11[0]->Draw("Psame");
	temp_11[1]->Draw("Psame");
	temp_22[1]->Draw("Psame");


	TLegend *w2 = new TLegend(0.45,0.2,0.6,0.4);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(45);
    w2->AddEntry(temp1[0], " True ", "P");
    w2->AddEntry(temp1[1], " After correction ", "P");
    w2->AddEntry(temp2[1], " Before correction ", "P");
    w2->Draw("same");

	// c1->cd(3);
	// gPad->SetLeftMargin(0.20);
	// gPad->SetBottomMargin(0.13);
	// gPad->SetTopMargin(0.06);
	// gPad->SetTicks();
	// base2->Draw();

	// temp3[0]->SetMarkerColor(kRed);
	// temp3[1]->SetMarkerColor(kBlack);
	// temp4[0]->SetMarkerStyle(25);
	// temp4[0]->SetMarkerColor(kBlue);

	// temp3[0]->Draw("Psame");
	// temp3[1]->Draw("Psame");
	// temp4[1]->Draw("Psame");

	// c1->cd(4);
	// gPad->SetLeftMargin(0.20);
	// gPad->SetBottomMargin(0.13);
	// gPad->SetTopMargin(0.06);
	// gPad->SetTicks();
	// base2->Draw();

	// temp_33[0]->SetMarkerColor(kRed);
	// temp_33[1]->SetMarkerStyle(20);
	// temp_33[1]->SetMarkerColor(kBlack);
	// temp_44[0]->SetMarkerColor(kBlue);

	// temp_33[0]->Draw("Psame");
	// temp_44[1]->Draw("Psame");
	// temp_33[1]->Draw("Psame");


}