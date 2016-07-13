#include "RiceStyle.h"

using namespace std;

void plotQ2vsV2(){
	
	TGaxis::SetMaxDigits(3);


	double Q2[] = {1,2,3,4,5,6,7,8,9,10,11};
	double pPb_v2_bincenter[] = {0.06167, 0.063072, 0.0648051, 0.0665957, 0.0682464, 0.0699678, 0.0723474, 0.07504, 0.0784839, 0.0821738, 0.0860937};

	TGraphErrors * gr1 = new TGraphErrors(11, Q2, pPb_v2_bincenter, 0,0);
	

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.18);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base = makeHist("base","base", "HF Q_{2} class", "v_{2,tracker}", 1000,0,12,kBlue);
	base->GetYaxis()->SetRangeUser(0.055,0.09);
	base->GetYaxis()->SetTitleSize(0.06.);
	base->Draw();


	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(1.4);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);
	gr1->Draw("PLsame");

	TLatex* r4 = new TLatex(0.22, 0.84, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);
    r4->Draw("same");

    TLatex* lmult = new TLatex(0.22, 0.78, "185 #leq N^{offline}_{trk} < 260");
    lmult->SetNDC();
    lmult->SetTextSize(24);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

}

