#include "RiceStyle.h"

using namespace std;

void plotSingleParticleClosure(){
	
	TFile* file1 = new TFile("../rootfiles/test.root");
	TFile* file2 = new TFile("../rootfiles/test_gen.root");

	TH1D* reco_pt = (TH1D*) file1->Get("ana/trkPt");
	TH1D* reco_eta = (TH1D*) file1->Get("ana/trk_eta");

	TH1D* gen_pt = (TH1D*) file2->Get("ana/trkPt");
	TH1D* gen_eta = (TH1D*) file2->Get("ana/trk_eta");

	TCanvas* c1 = makeMultiCanvas("c1","c1",2,1);
	c1->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.23);
	gPad->SetBottomMargin(0.16);

	reco_pt->SetStats(kFALSE);
	reco_pt->SetTitle("RECO_weighted/GEN");
	reco_pt->GetYaxis()->SetRangeUser(0.8,1.2);
	reco_pt->GetXaxis()->SetRangeUser(0.25,3.2);
	reco_pt->Divide( gen_pt );
	reco_pt->Draw();

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.23);
	gPad->SetBottomMargin(0.16);
	
	reco_eta->SetTitle("RECO_weighted/GEN");
	reco_eta->SetStats(kFALSE);
	reco_eta->GetYaxis()->SetRangeUser(0.8,1.2);
	reco_eta->Divide( gen_eta );
	reco_eta->Draw();


}