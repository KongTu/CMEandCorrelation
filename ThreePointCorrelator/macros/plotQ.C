#include "RiceStyle.h"

using namespace std;

void plotQ(){

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v2.root");
	TH2D* QPlusPlus = file->Get("ana/QvsdEtaPlusPlus");
	TH2D* QMinusMinus = file->Get("ana/QvsdEtaMinusMinus");
	TH2D* QPlusMinus = file->Get("ana/QvsdEtaPlusMinus");
	TH2D* QMinusPlus = file->Get("ana/QvsdEtaMinusPlus");



	TH1D* Qplusplus1D[48];
	TH1D* Qminusminus1D[48];
	TH1D* Qplusminus1D[48];
	TH1D* Qminusplus1D[48];


	for(int eta = 0; eta < 48; eta++){

		Qplusplus1D[eta] = (TH1D*) QPlusPlus->ProjectionY(Form("Qplusplus1D_%d", eta), eta+1, eta+2);
		Qminusminus1D[eta] = (TH1D*) QMinusMinus->ProjectionY(Form("Qminusminus1D_%d", eta), eta+1, eta+2);
		Qplusminus1D[eta] = (TH1D*) QPlusMinus->ProjectionY(Form("Qplusminus1D_%d", eta), eta+1, eta+2);
		Qminusplus1D[eta] = (TH1D*) QMinusPlus->ProjectionY(Form("Qminusplus1D_%d", eta), eta+1, eta+2);

	}

	TH1D* hist = makeHist("hist","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kBlack);
	TH1D* hist1 = makeHist("hist1","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kBlack);
	TH1D* hist2 = makeHist("hist2","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kRed);
	TH1D* hist3 = makeHist("hist3","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kBlue);
	TH1D* hist4 = makeHist("hist4","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>", 48,0,4.8, kGreen-3);

	for(int eta = 0; eta < 48; eta++){

		hist1->SetBinContent(eta+1, Qplusplus1D[eta]->GetMean() );
		hist2->SetBinContent(eta+1, Qminusminus1D[eta]->GetMean() );
		hist3->SetBinContent(eta+1, Qplusminus1D[eta]->GetMean() );
		hist4->SetBinContent(eta+1, Qminusplus1D[eta]->GetMean() );

	}

	hist->GetXaxis()->SetTitleColor(kBlack);
	hist->GetYaxis()->SetRangeUser(-0.00005,0.00005);
	hist->Draw();
	hist1->Draw("Psame");
	hist2->Draw("Psame");
	hist3->Draw("Psame");
	hist4->Draw("Psame");
	



}