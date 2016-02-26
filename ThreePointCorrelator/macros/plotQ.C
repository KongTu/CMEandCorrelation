#include "RiceStyle.h"

using namespace std;

void plotQ(){

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v12.root");
	TH2D* QPlusPlus = file->Get("ana/QvsdEtaPlusPlus");
	TH2D* QMinusMinus = file->Get("ana/QvsdEtaMinusMinus");
	TH2D* QPlusMinus = file->Get("ana/QvsdEtaPlusMinus");
	TH2D* QMinusPlus = file->Get("ana/QvsdEtaMinusPlus");
	
	TH1D* evtWeightedQp3 = file->Get("ana/evtWeightedQp3");
	TH1D* evtWeight = file->Get("ana/evtWeight");
	TH1D* Qp3 = file->Get("ana/Qp3");
	TH1D* Ntrk = file->Get("ana/Ntrk");

	double integral = 0;
	for(int i = 0; i < evtWeightedQp3->GetNbinsX(); i++ ){

		double temp = evtWeightedQp3->GetBinContent(i+1);
		double center = evtWeightedQp3->GetBinCenter(i+1);
		if(temp != 0 ){

			cout << "found in bin " << i+1 <<", bin center = " << evtWeightedQp3->GetBinCenter(i+1) <<  ", with entires = " << temp << endl;
			integral = integral + temp*center;
		}
	}

	int integral1 = 0;
	for(int i = 0; i < evtWeight->GetNbinsX(); i++ ){

		double temp = evtWeight->GetBinContent(i+1);
		double center = evtWeight->GetBinCenter(i+1);
		if(temp != 0 ){

			cout << "found in bin " << i+1 <<", bin center = " << evtWeight->GetBinCenter(i+1) <<  ", with entires = " << temp << endl;
			integral1 = integral1 + temp*center;
		}
	}

	double integral2 = 0;
	for(int i = 0; i < Qp3->GetNbinsX(); i++ ){

		double temp = Qp3->GetBinContent(i+1);
		double center = Qp3->GetBinCenter(i+1);
		if(temp != 0 ){

			cout << "found in bin " << i+1 <<", bin center = " << Qp3->GetBinCenter(i+1) <<  ", with entires = " << temp << endl;
			integral2 = integral2 + temp*center;
		}
	}

	cout << "integral of weightedQ: " << integral << endl;
	cout << "integral of event weight: " << integral1 << endl;
	cout << "integral of Qp3: " << integral2 << endl;

	double v2_3 = sqrt(integral/integral1);
	double v2_3_Q = sqrt( ( integral2/Qp3->GetEntries() ) );
	double v2_3_1 = sqrt( ( integral/evtWeightedQp3->GetEntries() ) ); 

	cout << "v2_3: " << v2_3 << endl;
	cout << "v2_3_1: " << v2_3_1 << endl;
	cout << "v2_3_Q: " << v2_3_Q << endl;
	
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

	TCanvas* c1 = makeCanvas("c1","");
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTicks();

	TH1D* hist = makeHist("hist","","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", 48,0,4.8, kBlack);
	
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
	hist->GetYaxis()->SetRangeUser(-0.0009,0.0009);
	hist->Draw();
	hist1->Scale(1.0/v2_3_Q);
	hist2->Scale(1.0/v2_3_Q);
	hist3->Scale(1.0/v2_3_Q);
	hist4->Scale(1.0/v2_3_Q);
	
	hist1->Draw("Psame");
	hist2->Draw("Psame");
	hist3->Draw("Psame");
	hist4->Draw("Psame");

	TLegend *w1 = new TLegend(0.2,0.20,0.45,0.43);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->AddEntry(hist1,"++","P");
    w1->AddEntry(hist2,"--","P");
    w1->AddEntry(hist4,"+-/-+","P");
    w1->Draw("same");
	



}