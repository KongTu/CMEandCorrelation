#include "RiceStyle.h"
#include "function.C"

using namespace std;

double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
//double dEtaBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};

void plotQ(){


	int nBins = sizeof(dEtaBins)/ sizeof(dEtaBins[0]) - 1;

	const vector<double> dEtaBinsVector;
	for(int i = 0; i < nBins + 1; i++){

		dEtaBinsVector.push_back( dEtaBins[i] );
		cout << "dEta: " << dEtaBins[i] << endl;
	}

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_PbPb_30_40_v12..root");
	TH2D* QPlusPlus = file->Get("ana/QvsdEtaPlusPlus");
	TH2D* QMinusMinus = file->Get("ana/QvsdEtaMinusMinus");
	TH2D* QPlusMinus = file->Get("ana/QvsdEtaPlusMinus");
	//TH2D* QMinusPlus = file->Get("ana/QvsdEtaMinusPlus");
	
	TH1D* QaQb = file->Get("ana/c2_ab");
	TH1D* QaQc = file->Get("ana/c2_ac");
	TH1D* QcQb = file->Get("ana/c2_cb");
	TH1D* Ntrk = file->Get("ana/Ntrk");

	// double integral = 0;
	// for(int i = 0; i < evtWeightedQp3->GetNbinsX(); i++ ){

	// 	double temp = evtWeightedQp3->GetBinContent(i+1);
	// 	double center = evtWeightedQp3->GetBinCenter(i+1);
	// 	if(temp != 0 ){

	// 		cout << "found in bin " << i+1 <<", bin center = " << evtWeightedQp3->GetBinCenter(i+1) <<  ", with entires = " << temp << endl;
	// 		integral = integral + temp*center;
	// 	}
	// }

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double v2_b = sqrt(meanQaQb*meanQcQb/meanQaQc);
	double v2_a = sqrt(meanQaQb*meanQaQc/meanQcQb );
	double v2_ab = sqrt( meanQaQb ); 

	cout << "v2_a: " << v2_a << endl;
	cout << "v2_b: " << v2_b << endl;
	cout << "v2_ab: " << v2_ab << endl;
	
	TH1D* Qplusplus1D[48];
	TH1D* Qminusminus1D[48];
	TH1D* Qplusminus1D[48];
	TH1D* Qminusplus1D[48];

	for(int eta = 0; eta < nBins; eta++){

		Qplusplus1D[eta] = (TH1D*) QPlusPlus->ProjectionY(Form("Qplusplus1D_%d", eta), eta+1, eta+2);
		Qminusminus1D[eta] = (TH1D*) QMinusMinus->ProjectionY(Form("Qminusminus1D_%d", eta), eta+1, eta+2);
		Qplusminus1D[eta] = (TH1D*) QPlusMinus->ProjectionY(Form("Qplusminus1D_%d", eta), eta+1, eta+2);
		//Qminusplus1D[eta] = (TH1D*) QMinusPlus->ProjectionY(Form("Qminusplus1D_%d", eta), eta+1, eta+2);

	}

	TCanvas* c1 = makeCanvas("c1","");
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTicks();

	TH1D* hist = makeHist("hist","","#Delta eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", 47,0.1,4.8, kBlack);
	
	TH1D* hist1 = makeHistDifferentBins("hist1","", "","",21,dEtaBins,kBlack);
	TH1D* hist2 = makeHistDifferentBins("hist2","", "","",21,dEtaBins,kRed);
	TH1D* hist3 = makeHistDifferentBins("hist3","", "","",21,dEtaBins,kBlue);
	TH1D* hist4 = makeHistDifferentBins("hist4","", "","",21,dEtaBins,kGreen-3);

	for(int eta = 0; eta < nBins; eta++){

		hist1->SetBinContent(eta+1, Qplusplus1D[eta]->GetMean() );
		hist2->SetBinContent(eta+1, Qminusminus1D[eta]->GetMean() );
		hist3->SetBinContent(eta+1, Qplusminus1D[eta]->GetMean() );

		hist1->SetBinError(eta+1, Qplusplus1D[eta]->GetMeanError() );
		hist2->SetBinError(eta+1, Qminusminus1D[eta]->GetMeanError() );
		hist3->SetBinError(eta+1, Qplusminus1D[eta]->GetMeanError() );

		//hist4->SetBinContent(eta+1, Qminusplus1D[eta]->GetMean() );

	}

	hist->GetXaxis()->SetTitleColor(kBlack);
	hist->GetYaxis()->SetRangeUser(-0.0009,0.0009);
	hist->Draw();
	hist1->Scale(1.0/v2_ab);
	hist2->Scale(1.0/v2_ab);
	hist3->Scale(1.0/v2_ab);
	//hist4->Scale(1.0/v2_ab);
	
	hist1->Draw("Psame");
	hist2->Draw("Psame");
	hist3->Draw("Psame");
	//hist4->Draw("Psame");

	TLegend *w1 = new TLegend(0.2,0.20,0.45,0.43);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->AddEntry(hist1,"++","P");
    w1->AddEntry(hist2,"--","P");
    w1->AddEntry(hist3,"+-/-+","P");
    w1->Draw("same");
	



}