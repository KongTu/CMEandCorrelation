#include "RiceStyle.h"
#include "function.h"

using namespace std;

void plotCorrectedQ(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v12.root");
	
	TH2D* QPlusPlus = file->Get("ana/QvsdEtaPlusPlus");
	TH2D* QMinusMinus = file->Get("ana/QvsdEtaMinusMinus");
	TH2D* QPlusMinus = file->Get("ana/QvsdEtaPlusMinus");
	TH2D* QMinusPlus = file->Get("ana/QvsdEtaMinusPlus");
	
	TH1D* evtWeightedQp3 = file->Get("ana/evtWeightedQp3");
	TH1D* evtWeight = file->Get("ana/evtWeight");
	TH1D* Qp3 = file->Get("ana/Qp3");
	TH1D* Ntrk = file->Get("ana/Ntrk");

	TH1D* HFsinSum = file->Get("ana/HFsinSum");
	TH1D* HFcosSum = file->Get("ana/HFcosSum");
	TH1D* weightSum = file->Get("ana/weightSum");

//correction from HF:

	double cosHFSum = getWeightedSum( HFcosSum, weightSum );
	double sinHFSum = getWeightedSum( HFsinSum, weightSum );

//correction table from TRK:

	ifstream ifile("v12_table.txt");
	if(ifile.is_open() ){

		vector< vector<double>> res;
		for(int line = 0; line < Nbins; line++){

			vector<double> t;
			double v1, v2, v3, v4;

			ifile >> v1 >> v2 >> v3 >> v4;
			t.push_back(v1); t.push_back(v2);t.push_back(v3);t.push_back(v4);

			res.push_back(t);
		}
	}

	for(int ieta = 0; ieta < Nbins; ieta++){
		for(int jeta = 0; jeta < Nbins; jeta++){

			if( ieta == jeta ) continue;

			double deltaEta = fabs( etabins[ieta] - etabins[jeta] );

			double cosPlusSum_a = res[ieta][0];
			double cosPlusSum_b = res[jeta][0];
			double sinPlusSum_a = res[ieta][1];
			double sinPlusSum_b = res[jeta][1];

			double cosMinusSum_a = res[ieta][2];
			double cosMinusSum_b = res[jeta][2];
			double sinMinusSum_a = res[ieta][3];
			double sinMinusSum_b = res[jeta][3];

			double PlusPlus = getReal(cosPlusSum_a, cosPlusSum_b, cosHFSum, sinPlusSum_a, sinPlusSum_b, sinHFSum);
			double MinusMinus = getReal(cosMinusSum_a, cosMinusSum_b, cosHFSum, sinMinusSum_a, sinMinusSum_b, sinHFSum);
			double PlusMinus = getReal(cosPlusSum_a, cosMinusSum_b, cosHFSum, sinPlusSum_a, sinMinusSum_b, sinHFSum);

			corrVsEtaPlusPlus->Fill( deltaEta, PlusPlus );
			corrVsEtaMinusMinus->Fill( deltaEta, MinusMinus );
			corrVsEtaPlusMinus->Fill( deltaEta, PlusMinus );
		}
	}

//projectionY from 2D
	for(int eta = 0; eta < 48; eta++){

		Qplusplus1D[eta] = (TH1D*) QPlusPlus->ProjectionY(Form("Qplusplus1D_%d", eta), eta+1, eta+2);
		Qminusminus1D[eta] = (TH1D*) QMinusMinus->ProjectionY(Form("Qminusminus1D_%d", eta), eta+1, eta+2);
		Qplusminus1D[eta] = (TH1D*) QPlusMinus->ProjectionY(Form("Qplusminus1D_%d", eta), eta+1, eta+2);
		Qminusplus1D[eta] = (TH1D*) QMinusPlus->ProjectionY(Form("Qminusplus1D_%d", eta), eta+1, eta+2);

		//correction
		corrVsEtaPlusPlus1D[eta] = (TH1D*)corrVsEtaPlusPlus->ProjectionY(Form("corrVsEtaPlusPlus1D_%d", eta), eta+1, eta+2);
		corrVsEtaMinusMinus1D[eta] = (TH1D*)corrVsEtaMinusMinus->ProjectionY(Form("corrVsEtaMinusMinus1D_%d", eta), eta+1, eta+2);
		corrVsEtaPlusMinus1D[eta] = (TH1D*)corrVsEtaPlusMinus->ProjectionY(Form("corrVsEtaPlusMinus1D_%d", eta), eta+1, eta+2);
		
	}

//fill histogram with averages
	for(int eta = 0; eta < 48; eta++){

		corrHistPlusPlus->SetBinContent(eta+1, corrVsEtaPlusPlus1D[eta]->GetMean() );
		corrHistMinusMinus->SetBinContent(eta+1, corrVsEtaMinusMinus1D[eta]->GetMean() );
		corrHistPlusMinus->SetBinContent(eta+1, corrVsEtaPlusMinus1D[eta]->GetMean() );

		double corrFactorPlusPlus = corrVsEtaPlusPlus1D[eta]->GetMean();
		double corrFactorMinusMinus = corrVsEtaMinusMinus1D[eta]->GetMean();
		double corrFactorPlusMinus = corrVsEtaPlusMinus1D[eta]->GetMean();

		hist1->SetBinContent(eta+1, Qplusplus1D[eta]->GetMean() + corrFactorPlusPlus );
		hist2->SetBinContent(eta+1, Qminusminus1D[eta]->GetMean() + corrFactorMinusMinus );
		hist3->SetBinContent(eta+1, Qplusminus1D[eta]->GetMean() + corrFactorPlusMinus );
		hist4->SetBinContent(eta+1, Qminusplus1D[eta]->GetMean() + corrFactorPlusMinus );

		
	}

	TCanvas* c2 = makeCanvas("c2","");
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTicks();
	
	corrHistPlusPlus->GetXaxis()->SetTitleColor(kBlack);
	corrHistPlusPlus->Draw("");
	corrHistMinusMinus->Draw("same");
	corrHistPlusMinus->Draw("same");

	TLegend *w2 = new TLegend(0.5,0.60,0.75,0.8);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->AddEntry(corrHistPlusPlus,"++","P");
    w2->AddEntry(corrHistMinusMinus,"--","P");
    w2->AddEntry(corrHistPlusMinus,"+-/-+","P");
    w2->Draw("same");

// plotting:
	TCanvas* c1 = makeCanvas("c1","");
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTicks();

	hist->GetXaxis()->SetTitleColor(kBlack);
	hist->GetYaxis()->SetRangeUser(-0.0009,0.0009);
	hist->Draw();
	double v2_3_Q = getV2( Qp3 );
	hist1->Scale(1.0/v2_3_Q);
	hist2->Scale(1.0/v2_3_Q);
	hist3->Scale(1.0/v2_3_Q);
	hist4->Scale(1.0/v2_3_Q);
	
	hist1->Draw("Psame");
	hist2->Draw("Psame");
	hist3->Draw("Psame");
	hist4->Draw("Psame");

	TLegend *w1 = new TLegend(0.5,0.20,0.75,0.43);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->AddEntry(hist1,"++","P");
    w1->AddEntry(hist2,"--","P");
    w1->AddEntry(hist4,"+-/-+","P");
    w1->Draw("same");
	



}