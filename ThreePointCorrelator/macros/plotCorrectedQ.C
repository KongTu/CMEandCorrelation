#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

/*****run instructions:
1. run makeCorrectionTable.C for correction tables
2. load correctionTable txt file into this macros
3. run this macros for corrected results
*/

vector<TH2D*> loading2DHistogram( TFile * file, TString hName, int Nbins ){

	vector<TH2D *> spectra;

  	stringstream HistName;

  	for (int mult = 0; mult < Nbins; mult++){

	    HistName.str("");

	    HistName << hName;
	    HistName << mult;

	    TH2D* temp2 = (TH2D*)file->Get( HistName.str().c_str() );
	    temp2->SetMarkerSize(1.3);
		temp2->SetMarkerStyle(20);
		temp2->SetStats(kFALSE);
		temp2->SetTitle("");
		temp2->SetXTitle("");
		temp2->SetYTitle("");
	    
	    spectra.push_back( temp2 );
  	}

  	return spectra;

}

void plotCorrectedQ(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v20.root");

	vector<TH2D*> QPlusPlus = loading2DHistogram(file, "ana/QvsdEtaPlusPlus_", 3);
	vector<TH2D*> QMinusMinus = loading2DHistogram(file, "ana/QvsdEtaMinusMinus_", 3);
	vector<TH2D*> QPlusMinus = loading2DHistogram(file, "ana/QvsdEtaPlusMinus_", 3);
	
	vector<TH1D*> averageCosHF = loadingHistogram(file, "ana/averageCosHF_", 2);
	vector<TH1D*> averageSinHF = loadingHistogram(file, "ana/averageSinHF_", 2);

	TH1D* QaQb = file->Get("ana/c2_ab");
	TH1D* QaQc = file->Get("ana/c2_ac");
	TH1D* QcQb = file->Get("ana/c2_cb");
	TH1D* Ntrk = file->Get("ana/Ntrk");

//correction table from TRK:

	ifstream ifile("PbPb_v10.txt");
	if(ifile.is_open() ){

		vector< vector<double>> res;
		for(int line = 0; line < Nbins; line++){

			vector<double> t;
			double v1, v2, v3, v4, v5, v6, v7, v8;

			ifile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8;
			t.push_back(v1); t.push_back(v2);t.push_back(v3);t.push_back(v4);
			t.push_back(v5); t.push_back(v6);t.push_back(v7);t.push_back(v8);

			res.push_back(t);
		}
	}

//all HF correction factors are repeated in each line with different eta bins of TRK
//only one line is needed. 

	double HFplusCos = res[0][4];
	double HFplusSin = res[0][5];
	double HFminusCos = res[0][6];
	double HFminusSin = res[0][7];

	for(int type = 0; type < 3; type++){
		
		corrVsEtaPlusPlus[type] = new TH2D(Form("corrVsEtaPlusPlus_%d", type), "corrVsEtaPlusPlus", NdEtaBins, dEtaBins, 100000, -0.00001, 0.00001);
		corrVsEtaMinusMinus[type] = new TH2D(Form("corrVsEtaMinusMinus_%d", type), "corrVsEtaMinusMinus", NdEtaBins, dEtaBins, 100000, -0.00001, 0.00001);
		corrVsEtaPlusMinus[type] = new TH2D(Form("corrVsEtaPlusMinus_%d", type), "corrVsEtaPlusMinus", NdEtaBins, dEtaBins, 100000, -0.00001, 0.00001);
	}

	for(int ieta = 0; ieta < Nbins; ieta++){
		for(int jeta = 0; jeta < Nbins; jeta++){

			if( ieta == jeta ) continue;

			double deltaEta = fabs( etabins[ieta] - etabins[jeta] );
			
			double TRKcosPlus_a = res[ieta][0];
			double TRKsinPlus_a = res[ieta][1];
			double TRKcosMinus_a = res[ieta][2];
			double TRKsinMinus_a = res[ieta][3];
			
			double TRKcosPlus_b = res[jeta][0];
			double TRKsinPlus_b = res[jeta][1];
			double TRKcosMinus_b = res[jeta][2];
			double TRKsinMinus_b = res[jeta][3];

			for(int type = 0; type < 3; type++){
				
				double HFCos = 0.;
				double HFSin = 0.;

				if( type == 0 ){//HF plus side

					HFCos = res[0][4];
					HFSin = res[0][5];
				}
				else if( type == 1 ){// HF mimus side
					
					HFCos = res[0][6];
					HFSin = res[0][7];
				}
				else if( type == 2 ){// HF both 

					HFCos = (res[0][4] + res[0][6]) / 2.0;
					HFSin = (res[0][5] + res[0][7]) / 2.0;
				}
				else{
					return;
				}
				
				double PlusPlus = getReal(TRKcosPlus_a, TRKcosPlus_b, HFCos, TRKsinPlus_a, TRKsinPlus_b, HFSin);
				double MinusMinus = getReal(TRKcosMinus_a, TRKcosMinus_b, HFCos, TRKsinMinus_a, TRKsinMinus_b, HFSin);
				double PlusMinus = getReal(TRKcosPlus_a, TRKcosMinus_b, HFCos, TRKsinPlus_a, TRKsinMinus_b, HFSin);

				corrVsEtaPlusPlus[type]->Fill( deltaEta, PlusPlus );
				corrVsEtaMinusMinus[type]->Fill( deltaEta, MinusMinus );
				corrVsEtaPlusMinus[type]->Fill( deltaEta, PlusMinus );
			}
		}
	}

//projectionY from 2D

	for( int type = 0; type < 3; type++ ){

		corrHistPlusPlus[type] = makeHistDifferentBins(Form("corrHistPlusPlus_%d", type),"","#Delta#eta", "correction factor", NdEtaBins, dEtaBins, kBlack);
		corrHistMinusMinus[type] = makeHistDifferentBins(Form("corrHistMinusMinus_%d", type),"","#Delta#eta", "correction factor", NdEtaBins, dEtaBins, kRed);
		corrHistPlusMinus[type] = makeHistDifferentBins(Form("corrHistPlusMinus_%d", type),"","#Delta#eta", "correction factor", NdEtaBins, dEtaBins, kBlue);

		show1[type] = makeHistDifferentBins(Form("show1_%d",type),"","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", NdEtaBins, dEtaBins, kBlack);
		show2[type] = makeHistDifferentBins(Form("show2_%d",type),"","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", NdEtaBins, dEtaBins, kRed);
		show3[type] = makeHistDifferentBins(Form("show3_%d",type),"","#Delta#eta", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", NdEtaBins, dEtaBins, kBlue);
		
		for(int eta = 0; eta < NdEtaBins; eta++){

			Qplusplus1D[type][eta] = (TH1D*) QPlusPlus[type]->ProjectionY(Form("Qplusplus1D_%d_%d", type, eta), eta+1, eta+1);
			Qminusminus1D[type][eta] = (TH1D*) QMinusMinus[type]->ProjectionY(Form("Qminusminus1D_%d_%d", type, eta), eta+1, eta+1);
			Qplusminus1D[type][eta] = (TH1D*) QPlusMinus[type]->ProjectionY(Form("Qplusminus1D_%d_%d", type, eta), eta+1, eta+1);

			//correction
			corrVsEtaPlusPlus1D[type][eta] = (TH1D*)corrVsEtaPlusPlus[type]->ProjectionY(Form("corrVsEtaPlusPlus1D_%d_%d", type, eta), eta+1, eta+1);
			corrVsEtaMinusMinus1D[type][eta] = (TH1D*)corrVsEtaMinusMinus[type]->ProjectionY(Form("corrVsEtaMinusMinus1D_%d_%d", type, eta), eta+1, eta+1);
			corrVsEtaPlusMinus1D[type][eta] = (TH1D*)corrVsEtaPlusMinus[type]->ProjectionY(Form("corrVsEtaPlusMinus1D_%d_%d", type, eta), eta+1, eta+1);
			
			double corrFactorPlusPlus = corrVsEtaPlusPlus1D[type][eta]->GetMean();
			double corrFactorMinusMinus = corrVsEtaMinusMinus1D[type][eta]->GetMean();
			double corrFactorPlusMinus = corrVsEtaPlusMinus1D[type][eta]->GetMean();

			double corrFactorPlusPlusErr = corrVsEtaPlusPlus1D[type][eta]->GetMeanError();
			double corrFactorMinusMinusErr = corrVsEtaMinusMinus1D[type][eta]->GetMeanError();
			double corrFactorPlusMinusErr = corrVsEtaPlusMinus1D[type][eta]->GetMeanError();		

			corrHistPlusPlus[type]->SetBinContent(eta+1, corrFactorPlusPlus );
			corrHistPlusPlus[type]->SetBinError(eta+1, corrFactorPlusPlusErr );

			corrHistMinusMinus[type]->SetBinContent(eta+1, corrFactorMinusMinus );
			corrHistMinusMinus[type]->SetBinError(eta+1, corrFactorMinusMinusErr );

			corrHistPlusMinus[type]->SetBinContent(eta+1, corrFactorPlusMinus );
			corrHistPlusMinus[type]->SetBinError(eta+1, corrFactorPlusMinusErr );

			double Qplusplus1DMean = Qplusplus1D[type][eta]->GetMean();


			cout << "Qplusplus1DMean " << eta+1 << " has average: " << Qplusplus1DMean << endl;

			double Qminusminus1DMean = Qminusminus1D[type][eta]->GetMean();
			double Qplusminus1DMean = Qplusminus1D[type][eta]->GetMean();

			double Qplusplus1DMeanErr = Qplusplus1D[type][eta]->GetMeanError();
			double Qminusminus1DMeanErr = Qminusminus1D[type][eta]->GetMeanError();
			double Qplusminus1DMeanErr = Qplusminus1D[type][eta]->GetMeanError();

			show1[type]->SetBinContent(eta+1, Qplusplus1DMean + corrFactorPlusPlus );
			show1[type]->SetBinError(eta+1, sqrt(Qplusplus1DMeanErr*Qplusplus1DMeanErr + corrFactorPlusPlusErr*corrFactorPlusPlusErr) );
			
			show2[type]->SetBinContent(eta+1, Qminusminus1DMean + corrFactorMinusMinus );
			show2[type]->SetBinError(eta+1, sqrt(Qminusminus1DMeanErr*Qminusminus1DMeanErr + corrFactorMinusMinusErr*corrFactorMinusMinusErr) );

			show3[type]->SetBinContent(eta+1, Qplusminus1DMean + corrFactorPlusMinus );
			show3[type]->SetBinError(eta+1, sqrt(Qplusminus1DMeanErr*Qplusminus1DMeanErr + corrFactorPlusMinusErr*corrFactorPlusMinusErr) );
		}
	}

	string txt[3] = {"HF+", "HF-", "HF+/-"};

	//we say a is HF-, b is HF+, ab is combined

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;

	double bCorr = (averageCosHF[0]->GetMean() * averageCosHF[0]->GetMean()) +  ( averageSinHF[0]->GetMean() * averageSinHF[0]->GetMean() );
	double aCorr = (averageCosHF[1]->GetMean() * averageCosHF[1]->GetMean()) +  ( averageSinHF[1]->GetMean() * averageSinHF[1]->GetMean() );
	
	double m1 = (averageCosHF[0]->GetMean() + averageCosHF[1]->GetMean())/2.0;
	double m2 = (averageSinHF[0]->GetMean() + averageSinHF[1]->GetMean())/2.0;
	double abCorr = m1*m1 + m2*m2;

//correction factor for the V2 values from HF. 
	double v2[3];
	v2[0] = sqrt(c2_b - bCorr);
	v2[1] = sqrt(c2_a - aCorr );
	v2[2] = sqrt(c2_ab - abCorr );

	cout << "v2[2]: " << v2[2] << endl; 

	TCanvas* c1[3];
	TCanvas* c2[3];

	for(int type = 0; type < 3; type++ ){

		c1[type] = makeCanvas(Form("c1_%d", type), txt[type].c_str());
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetTicks();

		hist->GetXaxis()->SetTitleColor(kBlack);
		hist->GetYaxis()->SetRangeUser(-0.009,0.009);
		hist->Draw();

		show1[type]->Scale( 1.0/v2[type] );
		show2[type]->Scale( 1.0/v2[type] );
		show3[type]->Scale( 1.0/v2[type] );

		show1[type]->Draw("Psame");
		show2[type]->Draw("Psame");
		show3[type]->Draw("Psame");

		TLegend *w1 = new TLegend(0.5,0.20,0.75,0.43);
		w1->SetLineColor(kWhite);
		w1->SetFillColor(0);
		w1->AddEntry(show1[type],"++","P");
		w1->AddEntry(show2[type],"--","P");
		w1->AddEntry(show3[type],"+-/-+","P");
		w1->Draw("same");
	}

	for(int type = 0; type < 3; type++ ){

		c2[type] = makeCanvas(Form("c2_%d", type), txt[type].c_str());
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetTicks();

		corrHistPlusPlus[type]->Draw("");
		corrHistMinusMinus[type]->Draw("same");
		corrHistPlusMinus[type]->Draw("same");

		TLegend *w2 = new TLegend(0.5,0.60,0.75,0.8);
		w2->SetLineColor(kWhite);
		w2->SetFillColor(0);
		w2->AddEntry(corrHistPlusPlus[type],"++","P");
		w2->AddEntry(corrHistMinusMinus[type],"--","P");
		w2->AddEntry(corrHistPlusMinus[type],"+-/-+","P");
		w2->Draw("same");

	}

	TCanvas* c3 = makeCanvas("c3", "");
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetTicks();

	hist->GetXaxis()->SetTitleColor(kBlack);
	hist->GetYaxis()->SetRangeUser(-0.009,0.009);
	hist->GetXaxis()->SetRangeUser(0.1,4.8);
	hist->Draw();

	for(int type = 0; type < 3; type++){

		show1[type]->Add(show2[type], +1);
		show1[type]->Scale(0.5);
	}

	show1[0]->SetMarkerColor(kRed);
	show1[1]->SetMarkerColor(kRed);
	show1[0]->SetMarkerStyle(20);
	show3[0]->SetMarkerStyle(21);
	show1[1]->SetMarkerStyle(24);
	show3[1]->SetMarkerStyle(25);

	show1[0]->Draw("Psame");
	show3[0]->Draw("Psame");

	show1[1]->Draw("Psame");
	show3[1]->Draw("Psame");

	TLegend *w3 = new TLegend(0.5,0.60,0.75,0.8);
	w3->SetLineColor(kWhite);
	w3->SetFillColor(0);
	w3->AddEntry(show1[0],"like sign, Pb-going","P");
	w3->AddEntry(show3[0],"unlike sign, Pb-going","P");
	w3->AddEntry(show1[1],"like sign, p-going","P");
	w3->AddEntry(show3[1],"unlike sign, p-going","P");

	w3->Draw("same");
	
}