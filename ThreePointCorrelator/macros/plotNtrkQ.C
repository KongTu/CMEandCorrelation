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


void plotNtrkQ(){

	TFile* file = new TFile("../rootfiles/CME_Qvsntrk_pPb_HM_185_220_v3.root");
	
	vector<TH2D*> QPlusPlus = loading2DHistogram(file, "ana/QvsNtrkPlusPlus_", 3);
	vector<TH2D*> QMinusMinus = loading2DHistogram(file, "ana/QvsNtrkMinusMinus_", 3);
	vector<TH2D*> QPlusMinus = loading2DHistogram(file, "ana/QvsNtrkPlusMinus_", 3);
	
	TH1D* QaQb = file->Get("ana/c2_ab");
	TH1D* QaQc = file->Get("ana/c2_ac");
	TH1D* QcQb = file->Get("ana/c2_cb");
	TH1D* Ntrk = file->Get("ana/Ntrk");

	ifstream ifile("ntrk_185_220.txt");
	if(ifile.is_open() ){

		vector<double> t;
		double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14;
		ifile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9 >> v10 >> v11 >> v12 >> v13 >> v14;
		t.push_back(v1); t.push_back(v2);t.push_back(v3);t.push_back(v4);
		t.push_back(v5); t.push_back(v6);t.push_back(v7);t.push_back(v8);
		t.push_back(v9); t.push_back(v10);t.push_back(v11);t.push_back(v12);
		t.push_back(v13);t.push_back(v14);
	}



	for( int type = 0; type < 3; type ++ ){

		ntrk1[type] = makeHistDifferentBins(Form("ntrk1_%d",type),"","", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", Nntrkbins, ntrkbins, kBlack);
		ntrk2[type] = makeHistDifferentBins(Form("ntrk2_%d",type),"","", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", Nntrkbins, ntrkbins, kRed);
		ntrk3[type] = makeHistDifferentBins(Form("ntrk3_%d",type),"","", "<cos(#phi_{1}+#phi_{2}-2#phi_{3})>/v_{2,3}", Nntrkbins, ntrkbins, kBlue);
		
		for(int trk = 0; trk < Nntrkbins; trk++ ){

			Qplusplus1Dntrk[type][trk] = (TH1D*) QPlusPlus[type]->ProjectionY(Form("Qplusplus1Dtrk_%d_%d", type, trk), trk+1, trk+2);
			Qminusminus1Dntrk[type][trk] = (TH1D*) QMinusMinus[type]->ProjectionY(Form("Qminusminus1Dtrk_%d_%d", type, trk), trk+1, trk+2);
			Qplusminus1Dntrk[type][trk] = (TH1D*) QPlusMinus[type]->ProjectionY(Form("Qplusminus1Dtrk_%d_%d", type, trk), trk+1, trk+2);

			double Qplusplus1DMean = Qplusplus1Dntrk[type][trk]->GetMean();
			double Qminusminus1DMean = Qminusminus1Dntrk[type][trk]->GetMean();
			double Qplusminus1DMean = Qplusminus1Dntrk[type][trk]->GetMean();

			double Qplusplus1DMeanErr = Qplusplus1Dntrk[type][trk]->GetMeanError();
			double Qminusminus1DMeanErr = Qminusminus1Dntrk[type][trk]->GetMeanError();
			double Qplusminus1DMeanErr = Qplusminus1Dntrk[type][trk]->GetMeanError();

			ntrk1[type]->SetBinContent(trk+1, Qplusplus1DMean );
			ntrk1[type]->SetBinError(trk+1, Qplusplus1DMeanErr );
			
			ntrk2[type]->SetBinContent(trk+1, Qminusminus1DMean );
			ntrk2[type]->SetBinError(trk+1, Qminusminus1DMeanErr );

			ntrk3[type]->SetBinContent(trk+1, Qplusminus1DMean );
			ntrk3[type]->SetBinError(trk+1, Qplusminus1DMeanErr );
		}

	}

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;

	double v2[3];
	v2[0] = sqrt( c2_b );
	v2[1] = sqrt( c2_a );
	v2[2] = sqrt( c2_ab );
	
	string txt[3] = {"HF+", "HF-", "HF+/-"};


	TCanvas* c1[3];
	for(int type = 0; type < 3; type++ ){

		c1[type] = makeCanvas(Form("c1_%d", type), txt[type].c_str());
		gPad->SetLeftMargin(0.15);
		gPad->SetBottomMargin(0.15);
		gPad->SetTicks();

		hist2->GetXaxis()->SetTitleColor(kBlack);
		hist2->GetYaxis()->SetRangeUser(-0.0009,0.0009);
		hist2->Draw();

		ntrk1[type]->Scale( 1.0/v2[type] );
		ntrk2[type]->Scale( 1.0/v2[type] );
		ntrk3[type]->Scale( 1.0/v2[type] );

		ntrk1[type]->Draw("Psame");
		ntrk2[type]->Draw("Psame");
		ntrk3[type]->Draw("Psame");

		TLegend *w1 = new TLegend(0.5,0.20,0.75,0.43);
		w1->SetLineColor(kWhite);
		w1->SetFillColor(0);
		w1->AddEntry(ntrk1[type],"++","P");
		w1->AddEntry(ntrk2[type],"--","P");
		w1->AddEntry(ntrk3[type],"+-/-+","P");
		w1->Draw("same");
	}















}