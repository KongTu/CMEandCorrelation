#include "RiceStyle.h"

using namespace std;

int Ntrkbins[] = {10,20,30,40,50,60,70,80,90,100,110,120,135,150,165,185,200,220,240,260,300};
double Ntrkbins_fill[] = {10,20,30,40,50,60,70,80,90,100,110,120,135,150,165,185,200,220,240,260,300};
const int nNtrkBins = sizeof(Ntrkbins)/ sizeof(Ntrkbins[0]) - 1;

double get3Real(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;

}

double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}


double get2Real( double R1, double R2, double I1, double I2){

	double real = R1*R2 - I1*I2;
	return real;

}
double get2Imag( double R1, double R2, double I1, double I2){

	double imag = R1*I2 + R2*I1;
	return imag;
}

void plotCorrectionNtrkHistogram(){

	TFile* file = new TFile("../rootfiles/test_trk.root");

	TH1D* QvsNtrk[48][3][2];
	TH1D* XY_real[48][3][2];TH1D* XY_imag[48][3][2];
	TH1D* XZ_real[48][3][2];TH1D* XZ_imag[48][3][2];
	TH1D* YZ_real[48][3][2];TH1D* YZ_imag[48][3][2];
	TH1D* X_real[48][3][2]; TH1D* X_imag[48][3][2];
	TH1D* Y_real[48][3][2]; TH1D* Y_imag[48][3][2];
	TH1D* Z_real[48][3][2]; TH1D* Z_imag[48][3][2];

	TH1D* QaQb = (TH1D*) file->Get("ana/c2_ab");
	TH1D* QaQc = (TH1D*) file->Get("ana/c2_ac");
	TH1D* QcQb = (TH1D*) file->Get("ana/c2_cb");
	TH1D* Ntrk = (TH1D*) file->Get("ana/Ntrk");

	TH1D* aveQ3[2][2];

	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){

			aveQ3[i][j] = (TH1D*)file->Get(Form("ana/aveQ3_%d_%d",i,j) );
		}
	}

	for(int trk = 0; trk < nNtrkBins; trk++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
		  
			  QvsNtrk[trk][sign][HF] = (TH1D*) file->Get( Form("ana/QvsNtrk_%d_%d_%d",trk,sign,HF) );
			  
			  XY_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/XY_real_%d_%d_%d",trk,sign,HF) );
			  XZ_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/XZ_real_%d_%d_%d",trk,sign,HF) );
			  YZ_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/YZ_real_%d_%d_%d",trk,sign,HF) );
			  
			  XY_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/XY_imag_%d_%d_%d",trk,sign,HF) );
			  XZ_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/XZ_imag_%d_%d_%d",trk,sign,HF) );
			  YZ_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/YZ_imag_%d_%d_%d",trk,sign,HF) );
			  
			  X_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/X_real_%d_%d_%d",trk,sign,HF) );
			  Y_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/Y_real_%d_%d_%d",trk,sign,HF) );
			  Z_real[trk][sign][HF] = (TH1D*) file->Get( Form("ana/Z_real_%d_%d_%d",trk,sign,HF) );
			  
			  X_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/X_imag_%d_%d_%d",trk,sign,HF) );
			  Y_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/Y_imag_%d_%d_%d",trk,sign,HF) );
			  Z_imag[trk][sign][HF] = (TH1D*) file->Get( Form("ana/Z_imag_%d_%d_%d",trk,sign,HF) );
			}
		}
	}

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;

	double bCorr = (aveQ3[0][0]->GetMean() * aveQ3[0][0]->GetMean()) +  ( aveQ3[0][1]->GetMean() * aveQ3[0][1]->GetMean() );
	double aCorr = (aveQ3[1][0]->GetMean() * aveQ3[1][0]->GetMean()) +  ( aveQ3[1][1]->GetMean() * aveQ3[1][1]->GetMean() );
	
	double m1 = (aveQ3[0][0]->GetMean() + aveQ3[1][0]->GetMean())/2.0;
	double m2 = (aveQ3[0][1]->GetMean() + aveQ3[1][1]->GetMean())/2.0;
	double abCorr = m1*m1 + m2*m2;

//correction factor for the V2 values from HF. 
	double v2[3];
	v2[0] = sqrt(c2_b - bCorr);
	v2[1] = sqrt(c2_a - aCorr );
	v2[2] = sqrt(c2_ab - abCorr );

	// v2[0] = sqrt(c2_b );
	// v2[1] = sqrt(c2_a  );
	// v2[2] = sqrt(c2_ab );
	cout << "v2[0]: " << v2[0] << endl;
	cout << "v2[1]: " << v2[1] << endl;
	cout << "v2[2]: " << v2[2] << endl;

	TH1D* hist1[3][2];
	TH1D* hist2[3][2];
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", nNtrkBins, Ntrkbins_fill);
			hist2[sign][HF] = new TH1D(Form("hist2_%d_%d",sign,HF),"test", nNtrkBins, Ntrkbins_fill);
		}
	}

	for(int trk = 0; trk < nNtrkBins; trk++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				double Q_total_real_Ntrk = QvsNtrk[trk][sign][HF]->GetMean();
				double Q_total_real_Ntrk_error = QvsNtrk[trk][sign][HF]->GetMeanError();

				cout << "Q_total_real_Ntrk_error: " <<Q_total_real_Ntrk_error<<endl;
				
				double XY_real_temp = 0.; double XZ_real_temp = 0.; double YZ_real_temp = 0.; double X_real_temp = 0.; double Y_real_temp = 0.; double Z_real_temp = 0.;
				double XY_imag_temp = 0.; double XZ_imag_temp = 0.; double YZ_imag_temp = 0.; double X_imag_temp = 0.; double Y_imag_temp = 0.; double Z_imag_temp = 0.;

				XY_real_temp = XY_real[trk][sign][HF]->GetMean();
				XZ_real_temp = XZ_real[trk][sign][HF]->GetMean();
				YZ_real_temp = YZ_real[trk][sign][HF]->GetMean();

				XY_imag_temp = XY_imag[trk][sign][HF]->GetMean();
				XZ_imag_temp = XZ_imag[trk][sign][HF]->GetMean();
				YZ_imag_temp = YZ_imag[trk][sign][HF]->GetMean();

				X_real_temp = X_real[trk][sign][HF]->GetMean();
				Y_real_temp = Y_real[trk][sign][HF]->GetMean();
				Z_real_temp = Z_real[trk][sign][HF]->GetMean();

				X_imag_temp = X_imag[trk][sign][HF]->GetMean();
				Y_imag_temp = Y_imag[trk][sign][HF]->GetMean();
				Z_imag_temp = Z_imag[trk][sign][HF]->GetMean();

				double term1_real = get2Real(XY_real_temp, Z_real_temp, XY_imag_temp, Z_imag_temp);
				double term2_real = get2Real(XZ_real_temp, Y_real_temp, XZ_imag_temp, Y_imag_temp);
				double term3_real = get2Real(YZ_real_temp, X_real_temp, YZ_imag_temp, X_imag_temp);
				double term4_real = get3Real(X_real_temp, Y_real_temp, Z_real_temp, X_imag_temp, Y_imag_temp, Z_imag_temp);
				double all_real = term1_real+term2_real+term3_real-2*term4_real;

				if( Q_total_real_Ntrk == 0.00000 ) continue;

				cout << "before correction is " << Q_total_real_Ntrk << endl;
				cout << "correction factor is: " << all_real << endl;
				cout << "after correction is " << Q_total_real_Ntrk - all_real << endl;
				cout << "-------- " << endl;
				hist1[sign][HF]->SetBinContent(trk+1, Q_total_real_Ntrk );
				hist1[sign][HF]->SetBinError(trk+1,  Q_total_real_Ntrk_error);
				hist2[sign][HF]->SetBinContent(trk+1, Q_total_real_Ntrk - all_real);
				hist2[sign][HF]->SetBinError(trk+1, Q_total_real_Ntrk_error);
			}
		}
	}

	TH1D* base3 = makeHist("base3","","N^{offline}_{trk}", "cos(#phi_{1}+#phi_{2}-2#phi_{3})/v2_{3}", 300,0,300);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.01,0.01);
    base3->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base4 = (TH1D*) base3->Clone("base4");

    TCanvas* c4 = makeMultiCanvas("c4","c4",2,1);
    for(int HF = 0; HF < 2; HF++){
		c4->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		
		if(HF == 0) {base3->SetTitle("Pb-going");base3->Draw();}
		if(HF == 1) {base4->SetTitle("p-going");base4->Draw();}

		TH1D* temp = (TH1D*)hist2[0][HF]->Clone("temp");
		temp->Add(hist2[1][HF], +1);
		temp->Scale(0.5);

		temp->Scale( 1.0/v2[HF] );
		temp->SetMarkerColor(kRed);
		temp->SetLineColor(kRed);
		temp->SetMarkerStyle(20);
		temp->Draw("Psame");

		TH1D* temp2 = (TH1D*)hist2[2][HF]->Clone("temp2");
		temp2->Scale( 1.0/v2[HF] );
		temp2->SetMarkerColor(kBlue);
		temp2->SetLineColor(kBlue);
		temp2->SetMarkerStyle(24);
		temp2->Draw("Psame");

	}

	TLegend *w2 = new TLegend(0.50,0.65,0.8,0.80);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp, "like sign");
    w2->AddEntry(temp2, "unlike sign");
    w2->Draw("same");


	TLegend *w1 = new TLegend(0.50,0.65,0.8,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(hist2[0][0],"corrected","P");
    w1->AddEntry(hist1[0][0],"uncorrected","P");

    TH1D* base = makeHist("base","", "N^{offline}_{trk}", "cos(#phi_{1}+#phi_{2}-2#phi_{3})", 300,0,300);
    base->GetXaxis()->SetTitleColor(kBlack);
    base->GetYaxis()->SetRangeUser(-0.05,0.05);
    base->GetYaxis()->SetTitleOffset(1.9);

	TCanvas* c1 = makeMultiCanvas("c1","c1",2,1);
	for(int HF = 0; HF < 2; HF++){
		c1->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		base->SetTitle("like sign(++)");
		base->Draw();

		hist1[0][HF]->SetMarkerColor(kBlack);
		hist1[0][HF]->SetLineColor(kBlack);
		hist1[0][HF]->SetMarkerStyle(20);
		hist1[0][HF]->Draw("Psame");

		hist2[0][HF]->SetMarkerColor(kRed);
		hist2[0][HF]->SetLineColor(kRed);
		hist2[0][HF]->SetMarkerStyle(20);
		hist2[0][HF]->Draw("Psame");

	}

	w1->Draw("same");


	TCanvas* c2 = makeMultiCanvas("c2","c2",2,1);
	TH1D* base1 = (TH1D*) base->Clone("base1");
	for(int HF = 0; HF < 2; HF++){
		c2->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		base1->SetTitle("like sign(--)");
		base1->Draw();

		hist1[1][HF]->SetMarkerColor(kBlack);
		hist1[1][HF]->SetLineColor(kBlack);
		hist1[1][HF]->SetMarkerStyle(20);
		hist1[1][HF]->Draw("Psame");
		
		hist2[1][HF]->SetMarkerColor(kRed);
		hist2[1][HF]->SetLineColor(kRed);
		hist2[1][HF]->SetMarkerStyle(20);
		hist2[1][HF]->Draw("Psame");
	
	}
	w1->Draw("same");

	TCanvas* c3 = makeMultiCanvas("c3","c3",2,1);
	TH1D* base2 = (TH1D*) base->Clone("base2");
	for(int HF = 0; HF < 2; HF++){
		c3->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		base2->SetTitle("unlike sign(+-/-+)");
		base2->Draw();

		hist1[2][HF]->SetMarkerColor(kBlack);
		hist1[2][HF]->SetLineColor(kBlack);
		hist1[2][HF]->SetMarkerStyle(20);
		hist1[2][HF]->Draw("Psame");
		
		hist2[2][HF]->SetMarkerColor(kRed);
		hist2[2][HF]->SetLineColor(kRed);
		hist2[2][HF]->SetMarkerStyle(20);
		hist2[2][HF]->Draw("Psame");
	}	
	w1->Draw("same");

}