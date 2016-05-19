#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;

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

void plotCorrectionCrosschecks(){

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_PbPb_30_40_v31.root");

	TH1D* QvsdEta[48][3][2];
	// TH1D* XY_real[48][3][2];TH1D* XY_imag[48][3][2];
	// TH1D* XZ_real[48][3][2];TH1D* XZ_imag[48][3][2];
	// TH1D* YZ_real[48][3][2];TH1D* YZ_imag[48][3][2];
	// TH1D* X_real[48][3][2]; TH1D* X_imag[48][3][2];
	// TH1D* Y_real[48][3][2]; TH1D* Y_imag[48][3][2];
	// TH1D* Z_real[48][3][2]; TH1D* Z_imag[48][3][2];

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

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
		  
			  QvsdEta[deta][sign][HF] = (TH1D*) file->Get( Form("ana/QvsdEta_%d_%d_%d",deta,sign,HF) );
			  
			  // XY_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/XY_real_%d_%d_%d",deta,sign,HF) );
			  // XZ_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/XZ_real_%d_%d_%d",deta,sign,HF) );
			  // YZ_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/YZ_real_%d_%d_%d",deta,sign,HF) );
			  
			  // XY_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/XY_imag_%d_%d_%d",deta,sign,HF) );
			  // XZ_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/XZ_imag_%d_%d_%d",deta,sign,HF) );
			  // YZ_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/YZ_imag_%d_%d_%d",deta,sign,HF) );
			  
			  // X_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/X_real_%d_%d_%d",deta,sign,HF) );
			  // Y_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/Y_real_%d_%d_%d",deta,sign,HF) );
			  // Z_real[deta][sign][HF] = (TH1D*) file->Get( Form("ana/Z_real_%d_%d_%d",deta,sign,HF) );
			  
			  // X_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/X_imag_%d_%d_%d",deta,sign,HF) );
			  // Y_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/Y_imag_%d_%d_%d",deta,sign,HF) );
			  // Z_imag[deta][sign][HF] = (TH1D*) file->Get( Form("ana/Z_imag_%d_%d_%d",deta,sign,HF) );
			
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
	//v2[0] = sqrt(meanQcQb);
	// v2[0] = sqrt(c2_b - bCorr);
	// v2[1] = sqrt(c2_a - aCorr );
	// v2[2] = sqrt(c2_ab - abCorr );

	v2[0] = sqrt(c2_b );
	v2[1] = sqrt(c2_a  );
	v2[2] = sqrt(c2_ab );
	cout << "v2[0]: " << v2[0] << endl;
	cout << "v2[1]: " << v2[1] << endl;
	cout << "v2[2]: " << v2[2] << endl;

	TH1D* hist1[3][1];
	TH1D* hist2[3][1];
	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", NdEtaBins, dEtaBins);
			hist2[sign][HF] = new TH1D(Form("hist2_%d_%d",sign,HF),"test", NdEtaBins, dEtaBins);
		}
	}

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				double Q_total_real_dEta = QvsdEta[deta][sign][HF]->GetMean();
				double Q_total_real_dEta_error = QvsdEta[deta][sign][HF]->GetMeanError();



				// cout << "Q_total_real_dEta_error: " <<Q_total_real_dEta_error<<endl;
				
				// double XY_real_temp = 0.; double XZ_real_temp = 0.; double YZ_real_temp = 0.; double X_real_temp = 0.; double Y_real_temp = 0.; double Z_real_temp = 0.;
				// double XY_imag_temp = 0.; double XZ_imag_temp = 0.; double YZ_imag_temp = 0.; double X_imag_temp = 0.; double Y_imag_temp = 0.; double Z_imag_temp = 0.;

				// XY_real_temp = XY_real[deta][sign][HF]->GetMean();
				// XZ_real_temp = XZ_real[deta][sign][HF]->GetMean();
				// YZ_real_temp = YZ_real[deta][sign][HF]->GetMean();

				// XY_imag_temp = XY_imag[deta][sign][HF]->GetMean();
				// XZ_imag_temp = XZ_imag[deta][sign][HF]->GetMean();
				// YZ_imag_temp = YZ_imag[deta][sign][HF]->GetMean();

				// X_real_temp = X_real[deta][sign][HF]->GetMean();
				// Y_real_temp = Y_real[deta][sign][HF]->GetMean();
				// Z_real_temp = Z_real[deta][sign][HF]->GetMean();

				// X_imag_temp = X_imag[deta][sign][HF]->GetMean();
				// Y_imag_temp = Y_imag[deta][sign][HF]->GetMean();
				// Z_imag_temp = Z_imag[deta][sign][HF]->GetMean();

				// double term1_real = get2Real(XY_real_temp, Z_real_temp, XY_imag_temp, Z_imag_temp);
				// double term2_real = get2Real(XZ_real_temp, Y_real_temp, XZ_imag_temp, Y_imag_temp);
				// double term3_real = get2Real(YZ_real_temp, X_real_temp, YZ_imag_temp, X_imag_temp);
				// double term4_real = get3Real(X_real_temp, Y_real_temp, Z_real_temp, X_imag_temp, Y_imag_temp, Z_imag_temp);
				// double all_real = term1_real+term2_real+term3_real-2*term4_real;

				// cout << "before correction is " << Q_total_real_dEta << endl;
				// cout << "correction factor is: " << all_real << endl;
				// cout << "after correction is " << Q_total_real_dEta - all_real << endl;
				// cout << "-------- " << endl;
				hist1[sign][HF]->SetBinContent(deta+1, Q_total_real_dEta );
				hist1[sign][HF]->SetBinError(deta+1,  Q_total_real_dEta_error);


				// hist2[sign][HF]->SetBinContent(deta+1, Q_total_real_dEta - all_real);
				// hist2[sign][HF]->SetBinError(deta+1, Q_total_real_dEta_error);
			}
		}
	}


	TH1D* base3 = makeHist("base3","","#Delta#eta", "cos(#phi_{1}+#phi_{2}-2#phi_{3})/v2_{3}", 48,0,4.8);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.0006,0.0003);
    base3->GetXaxis()->SetRangeUser(0.1,2.0);
    base3->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base4 = (TH1D*) base3->Clone("base4");

    TCanvas* c4 = makeCanvas("c4","c4");

    //TCanvas* c4 = makeMultiCanvas("c4","c4",2,1);
    for(int HF = 0; HF < 1; HF++){
		c4->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		
		if(HF == 0) {base3->SetTitle("PbPb 30-40%");base3->Draw();}
		if(HF == 1) {base4->SetTitle("p-going");base4->Draw();}

		TH1D* temp1 = (TH1D*)hist1[0][0]->Clone("temp1");
		TH1D* temp2 = (TH1D*)hist1[0][1]->Clone("temp2");
		TH1D* temp3 = (TH1D*)hist1[1][0]->Clone("temp3");
		TH1D* temp4 = (TH1D*)hist1[1][1]->Clone("temp4");

		temp1->Add(hist1[0][1], +1);
		temp1->Add(hist1[1][0], +1);
		temp1->Add(hist1[1][1], +1);

		temp1->Scale(0.25);

		temp1->Scale( 1.0/v2[2] );
		temp1->SetMarkerColor(kRed);
		temp1->SetLineColor(kRed);
		temp1->SetMarkerStyle(20);
		temp1->Draw("Psame");

		TH1D* temp5 = (TH1D*)hist1[2][0]->Clone("temp5");
		TH1D* temp6 = (TH1D*)hist1[2][1]->Clone("temp6");
		temp5->Add(temp6, +1);
		temp5->Scale(0.5);
		temp5->Scale( 1.0/v2[2] );
		temp5->SetMarkerColor(kBlue);
		temp5->SetLineColor(kBlue);
		temp5->SetMarkerStyle(24);
		temp5->Draw("Psame");

	}

	TLegend *w2 = new TLegend(0.50,0.65,0.8,0.80);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp1, "like sign");
    w2->AddEntry(temp5, "unlike sign");
    w2->Draw("same");


	TLegend *w1 = new TLegend(0.50,0.65,0.8,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(hist2[0][0],"corrected","P");
    w1->AddEntry(hist1[0][0],"uncorrected","P");


}