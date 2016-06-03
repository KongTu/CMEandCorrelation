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

	gStyle->SetErrorX(0);
	TGaxis::SetMaxDigits(3);

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}

	double dEtaBinsCenter[50];
	for(int deta = 0; deta < NdEtaBins; deta++){

		double bincenter = (dEtaBins[deta+1] - dEtaBins[deta])/2.0 + dEtaBins[deta];
		dEtaBinsCenter[deta] = bincenter;

	}


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_PbPb_30_40_v38.root");

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

	TH1D* base3 = makeHist("base3","","#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.0007,0.0005);
    base3->GetXaxis()->SetRangeUser(0,1.6);
    base3->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base4 = (TH1D*) base3->Clone("base4");

    TCanvas* c4 = makeCanvas("c4","c4");

    //TCanvas* c4 = makeMultiCanvas("c4","c4",2,1);
    for(int HF = 0; HF < 1; HF++){
		c4->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		
		if(HF == 0) {base3->SetTitle("");base3->Draw();}
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
		temp1->SetBinContent(1,100);
		temp1->SetMarkerColor(kRed);
		temp1->SetLineColor(kRed);
		temp1->SetMarkerStyle(20);
		temp1->Draw("Psame");

		TH1D* temp5 = (TH1D*)hist1[2][0]->Clone("temp5");
		TH1D* temp6 = (TH1D*)hist1[2][1]->Clone("temp6");
		temp5->Add(temp6, +1);
		temp5->Scale(0.5);
		temp5->Scale( 1.0/v2[2] );
		temp5->SetBinContent(1,100);
		temp5->SetMarkerColor(kBlue);
		temp5->SetLineColor(kBlue);
		temp5->SetMarkerStyle(24);
		temp5->Draw("Psame");

	}

	TLegend *w1 = new TLegend(0.50,0.65,0.8,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(hist2[0][0],"corrected","P");
    w1->AddEntry(hist1[0][0],"uncorrected","P");

    double value1[50];
    double value1_error[50];
    double value2[50];
    double value2_error[50];

    for(int deta = 0; deta < NdEtaBins; deta++){


    	value1[deta] = temp1->GetBinContent(deta+1);
    	value1_error[deta] = temp1->GetBinContent(deta+1);

    	value2[deta] = temp5->GetBinContent(deta+1);
    	value2_error[deta] = temp5->GetBinContent(deta+1);
    }

    TBox *box1[50];
    TBox *box2[50];

    for(int deta = 1; deta < 16; deta++){

    	double xe = 0.005;
    	double ye = 0.00006;

    	box1[deta] = new TBox(dEtaBinsCenter[deta]-xe,value1[deta]-ye,dEtaBinsCenter[deta]+xe,value1[deta]+ye);
		box1[deta]->SetFillColor(kRed);
        box1[deta]->SetFillStyle(0);
    	box1[deta]->SetLineWidth(1);
    	box1[deta]->SetLineColor(kRed);
        box1[deta]->Draw("SAME");

		box2[deta] = new TBox(dEtaBinsCenter[deta]-xe,value2[deta]-ye,dEtaBinsCenter[deta]+xe,value2[deta]+ye);
		box2[deta]->SetFillColor(kBlue);
        box2[deta]->SetFillStyle(0);
    	box2[deta]->SetLineWidth(1);
    	box2[deta]->SetLineColor(kBlue);
        box2[deta]->Draw("SAME");
    }

	double alicebins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6};
	double aliceBincenter[] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5};
	double y1[] = {0.000013,-0.000044,-0.000037,-0.000033, -0.000025,-0.000063,-0.000011,-0.000208};
	double y1_error[] = {0.000023, 0.000024, 0.000026, 0.000029,0.000033,0.000039,0.000051,0.000101};
	double y2[] = {-0.000270, -0.000245, -0.000198, -0.000138, -0.000070, -0.000047,-0.000064, 0.000018 };
	double y2_error[] = {0.000016, 0.000017, 0.000019, 0.000021, 0.000023, 0.000028, 0.000037, 0.000073 };


	TH1D* alice1 = new TH1D("alice1", "alice1", 8, alicebins );
	TH1D* alice2 = new TH1D("alice2", "alice2", 8, alicebins );

	for(int deta = 0; deta < 8; deta++){

		alice1->SetBinContent(deta+1, y1[deta]);
		alice1->SetBinError(deta+1, y1_error[deta]);
		
		alice2->SetBinContent(deta+1, y2[deta]);
		alice2->SetBinError(deta+1, y2_error[deta]);

	}

	alice1->SetMarkerColor(kGreen+3);
	alice1->SetLineColor(kGreen+3);
	alice1->SetMarkerStyle(21);

	alice2->SetMarkerColor(kBlack);
	alice2->SetLineColor(kBlack);
	alice2->SetMarkerStyle(21);

	alice1->Draw("Psame");
	alice2->Draw("Psame");


    TBox *box3[50];
    TBox *box4[50];

    for(int deta = 0; deta < 8; deta++){

    	double xe = 0.005;
    	double ye = 0.0001;

    	box3[deta] = new TBox(aliceBincenter[deta]-xe,y1[deta]-ye,aliceBincenter[deta]+xe,y1[deta]+ye);
		box3[deta]->SetFillColor(kGreen+3);
        box3[deta]->SetFillStyle(0);
    	box3[deta]->SetLineWidth(1);
    	box3[deta]->SetLineColor(kGreen+3);
        box3[deta]->Draw("SAME");

		box4[deta] = new TBox(aliceBincenter[deta]-xe,y2[deta]-ye,aliceBincenter[deta]+xe,y2[deta]+ye);
		box4[deta]->SetFillColor(kBlack);
        box4[deta]->SetFillStyle(0);
    	box4[deta]->SetLineWidth(1);
    	box4[deta]->SetLineColor(kBlack);
        box4[deta]->Draw("SAME");
    }

	TLegend *w2 = new TLegend(0.50,0.2,0.8,0.40);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp1, "CMS same");
    w2->AddEntry(temp5, "CMS opposite");
 	w2->AddEntry(alice2, "ALICE same");
    w2->AddEntry(alice1, "ALICE opposite");   
    w2->Draw("same");

	TLatex* r4 = new TLatex(0.24, 0.82, "PbPb #sqrt{s_{NN}} = 2.76 TeV");
	r4->SetNDC();
	r4->SetTextSize(23);
	r4->SetTextFont(43);
	r4->SetTextColor(kBlack);
	r4->Draw("same");

	TLatex* r5 = new TLatex(0.74, 0.82, "30-40%");
	r5->SetNDC();
	r5->SetTextSize(23);
	r5->SetTextFont(43);
	r5->SetTextColor(kBlack);
	r5->Draw("same");

   	TLatex* r1 = new TLatex(0.60,0.92, "CMS");
    r1->SetNDC();
    r1->SetTextSize(0.05);
    r1->Draw("same");

    TLatex* r2 = new TLatex(0.73,0.92, "Preliminary");
    r2->SetNDC();
    r2->SetTextSize(22);
    r2->SetTextFont(53);
    r2->Draw("same");
}