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

void plotGenCorrelation(){

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_GEN_v4.root");

	TH1D* QvsdEta[48][3][2];
	TH1D* PvsdEta[48][3];

	TH1D* QaQb = (TH1D*) file->Get("ana/c2_ab");
	TH1D* QaQc = (TH1D*) file->Get("ana/c2_ac");
	TH1D* QcQb = (TH1D*) file->Get("ana/c2_cb");
	TH1D* Ntrk = (TH1D*) file->Get("ana/Ntrk");

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){

		  PvsdEta[deta][sign] = (TH1D*) file->Get( Form("ana/PvsdEta_%d_%d",deta,sign) );

			for(int HF = 0; HF < 2; HF++){
		  
			  QvsdEta[deta][sign][HF] = (TH1D*) file->Get( Form("ana/QvsdEta_%d_%d_%d",deta,sign,HF) );
		
			}
		}
	}

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;


//correction factor for the V2 values from HF. 
	double v2[3];
	v2[0] = sqrt(c2_b  );
	v2[1] = sqrt(c2_a  );
	v2[2] = sqrt(c2_ab );

	cout << "v2[0]: " << v2[0] << endl;
	cout << "v2[1]: " << v2[1] << endl;
	cout << "v2[2]: " << v2[2] << endl;

	TH1D* hist1[3][2];
	TH1D* hist2[3];
	for(int sign = 0; sign < 3; sign++){

		hist2[sign] = new TH1D(Form("hist2_%d",sign),"test", NdEtaBins, dEtaBins);

		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", NdEtaBins, dEtaBins);
		}
	}

	for(int deta = 0; deta < NdEtaBins; deta++){

		for(int sign = 0; sign < 3; sign++){
				
			double P_total_real_dEta = PvsdEta[deta][sign]->GetMean();
			double P_total_real_dEta_error = PvsdEta[deta][sign]->GetMeanError();
			
			hist2[sign]->SetBinContent(deta+1, P_total_real_dEta);
			hist2[sign]->SetBinError(deta+1, P_total_real_dEta_error);

			for(int HF = 0; HF < 2; HF++){

				double Q_total_real_dEta = QvsdEta[deta][sign][HF]->GetMean();
				double Q_total_real_dEta_error = QvsdEta[deta][sign][HF]->GetMeanError();

				hist1[sign][HF]->SetBinContent(deta+1, Q_total_real_dEta);
				hist1[sign][HF]->SetBinError(deta+1,  Q_total_real_dEta_error);

			}
		}
	}
	
	TH1D* base1 = makeHist("base1","","#Delta#eta", "#LTcos(#phi_{1}-#phi_{2})#GT", 48,0,4.8);
    base1->GetXaxis()->SetTitleColor(kBlack);
    base1->GetYaxis()->SetRangeUser(-0.01,0.02);
    base1->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base2 = (TH1D*) base1->Clone("base4");

	TH1D* base3 = makeHist("base3","","#Delta#eta", "#LTcos(#phi_{1}+#phi_{2}-2#phi_{3})/v_{2,3}#GT", 48,0,4.8);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.002,0.004);
    base3->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base4 = (TH1D*) base3->Clone("base4");

    TCanvas* c1 = makeMultiCanvas("c1","c1",2,1);
    for(int HF = 0; HF < 2; HF++){
		c1->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		
		if(HF == 0) {base3->SetTitle("Pb-going");base3->Draw();}
		if(HF == 1) {base4->SetTitle("p-going");base4->Draw();}

		TH1D* temp = (TH1D*)hist1[0][HF]->Clone("temp");
		temp->Add(hist1[1][HF], +1);
		temp->Scale(0.5);

		temp->Scale( 1.0/v2[HF] );
		temp->SetMarkerColor(kRed);
		temp->SetLineColor(kRed);
		temp->SetMarkerStyle(20);
		temp->Draw("Psame");

		TH1D* temp2 = (TH1D*)hist1[2][HF]->Clone("temp2");
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

    TCanvas* c2 = makeCanvas("c2","c2");
	gPad->SetTicks();
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);

	base1->Draw();
	
	TH1D* temp = (TH1D*)hist2[0]->Clone("temp");
	temp->Add(hist2[1], +1);
	temp->Scale(0.5);

	temp->SetMarkerColor(kRed);
	temp->SetLineColor(kRed);
	temp->SetMarkerStyle(20);
	temp->Draw("Psame");

	TH1D* temp2 = (TH1D*)hist2[2]->Clone("temp2");
	temp2->SetMarkerColor(kBlue);
	temp2->SetLineColor(kBlue);
	temp2->SetMarkerStyle(24);
	temp2->Draw("Psame");

	
	
	w2->Draw("same");

}