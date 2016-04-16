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

void plotCorrectionGenNested(){

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_GEN_NestedLoop_v3.root");

	TH1D* QvsdEta[48];

	for(int deta = 0; deta < NdEtaBins; deta++){
		  
	  QvsdEta[deta] = (TH1D*) file->Get( Form("ana/QvsdEta_%d",deta) );
			
	}


	TH1D* hist1[3][2];
	TH1D* hist2[3][2];
	for(int sign = 0; sign < 1; sign++){
		for(int HF = 0; HF < 1; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", NdEtaBins, dEtaBins);
			hist2[sign][HF] = new TH1D(Form("hist2_%d_%d",sign,HF),"test", NdEtaBins, dEtaBins);
		}
	}

	for(int deta = 0; deta < NdEtaBins; deta++){

		double Q_total_real_dEta = QvsdEta[deta]->GetMean();
		double Q_total_real_dEta_error = QvsdEta[deta]->GetMeanError();

		hist1[0][0]->SetBinContent(deta+1, Q_total_real_dEta );
		hist1[0][0]->SetBinError(deta+1,  Q_total_real_dEta_error);
	
	}

	TH1D* base3 = makeHist("base3","like sign(++)","#Delta#eta", "cos(#phi_{1}+#phi_{2}-2#phi_{3})/v2_{3}", 48,0,4.8);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.0015,0.001);
    base3->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base4 = (TH1D*) base3->Clone("base4");

    TCanvas* c4 = makeCanvas("c4","c4");
	gPad->SetTicks();
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);
	base3->Draw();
	TH1D* temp = (TH1D*)hist1[0][0]->Clone("temp");
	temp->Scale(1.0/0.20);
	temp->SetMarkerColor(kRed);
	temp->SetLineColor(kRed);
	temp->SetMarkerStyle(20);
	temp->Draw("Psame");


}