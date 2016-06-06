#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double dEtaReBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,4.8};
const int NdEtaReBins = sizeof(dEtaReBins) / sizeof(dEtaReBins[0]) - 1;
double dEtaReBinCenter[] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.6,3.1,3.8,4.5};

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double weightedAverage(double a1, double a2, double a3, double eta1, double eta2, double eta3){

	double temp1 = a1*eta1 + a2*eta2 + a3*eta3;
	double temp2 = (a1+a2+a3);

	return temp1/temp2;
}

double weightedAverageError(double a1, double a2, double a3, double etaError1, double etaError2, double etaError3){

	double temp1 = (a1/(a1+a2+a3))*(a1/(a1+a2+a3));
	double temp2 = etaError1*etaError1;
	double temp3 = (a2/(a1+a2+a3))*(a2/(a1+a2+a3));
	double temp4 = etaError2*etaError2;
	double temp5 = (a3/(a1+a2+a3))*(a3/(a1+a2+a3));
	double temp6 = etaError3*etaError3;

	double total = temp1*temp2 + temp3*temp4 + temp5*temp6;

	return sqrt(total);

}

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

	TGaxis::SetMaxDigits(3);

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_GEN_v12.root");
	TFile* file1 = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_v51.root");
	TFile* file2 = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3and4.root");

	TH1D* delEta3p[3][2];

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){

			delEta3p[sign][HF] = (TH1D*) file2->Get(Form("ana/delEta3p_%d_%d",sign,HF));
		}
	}
	
	TH1D* QvsdEta1[48][3][2];
	TH1D* QvsdEta2[48][3][2];

	TH1D* PvsdEta1[48][3];
	TH1D* PvsdEta2[48][3];

	TH1D* QaQb1 = (TH1D*) file->Get("ana/c2_ab");
	TH1D* QaQc1 = (TH1D*) file->Get("ana/c2_ac");
	TH1D* QcQb1 = (TH1D*) file->Get("ana/c2_cb");
	TH1D* Ntrk1 = (TH1D*) file->Get("ana/Ntrk");

	TH1D* QaQb2 = (TH1D*) file1->Get("ana/c2_ab");
	TH1D* QaQc2 = (TH1D*) file1->Get("ana/c2_ac");
	TH1D* QcQb2 = (TH1D*) file1->Get("ana/c2_cb");
	TH1D* Ntrk2 = (TH1D*) file1->Get("ana/Ntrk");

	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){

		  PvsdEta1[deta][sign] = (TH1D*) file->Get( Form("ana/PvsdEta_%d_%d",deta,sign) );
		  PvsdEta2[deta][sign] = (TH1D*) file1->Get( Form("ana/PvsdEta_%d_%d",deta,sign) );

			for(int HF = 0; HF < 2; HF++){
		  
			  QvsdEta1[deta][sign][HF] = (TH1D*) file->Get( Form("ana/QvsdEta_%d_%d_%d",deta,sign,HF) );
			  QvsdEta2[deta][sign][HF] = (TH1D*) file1->Get( Form("ana/QvsdEta_%d_%d_%d",deta,sign,HF) );
		
			}
		}
	}

	double meanQaQb1 = QaQb1->GetMean();
	double meanQaQc1 = QaQc1->GetMean();
	double meanQcQb1 = QcQb1->GetMean();

	double c2_a1 = meanQaQb1*meanQaQc1/meanQcQb1;
	double c2_b1 = meanQaQb1*meanQcQb1/meanQaQc1;
	double c2_ab1 = meanQaQb1;

	double meanQaQb2 = QaQb2->GetMean();
	double meanQaQc2 = QaQc2->GetMean();
	double meanQcQb2 = QcQb2->GetMean();

	double c2_a2 = meanQaQb2*meanQaQc2/meanQcQb2;
	double c2_b2 = meanQaQb2*meanQcQb2/meanQaQc2;
	double c2_ab2 = meanQaQb2;

//correction factor for the V2 values from HF. 
	double v2_1[3];
	v2_1[0] = sqrt(c2_b1  );
	v2_1[1] = sqrt(c2_a1  );
	v2_1[2] = sqrt(c2_ab1 );
	
	double v2_2[3];
	v2_2[0] = sqrt(c2_b2  );
	v2_2[1] = sqrt(c2_a2  );
	v2_2[2] = sqrt(c2_ab2 );

	cout << "v2_1[0]: " << v2_1[0] << endl;
	cout << "v2_1[1]: " << v2_1[1] << endl;
	cout << "v2_1[2]: " << v2_1[2] << endl;

	cout << "v2_2[0]: " << v2_2[0] << endl;
	cout << "v2_2[1]: " << v2_2[1] << endl;
	cout << "v2_2[2]: " << v2_2[2] << endl;

	TH1D* hist1[3][2];
	TH1D* hist2[3][2];

	TH1D* hist3[3];
	TH1D* hist4[3];
	for(int sign = 0; sign < 3; sign++){

		hist3[sign] = new TH1D(Form("hist3_%d",sign),"test", NdEtaBins, dEtaBins);
		hist4[sign] = new TH1D(Form("hist4_%d",sign),"test", NdEtaBins, dEtaBins);

		for(int HF = 0; HF < 2; HF++){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", NdEtaReBins2, dEtaReBins2);
			hist2[sign][HF] = new TH1D(Form("hist2_%d_%d",sign,HF),"test", NdEtaReBins2, dEtaReBins2);
		
		}
	}

	for(int deta = 0; deta < NdEtaReBins2; deta++){
		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				if(deta < 8){

					double weight1 = delEta3p[sign][HF]->GetBinContent( 3*deta+1 );
					double weight2 = delEta3p[sign][HF]->GetBinContent( 3*deta+2 );
					double weight3 = delEta3p[sign][HF]->GetBinContent( 3*deta+3 );
					
					double Q_total_real_dEta1 = QvsdEta1[3*deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error1 = QvsdEta1[3*deta][sign][HF]->GetMeanError();

					double Q_total_real_dEta2 = QvsdEta1[3*deta+1][sign][HF]->GetMean();
					double Q_total_real_dEta_error2 = QvsdEta1[3*deta+1][sign][HF]->GetMeanError();

					double Q_total_real_dEta3 = QvsdEta1[3*deta+2][sign][HF]->GetMean();
					double Q_total_real_dEta_error3 = QvsdEta1[3*deta+2][sign][HF]->GetMeanError();	
					
					double value = weightedAverage(weight1, weight2, weight1, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3 );
					double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );
					
					hist1[sign][HF]->SetBinContent(deta+1, value );
					hist1[sign][HF]->SetBinError(deta+1, error );

					double Q_total_real_dEta1 = QvsdEta2[3*deta][sign][HF]->GetMean();
					double Q_total_real_dEta_error1 = QvsdEta2[3*deta][sign][HF]->GetMeanError();

					double Q_total_real_dEta2 = QvsdEta2[3*deta+1][sign][HF]->GetMean();
					double Q_total_real_dEta_error2 = QvsdEta2[3*deta+1][sign][HF]->GetMeanError();

					double Q_total_real_dEta3 = QvsdEta2[3*deta+2][sign][HF]->GetMean();
					double Q_total_real_dEta_error3 = QvsdEta2[3*deta+2][sign][HF]->GetMeanError();	
					
					double value = weightedAverage(weight1, weight2, weight1, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3 );
					double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );
					
					hist2[sign][HF]->SetBinContent(deta+1, value );
					hist2[sign][HF]->SetBinError(deta+1, error );

				}
				else{

					double Q_total_real_dEta1 = QvsdEta1[27][sign][HF]->GetMean();
					double Q_total_real_dEta_error1 = QvsdEta1[27][sign][HF]->GetMeanError();
					
					double Q_total_real_dEta2 = QvsdEta1[28][sign][HF]->GetMean();
					double Q_total_real_dEta_error2 = QvsdEta1[28][sign][HF]->GetMeanError();

					double weight1 = delEta3p[sign][HF]->GetBinContent( 28 );
					double weight2 = delEta3p[sign][HF]->GetBinContent( 29 );

					double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
					double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
					
					hist1[sign][HF]->SetBinContent(deta+1, value );
					hist1[sign][HF]->SetBinError(deta+1,  error);

					double Q_total_real_dEta1 = QvsdEta2[27][sign][HF]->GetMean();
					double Q_total_real_dEta_error1 = QvsdEta2[27][sign][HF]->GetMeanError();
					
					double Q_total_real_dEta2 = QvsdEta2[28][sign][HF]->GetMean();
					double Q_total_real_dEta_error2 = QvsdEta2[28][sign][HF]->GetMeanError();

					double weight1 = delEta3p[sign][HF]->GetBinContent( 28 );
					double weight2 = delEta3p[sign][HF]->GetBinContent( 29 );

					double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
					double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
					
					hist2[sign][HF]->SetBinContent(deta+1, value );
					hist2[sign][HF]->SetBinError(deta+1,  error);
				}
			}
		}
	}

	for(int deta = 0; deta < NdEtaBins; deta++){

		for(int sign = 0; sign < 3; sign++){
				
			double P_total_real_dEta = PvsdEta1[deta][sign]->GetMean();
			double P_total_real_dEta_error = PvsdEta1[deta][sign]->GetMeanError();
			
			hist3[sign]->SetBinContent(deta+1, P_total_real_dEta);
			hist3[sign]->SetBinError(deta+1, P_total_real_dEta_error);

			double P_total_real_dEta = PvsdEta2[deta][sign]->GetMean();
			double P_total_real_dEta_error = PvsdEta2[deta][sign]->GetMeanError();
			
			hist4[sign]->SetBinContent(deta+1, P_total_real_dEta);
			hist4[sign]->SetBinError(deta+1, P_total_real_dEta_error);

			for(int HF = 0; HF < 2; HF++){

				// double Q_total_real_dEta = QvsdEta1[deta][sign][HF]->GetMean();
				// double Q_total_real_dEta_error = QvsdEta1[deta][sign][HF]->GetMeanError();

				// hist1[sign][HF]->SetBinContent(deta+1, Q_total_real_dEta);
				// hist1[sign][HF]->SetBinError(deta+1,  Q_total_real_dEta_error);
				
				// double Q_total_real_dEta = QvsdEta2[deta][sign][HF]->GetMean();
				// double Q_total_real_dEta_error = QvsdEta2[deta][sign][HF]->GetMeanError();
	
				// hist2[sign][HF]->SetBinContent(deta+1, Q_total_real_dEta);
				// hist2[sign][HF]->SetBinError(deta+1,  Q_total_real_dEta_error);
			}
		}
	}
	
	TH1D* base1 = makeHist("base1","","#Delta#eta", "#LTcos(#phi_{1}-#phi_{2})#GT", 48,0,4.8);
    base1->GetXaxis()->SetTitleColor(kBlack);
    base1->GetYaxis()->SetRangeUser(-0.01,0.02);
    base1->GetYaxis()->SetTitleOffset(1.9);

    TH1D* base2 = (TH1D*) base1->Clone("base4");

	TH1D* base3 = makeHist("base3","","#Delta#eta", "#LTcos(#phi_{1}+#phi_{2}-2#Psi_{RP})#GT", 48,0,4.8);
    base3->GetXaxis()->SetTitleColor(kBlack);
    base3->GetYaxis()->SetRangeUser(-0.0025,0.0045);
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

		TH1D* temp3 = (TH1D*)hist2[0][HF]->Clone("temp3");
		temp3->Add(hist2[1][HF], +1);
		temp3->Scale(0.5);
		temp3->Scale( 1.0/v2_2[HF] );
		temp3->SetMarkerColor(kRed);
		temp3->SetLineColor(kRed);
		temp3->SetMarkerStyle(24);
		
		TH1D* temp = (TH1D*)hist1[0][HF]->Clone("temp");
		temp->Add(hist1[1][HF], +1);
		temp->Scale(0.5);
		temp->Scale( 1.0/v2_1[HF] );
		temp->SetMarkerColor(kRed);
		temp->SetLineColor(kRed);
		temp->SetMarkerStyle(20);
		temp->Draw("Psame");
		temp3->Draw("Psame");

		TH1D* temp4 = (TH1D*)hist2[2][HF]->Clone("temp4");
		temp4->Scale( 1.0/v2_2[HF] );
		temp4->SetMarkerColor(kBlue);
		temp4->SetLineColor(kBlue);
		temp4->SetMarkerStyle(25);		

		TH1D* temp2 = (TH1D*)hist1[2][HF]->Clone("temp2");
		temp2->Scale( 1.0/v2_1[HF] );
		temp2->SetMarkerColor(kBlue);
		temp2->SetLineColor(kBlue);
		temp2->SetMarkerStyle(21);
		
		temp2->Draw("Psame");
		temp4->Draw("Psame");

	}

	TLegend *w2 = new TLegend(0.30,0.55,0.7,0.80);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(43);
    w2->AddEntry(temp, "same GEN");
    w2->AddEntry(temp2, "opposite GEN");
    w2->AddEntry(temp3, "same RECO");
    w2->AddEntry(temp4, "opposite RECO");
    w2->Draw("same");

    TLine* l1 = new TLine(0,0.00007,4.8,0.00007);
    l1->SetLineColor(kBlack);
    l1->SetLineStyle(2);

    TCanvas* c3 = makeMultiCanvas("c3","c3",2,1);
    for(int HF = 0; HF < 2; HF++){
		c3->cd(HF+1);
		gPad->SetTicks();
		gPad->SetLeftMargin(0.20);
		gPad->SetBottomMargin(0.16);
		
		if(HF == 0) {TH1D* base5 = (TH1D*) base3->Clone();base5->GetYaxis()->SetRangeUser(-0.002, 0.002);base5->SetTitle("Pb-going");base5->Draw();base5->GetYaxis()->SetTitle("RECO - GEN");}
		if(HF == 1) {TH1D* base6 = (TH1D*) base4->Clone();base6->GetYaxis()->SetRangeUser(-0.002, 0.002);base6->SetTitle("p-going");base6->Draw();base6->GetYaxis()->SetTitle("RECO - GEN");}

		TH1D* ratio3 = (TH1D*)hist2[0][HF]->Clone("ratio3");
		ratio3->Add(hist2[1][HF], +1);
		ratio3->Scale(0.5);
		ratio3->Scale( 1.0/v2_2[HF] );
		ratio3->SetMarkerColor(kRed);
		ratio3->SetLineColor(kRed);
		ratio3->SetMarkerStyle(24);
		
		TH1D* ratio = (TH1D*)hist1[0][HF]->Clone("ratio");
		ratio->Add(hist1[1][HF], +1);
		ratio->Scale(0.5);
		ratio->Scale( 1.0/v2_1[HF] );
		ratio->SetMarkerColor(kRed);
		ratio->SetLineColor(kRed);
		ratio->SetMarkerStyle(20);
		
		ratio3->Add(ratio, -1);
		ratio3->Draw("Psame");
		//l1->Draw("same");

		TH1D* ratio4 = (TH1D*)hist2[2][HF]->Clone("ratio4");
		ratio4->Scale( 1.0/v2_2[HF] );
		ratio4->SetMarkerColor(kBlue);
		ratio4->SetLineColor(kBlue);
		ratio4->SetMarkerStyle(25);		

		TH1D* ratio2 = (TH1D*)hist1[2][HF]->Clone("ratio2");
		ratio2->Scale( 1.0/v2_1[HF] );
		ratio2->SetMarkerColor(kBlue);
		ratio2->SetLineColor(kBlue);
		ratio2->SetMarkerStyle(21);
		
		ratio4->Add(ratio2, -1);
		ratio4->Draw("Psame");
		//l1->Draw("same");


	}    

	TLegend *w3 = new TLegend(0.50,0.65,0.8,0.80);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(ratio3, "same");
    w3->AddEntry(ratio4, "opposite");

    w3->Draw("same");


    TCanvas* c4 = makeMultiCanvas("c4","c4",2,1);

    for(int HF = 0; HF < 2; HF++){
	
	c4->cd(HF+1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);

	if(HF == 0) {TH1D* base7 = (TH1D*) base3->Clone();base7->GetYaxis()->SetRangeUser(-0.002, 0.002);base7->SetTitle("Pb-going");base7->Draw();base7->GetYaxis()->SetTitle("oppo - same");}
	if(HF == 1) {TH1D* base8 = (TH1D*) base4->Clone();base8->GetYaxis()->SetRangeUser(-0.002, 0.002);base8->SetTitle("p-going");base8->Draw();base8->GetYaxis()->SetTitle("oppo - same");}

	TH1D* ratio1 = (TH1D*)hist2[0][HF]->Clone("ratio1");
	ratio1->Add(hist2[1][HF], +1);
	ratio1->Scale(0.5);
	ratio1->Scale( 1.0/v2_2[HF] );
	ratio1->SetMarkerColor(kRed);
	ratio1->SetLineColor(kRed);
	ratio1->SetMarkerStyle(24);

	TH1D* ratio11 = (TH1D*)hist2[2][HF]->Clone("ratio11");

	ratio11->Add(ratio1, -1);
	ratio11->SetMarkerStyle(20);
	ratio11->SetMarkerColor(kRed);
	ratio11->Draw("Psame");
	
	TH1D* ratio5 = (TH1D*)hist1[0][HF]->Clone("ratio5");
	ratio5->Add(hist1[1][HF], +1);
	ratio5->Scale(0.5);
	ratio5->Scale( 1.0/v2_1[HF] );
	ratio5->SetMarkerColor(kRed);
	ratio5->SetLineColor(kRed);
	ratio5->SetMarkerStyle(20);
	
	TH1D* ratio22 = (TH1D*)hist1[2][HF]->Clone("ratio22");

	ratio22->Add(ratio5, -1);
	ratio22->SetMarkerStyle(21);
	ratio22->SetMarkerColor(kBlue);
	ratio22->Draw("Psame");

	}

	TLegend *w4 = new TLegend(0.50,0.65,0.8,0.80);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(20);
    w4->SetTextFont(43);
    w4->AddEntry(ratio11, "GEN");
    w4->AddEntry(ratio22, "RECO");

    w4->Draw("same");


    TCanvas* c5 = makeMultiCanvas("c5","c5",2,1);

    for(int HF = 0; HF < 2; HF++){
	
	c5->cd(HF+1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);

	if(HF == 0) {TH1D* base7 = (TH1D*) base3->Clone();base7->GetYaxis()->SetRangeUser(-0.002, 0.002);base7->SetTitle("Pb-going");base7->Draw();base7->GetYaxis()->SetTitle("oppo - same");}
	if(HF == 1) {TH1D* base8 = (TH1D*) base4->Clone();base8->GetYaxis()->SetRangeUser(-0.002, 0.002);base8->SetTitle("p-going");base8->Draw();base8->GetYaxis()->SetTitle("oppo - same");}

	TH1D* ratio1 = (TH1D*)hist2[0][HF]->Clone("ratio1");
	ratio1->Add(hist2[1][HF], +1);
	ratio1->Scale(0.5);
	ratio1->Scale( 1.0/v2_2[HF] );
	ratio1->SetMarkerColor(kRed);
	ratio1->SetLineColor(kRed);
	ratio1->SetMarkerStyle(24);

	TH1D* ratio11 = (TH1D*)hist2[2][HF]->Clone("ratio11");

	ratio11->Add(ratio1, -1);
	ratio11->SetMarkerStyle(20);
	ratio11->SetMarkerColor(kRed);
	
	TH1D* ratio5 = (TH1D*)hist1[0][HF]->Clone("ratio5");
	ratio5->Add(hist1[1][HF], +1);
	ratio5->Scale(0.5);
	ratio5->Scale( 1.0/v2_1[HF] );
	ratio5->SetMarkerColor(kRed);
	ratio5->SetLineColor(kRed);
	ratio5->SetMarkerStyle(20);
	
	TH1D* ratio22 = (TH1D*)hist1[2][HF]->Clone("ratio22");

	ratio22->Add(ratio5, -1);
	ratio22->SetMarkerStyle(21);
	ratio22->SetMarkerColor(kBlue);
	


	ratio11->Add(ratio22, -1);
	ratio11->Draw("Psame");

	}

	return;


    TCanvas* c2 = makeCanvas("c2","c2");
	gPad->SetTicks();
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.16);

	base1->Draw();
	
	TH1D* temp = (TH1D*)hist3[0]->Clone("temp");
	temp->Add(hist3[1], +1);
	temp->Scale(0.5);

	temp->SetMarkerColor(kRed);
	temp->SetLineColor(kRed);
	temp->SetMarkerStyle(20);
	temp->Draw("Psame");

	TH1D* temp2 = (TH1D*)hist3[2]->Clone("temp2");
	temp2->SetMarkerColor(kBlue);
	temp2->SetLineColor(kBlue);
	temp2->SetMarkerStyle(21);
	temp2->Draw("Psame");

	TH1D* temp3 = (TH1D*)hist4[0]->Clone("temp3");
	temp3->Add(hist4[1], +1);
	temp3->Scale(0.5);

	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);
	temp3->SetMarkerStyle(24);
	temp3->Draw("Psame");

	TH1D* temp4 = (TH1D*)hist4[2]->Clone("temp4");
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);
	temp4->SetMarkerStyle(25);
	temp4->Draw("Psame");
	
	w2->Draw("same");

    // c1->Print("../systematics/closure_40_1000_cumulant.pdf");
    // c3->Print("../systematics/closure_40_1000_diff_cumulant.pdf");

}