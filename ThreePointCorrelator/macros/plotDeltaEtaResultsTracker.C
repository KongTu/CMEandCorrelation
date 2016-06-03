#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
//rebin option1:
double dEtaReBins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.4,4.2,4.8};
const int NdEtaReBins = sizeof(dEtaReBins) / sizeof(dEtaReBins[0]) - 1;

double dEtaReBinCenter[] = {0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.6,3.1,3.8,4.5};

//rebin option2:
double dEtaReBins2[] = {0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.2,2.8,3.8,4.8};
const int NdEtaReBins2 = sizeof(dEtaReBins2) / sizeof(dEtaReBins2[0]) - 1;

double dEtaReBinCenter2[] = {0.15,0.45,0.75,1.05,1.35,1.65,2.0,2.5,3.3,4.3};

double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
const int Nmults = 3;

double total_systematics_pPb = 0.000077;
double total_systematics_PbPb = 0.000071;

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

void plotDeltaEtaResultsTracker(){

	gStyle->SetErrorX(0);

	TFile *file1 = new TFile("../dataPoints/diff.root");
	TH1D* temp1 = (TH1D*) file1->Get("temp1");
	TH1D* temp2 = (TH1D*) file1->Get("temp2");
	TH1D* temp3 = (TH1D*) file1->Get("temp3");
	TH1D* temp4 = (TH1D*) file1->Get("temp4");
	TH1D* temp5 = (TH1D*) file1->Get("temp5");
	TH1D* temp9 = (TH1D*) file1->Get("temp9");
	TH1D* diff1 = (TH1D*) file1->Get("diff1");
	TH1D* diff2 = (TH1D*) file1->Get("diff2");
	TH1D* diff3 = (TH1D*) file1->Get("diff3");


	TFile* file[2];

	//file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v32_3.root");
	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_Systematics_v9.root");

	TH1D* QvsdEta[48][3];
	TH1D* delEta3p[3];

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){

			delEta3p[sign] = (TH1D*) file[0]->Get(Form("ana/delEta3p_%d",sign));
		}
	}
	
	TH1D* QaQb;
	TH1D* QaQc; 
	TH1D* QcQb;
	TH1D* aveQ3[2][2];

	QaQb = (TH1D*)file[0]->Get("ana/c2_ab");
	QaQc = (TH1D*)file[0]->Get("ana/c2_ac");
	QcQb = (TH1D*)file[0]->Get("ana/c2_cb");

	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){

			aveQ3[i][j] = (TH1D*)file[0]->Get(Form("ana/aveQ3_%d_%d",i,j) );
		}
	}
	
	for(int deta = 0; deta < NdEtaBins; deta++){
		for(int sign = 0; sign < 3; sign++){
		  
		  QvsdEta[deta][sign] = (TH1D*) file[0]->Get( Form("ana/QvsdEta_%d_%d",deta,sign) );
			  
		}
	}
	
	double v2[4];//get corrected v2_3

	double meanQaQb = QaQb->GetMean();
	double meanQaQc = QaQc->GetMean();
	double meanQcQb = QcQb->GetMean();

	double c2_a = meanQaQb*meanQaQc/meanQcQb;
	double c2_b = meanQaQb*meanQcQb/meanQaQc;
	double c2_ab = meanQaQb;
	double c2_c = meanQaQc*meanQcQb/meanQaQb;

	double bCorr = (aveQ3[0][0]->GetMean() * aveQ3[0][0]->GetMean()) +  ( aveQ3[0][1]->GetMean() * aveQ3[0][1]->GetMean() );
	double aCorr = (aveQ3[1][0]->GetMean() * aveQ3[1][0]->GetMean()) +  ( aveQ3[1][1]->GetMean() * aveQ3[1][1]->GetMean() );

	double m1 = (aveQ3[0][0]->GetMean() + aveQ3[1][0]->GetMean())/2.0;
	double m2 = (aveQ3[0][1]->GetMean() + aveQ3[1][1]->GetMean())/2.0;

	double abCorr = m1*m1 + m2*m2;

	v2[0] = sqrt(c2_b - bCorr);
	v2[1] = sqrt(c2_a - aCorr );
	v2[2] = sqrt(c2_ab - abCorr );
	v2[3] = sqrt(c2_c);

	cout << "v2[0]: " << v2[0] << endl;
	cout << "v2[1]: " << v2[1] << endl;
	cout << "v2[2]: " << v2[2] << endl;
	cout << "v2[3]: " << v2[3] << endl;

	TH1D* hist1[3];
	for(int sign = 0; sign < 3; sign++){
		hist1[sign] = new TH1D(Form("hist1_%d",sign),"test", NdEtaReBins2, dEtaReBins2);
	}	
	
	for(int deta = 0; deta < NdEtaReBins2; deta++){
		for(int sign = 0; sign < 3; sign++){

			if(deta < 8){

				double Q_total_real_dEta1 = QvsdEta[3*deta][sign]->GetMean();
				double Q_total_real_dEta_error1 = QvsdEta[3*deta][sign]->GetMeanError();

				double Q_total_real_dEta2 = QvsdEta[3*deta+1][sign]->GetMean();
				double Q_total_real_dEta_error2 = QvsdEta[3*deta+1][sign]->GetMeanError();

				double Q_total_real_dEta3 = QvsdEta[3*deta+2][sign]->GetMean();
				double Q_total_real_dEta_error3 = QvsdEta[3*deta+2][sign]->GetMeanError();

				double weight1 = delEta3p[sign]->GetBinContent( 3*deta+1 );
				double weight2 = delEta3p[sign]->GetBinContent( 3*deta+2 );
				double weight3 = delEta3p[sign]->GetBinContent( 3*deta+3 );
				
				double value = weightedAverage(weight1, weight2, weight3, Q_total_real_dEta1, Q_total_real_dEta2, Q_total_real_dEta3);
				double error = weightedAverageError(weight1, weight2, weight3, Q_total_real_dEta_error1, Q_total_real_dEta_error2, Q_total_real_dEta_error3 );				

				hist1[sign]->SetBinContent(deta+1, value );
				hist1[sign]->SetBinError(deta+1, error );

			}
			if(deta >= 8){

				hist1[sign]->SetBinContent(deta+1, 10000.0 );
				hist1[sign]->SetBinError(deta+1, 0.0 );

			}
			if(deta == 9){

				// double Q_total_real_dEta1 = QvsdEta[27][sign]->GetMean();
				// double Q_total_real_dEta_error1 = QvsdEta[27][sign]->GetMeanError();
				
				// double Q_total_real_dEta2 = QvsdEta[28][sign]->GetMean();
				// double Q_total_real_dEta_error2 = QvsdEta[28][sign]->GetMeanError();

				// double weight1 = delEta3p[sign]->GetBinContent( 28 );
				// double weight2 = delEta3p[sign]->GetBinContent( 29 );

				// double value = weightedAverage(weight1, weight2, 0, Q_total_real_dEta1, Q_total_real_dEta2, 0);
				// double error = weightedAverageError(weight1, weight2, 0, Q_total_real_dEta_error1, Q_total_real_dEta_error2, 0 );
				// cout << "value: " << value << endl;
				// cout << "error: " << error << endl;
				// hist1[sign]->SetBinContent(deta+1, value );
				// hist1[sign]->SetBinError(deta+1,  error);
			}
		}
	}
	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8,kBlack);
	TH1D* base2 = makeHist("base2", "", "#Delta#eta", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 48,0,4.8,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0012, 0.0013);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.0012, 0.0013);
	base2->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetXaxis()->SetTitleOffset(0.95);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetXaxis()->SetTitleOffset(0.95);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.3);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT (opposite - same)");
	base3->GetYaxis()->SetRangeUser(-0.0007,0.0018);
	base3->GetYaxis()->SetTitleOffset(1.2);
	base3->GetXaxis()->SetTitleOffset(1.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.2);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.2);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.3);
	base3->GetXaxis()->SetNdivisions(5,6,0);
	
	TH1D* temp11 = (TH1D*)hist1[0]->Clone("temp11");
	temp11->Add(hist1[1], +1);
	temp11->Scale(0.5);
	temp11->Scale(1.0/v2[3]);
	temp11->SetMarkerStyle(24);
	temp11->SetMarkerSize(1.4);
	temp11->SetMarkerColor(kRed);
	temp11->SetLineColor(kRed);

	TH1D* temp22 = (TH1D*) hist1[2]->Clone("temp22");
	temp22->SetMarkerStyle(25);
	temp22->Scale(1.0/v2[3]);
	temp22->SetMarkerColor(kBlue);
	temp22->SetMarkerSize(1.4);
	temp22->SetLineColor(kBlue);

	TCanvas* c2 = new TCanvas("c2","c2",1200,500);
	c2->Divide(3,1,0,0);
	c2->cd(1);
	gPad->SetLeftMargin(0.2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.07);
	gPad->SetTicks();
	base1->Draw();

	temp1->SetMarkerStyle(24);
	temp1->SetMarkerSize(1.4);
	temp1->SetMarkerColor(kRed);
	temp1->SetLineColor(kRed);

	temp2->SetMarkerColor(kBlue);
	temp2->SetMarkerSize(1.4);
	temp2->SetLineColor(kBlue);

	temp3->SetMarkerSize(1.4);
	temp3->SetMarkerStyle(24);
	temp3->SetMarkerColor(kRed);
	temp3->SetLineColor(kRed);
	
	temp4->SetMarkerSize(1.4);
	temp4->SetMarkerStyle(25);	
	temp4->SetMarkerColor(kBlue);
	temp4->SetLineColor(kBlue);

	temp1->Draw("Psame");
	temp2->Draw("Psame");

    TLatex* r3 = new TLatex(0.25, 0.85, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r3->SetNDC();
    r3->SetTextSize(23);
    r3->SetTextFont(43);
    r3->SetTextColor(kBlack);
    r3->Draw("same");

    TLatex* lmult = new TLatex(0.25, 0.76, "185 #leq N^{offline}_{trk} < 260");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

    double value1[50];
    double value1_error[50];
    double value2[50];
    double value2_error[50];
    double value3[50];
    double value3_error[50];
    double value4[50];
    double value4_error[50];
    double value5[50];
    double value5_error[50];
    double value6[50];
    double value6_error[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	value1[deta] = temp1->GetBinContent(deta+1);
    	value1_error[deta] = temp1->GetBinError(deta+1);

    	value2[deta] = temp2->GetBinContent(deta+1);
    	value2_error[deta] = temp2->GetBinError(deta+1);

    	value3[deta] = temp3->GetBinContent(deta+1);
    	value3_error[deta] = temp3->GetBinError(deta+1);

    	value4[deta] = temp4->GetBinContent(deta+1);
    	value4_error[deta] = temp4->GetBinError(deta+1);

    	value5[deta] = temp11->GetBinContent(deta+1);
    	value5_error[deta] = temp11->GetBinError(deta+1);

    	value6[deta] = temp22->GetBinContent(deta+1);
    	value6_error[deta] = temp22->GetBinError(deta+1);
    }

    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.02;
    	double ye = total_systematics_pPb;

    	box1[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value1[deta]-ye,dEtaReBinCenter2[deta]+xe,value1[deta]+ye);
		box1[deta]->SetFillColor(kRed);
        box1[deta]->SetFillStyle(0);
    	box1[deta]->SetLineWidth(1);
    	box1[deta]->SetLineColor(kRed);
        box1[deta]->Draw("SAME");

		box2[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value2[deta]-ye,dEtaReBinCenter2[deta]+xe,value2[deta]+ye);
		box2[deta]->SetFillColor(kBlue);
        box2[deta]->SetFillStyle(0);
    	box2[deta]->SetLineWidth(1);
    	box2[deta]->SetLineColor(kBlue);
        box2[deta]->Draw("SAME");
    }

    TLatex* r6 = new TLatex(0.73, 0.2, "Pb-going");
    r6->SetNDC();
    r6->SetTextSize(23);
    r6->SetTextFont(43);
    r6->SetTextColor(kBlack);
    r6->Draw("same");

//PAD1:
	c2->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.07);
	gPad->SetTicks();
	base2->Draw();

	temp3->Draw("Psame");
	temp4->Draw("Psame");

    TBox *box3[50];
    TBox *box4[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.02;
    	double ye = total_systematics_pPb;

    	box3[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value3[deta]-ye,dEtaReBinCenter2[deta]+xe,value3[deta]+ye);
		box3[deta]->SetFillColor(kRed);
        box3[deta]->SetFillStyle(0);
    	box3[deta]->SetLineWidth(1);
    	box3[deta]->SetLineColor(kRed);
        box3[deta]->Draw("SAME");

		box4[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value4[deta]-ye,dEtaReBinCenter2[deta]+xe,value4[deta]+ye);
		box4[deta]->SetFillColor(kBlue);
        box4[deta]->SetFillStyle(0);
    	box4[deta]->SetLineWidth(1);
    	box4[deta]->SetLineColor(kBlue);
        box4[deta]->Draw("SAME");
    }

    TLatex* r32 = new TLatex(0.11, 0.84, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r32->SetNDC();
    r32->SetTextSize(23);
    r32->SetTextFont(43);
    r32->SetTextColor(kBlack);
    //r32->Draw("same");

    TLatex* r62 = new TLatex(0.73, 0.2, "p-going");
    r62->SetNDC();
    r62->SetTextSize(23);
    r62->SetTextFont(43);
    r62->SetTextColor(kBlack);
    r62->Draw("same");

	c2->cd(3);
	gPad->SetBottomMargin(0.13);
	gPad->SetTicks();
	gPad->SetTopMargin(0.07);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	base2->Draw();

    TBox *box1[50];
    TBox *box2[50];

    for(int deta = 0; deta < NdEtaReBins2; deta++){

    	double xe = 0.02;
    	double ye = total_systematics_pPb;

    	box1[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value5[deta]-ye,dEtaReBinCenter2[deta]+xe,value5[deta]+ye);
		box1[deta]->SetFillColor(kRed);
        box1[deta]->SetFillStyle(0);
    	box1[deta]->SetLineWidth(1);
    	box1[deta]->SetLineColor(kRed);
        box1[deta]->Draw("SAME");

		box2[deta] = new TBox(dEtaReBinCenter2[deta]-xe,value6[deta]-ye,dEtaReBinCenter2[deta]+xe,value6[deta]+ye);
		box2[deta]->SetFillColor(kBlue);
        box2[deta]->SetFillStyle(0);
    	box2[deta]->SetLineWidth(1);
    	box2[deta]->SetLineColor(kBlue);
        box2[deta]->Draw("SAME");
    }

	temp11->Draw("Psame");
	temp22->Draw("Psame");

	TLegend *w2 = new TLegend(0.60,0.75,0.75,0.9);
    w2->SetLineColor(kWhite);
    w2->SetFillColor(0);
    w2->SetTextSize(20);
    w2->SetTextFont(45);
    w2->AddEntry(temp5, "  same");
    w2->AddEntry(temp9, "  opposite");
    w2->Draw("same");

    TLatex* r63 = new TLatex(0.73, 0.2, "Tracker");
    r63->SetNDC();
    r63->SetTextSize(23);
    r63->SetTextFont(43);
    r63->SetTextColor(kBlack);
    r63->Draw("same");

   	TLatex* r1 = new TLatex(0.58,0.95, "CMS");
    r1->SetNDC();
    r1->SetTextSize(0.06);
    r1->Draw("same");

    TLatex* r2 = new TLatex(0.72,0.95, "Preliminary");
    r2->SetNDC();
    r2->SetTextSize(22);
    r2->SetTextFont(53);
    r2->Draw("same");

    TCanvas* c3 = new TCanvas("c3","c3", 600, 600);
 	gPad->SetLeftMargin(0.16);
 	gPad->SetBottomMargin(0.13);
 	gPad->SetTopMargin(0.1);
	gPad->SetTicks();
	base3->Draw();

	TH1D* diff1 = (TH1D*) temp2->Clone("diff1");
	diff1->SetMarkerStyle(20);
	diff1->SetMarkerSize(1.3);
	diff1->SetMarkerColor(kRed);
	diff1->SetLineColor(kRed);
	diff1->Add(temp1, -1);
	
	TH1D* diff2 = (TH1D*) temp4->Clone("diff2");
	diff2->SetMarkerStyle(21);
	diff2->SetMarkerSize(1.3);
	diff2->SetMarkerColor(kBlue);
	diff2->SetLineColor(kBlue);
	diff2->Add(temp3, -1);

	diff1->Draw("Psame");
	diff2->Draw("Psame");

    TLatex* r33 = new TLatex(0.25, 0.82, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r33->SetNDC();
    r33->SetTextSize(23);
    r33->SetTextFont(43);
    r33->SetTextColor(kBlack);
    //r33->Draw("same");

    TLatex* lmult = new TLatex(0.20, 0.82, "185 #leq N^{offline}_{trk} < 260");
    lmult->SetNDC();
    lmult->SetTextSize(26);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);
    lmult->Draw("same");

	TH1D* diff3 = (TH1D*) temp9->Clone("diff3");
	diff3->Add(temp5, -1);
	diff3->SetMarkerStyle(28);
	diff3->SetMarkerSize(1.4);
	diff3->SetMarkerColor(kBlack);
	diff3->SetLineColor(kBlack);
	diff3->Draw("Psame");

	TH1D* diff4 = (TH1D*) temp22->Clone("diff4");
	diff4->Add(temp11, -1);
	diff4->SetMarkerStyle(34);
	diff4->SetMarkerSize(1.6);
	diff4->SetMarkerColor(kGreen+3);
	diff4->Draw("Psame");

   	TLatex* r11 = new TLatex(0.62,0.91, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.71,0.91, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

	TLegend *w3 = new TLegend(0.40,0.5,0.6,0.7);
    w3->SetLineColor(kWhite);
    w3->SetFillColor(0);
    w3->SetTextSize(20);
    w3->SetTextFont(43);
    w3->AddEntry(diff1, "pPb #sqrt{s_{NN}} = 5.02 TeV, Pb-going");
    w3->AddEntry(diff2, "pPb #sqrt{s_{NN}} = 5.02 TeV, p-going");
    w3->AddEntry(diff4, "pPb #sqrt{s_{NN}} = 5.02 TeV, Tracker");
    w3->AddEntry(diff3, "PbPb #sqrt{s_{NN}} = 2.76 TeV");
    w3->Draw("same");










}