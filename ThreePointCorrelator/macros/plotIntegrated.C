#include "RiceStyle.h"

using namespace std;

const int MAXTRACKS = 40000;
//double etabins[] = {-2.4,-2.35,-2.3,-2.25,-2.2,-2.15,-2.1,-2.05,-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7,-1.65,-1.6,-1.55,-1.5,-1.45,-1.4,-1.35,-1.3,-1.25,-1.2,-1.15,-1.1,-1.05,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;
double ntrkBins[] = {0,35,60,90,120,150,185,220,260};
const int NntrkBins = sizeof(ntrkBins) / sizeof(ntrkBins[0]) - 1;
int ntrkBinCenter[] = {17.5, 47.5, 75, 105, 135, 167.5, 202.5, 240};

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double pPb_ntrkBinCenter[] = {16.29,46.1,74.22,101.7,131.3,162.1,196.7,231.5};
double PbPb_ntrkBinCenter[] ={13.8,46.15,73.67,103.9,134,167,202,239.1};

double PbPb_ntrkCentralityBinCenter[] = {151.6, 270.2, 441.9, 685.4,1024,1376,1721};
const int Nmults = 23;

double total_systematics_pPb = 0.00015;
double total_systematics_PbPb = 0.00014;

void plotIntegrated(){

	gStyle->SetErrorX(0);

	TFile* file[3];
	file[0] = new TFile("../dataPoints/pPb_data.root");
	file[1] = new TFile("../dataPoints/PbPb_data.root");
	file[2] = new TFile("../dataPoints/PbPb_centrality_data.root");

	TGraphErrors* gr1[4];
	for(int i = 0; i < 4; i++){

		gr1[i] = (TGraphErrors*) file[0]->Get(Form("Graph;%d", i+1));
	}

	TGraphErrors* gr2[2];
	for(int i = 0; i < 2; i++){

		gr2[i] = (TGraphErrors*) file[1]->Get(Form("Graph;%d", i+1));
	}
	
	TGraphErrors* gr3[2];
	for(int i = 0; i < 2; i++){

		gr3[i] = (TGraphErrors*) file[2]->Get(Form("Graph;%d", i+1));
	}
	
	for(int i = 0; i < 4; i++){
		for(int mult = 0; mult < 3; mult++){
			double x, y;
			gr1[i]->GetPoint(mult, x, y);
			gr1[i]->SetPoint(mult, x, 100);
		}
	}
	for(int i = 0; i < 2; i++){
		for(int mult = 0; mult < 3; mult++){
			double x, y;
			gr2[i]->GetPoint(mult, x, y);
			gr2[i]->SetPoint(mult, x, 100);
		}
	}


//start plotting

	TGaxis::SetMaxDigits(3);

	TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 5000,0.1,10000,kBlack);
	TH1D* base2 = makeHist("base2", "p-going", "N^{offline}_{trk}", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#Psi_{RP})#GT", 5000,0.1,10000,kBlack);

	base1->GetYaxis()->SetRangeUser(-0.0013, 0.0007);
	base1->GetXaxis()->SetRangeUser(70, 4000);

	base1->GetXaxis()->SetTitleColor(kBlack);
	
	base2->GetYaxis()->SetRangeUser(-0.001, 0.015);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);
	fixedFontHist1D(base2,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.4);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);

	base2->GetYaxis()->SetTitleOffset(1.9);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.4);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.4);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.4);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.4);
	
	TH1D* base3 = (TH1D*) base1->Clone("base3");
	base3->GetYaxis()->SetRangeUser(-0.0006,0.0012);
	base3->GetYaxis()->SetTitleOffset(1.9);
	base3->GetXaxis()->SetTitleOffset(3.1);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.0);
	base3->GetYaxis()->SetNdivisions(6);
	
	TH1D* base4 = (TH1D*) base2->Clone("base4");
	base4->GetYaxis()->SetRangeUser(-0.0006,0.0012);
	base4->GetYaxis()->SetTitleOffset(1.9);
	base4->GetXaxis()->SetTitleOffset(3.1);
	base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.0);
	base4->GetYaxis()->SetNdivisions(6);

	TCanvas* c1 = new TCanvas("c1","c1",1,1,650,650);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);
	gPad->SetLogx(1);

	base1->Draw();

	gr1[0]->SetMarkerStyle(20);
	gr1[0]->SetMarkerSize(1.4);
	gr1[0]->SetMarkerColor(kRed);
	gr1[0]->SetLineColor(kRed);
	gr1[0]->Draw("Psame");

	gr1[1]->SetMarkerStyle(21);
	gr1[1]->SetMarkerSize(1.4);
	gr1[1]->SetMarkerColor(kBlue);
	gr1[1]->SetLineColor(kBlue);
	gr1[1]->Draw("Psame");

	// gr1[2]->SetMarkerStyle(20);
	// gr1[2]->SetMarkerColor(kRed);
	// gr1[2]->SetLineColor(kRed);
	// gr1[2]->Draw("Psame");

	// gr1[3]->SetMarkerStyle(21);
	// gr1[3]->SetMarkerColor(kBlue);
	// gr1[3]->SetLineColor(kBlue);
	// gr1[3]->Draw("Psame");

	gr2[0]->SetMarkerStyle(28);
	gr2[0]->SetMarkerSize(1.4);
	gr2[0]->SetMarkerColor(kRed);
	gr2[0]->SetLineColor(kRed);
	gr2[0]->Draw("Psame");

	gr2[1]->SetMarkerStyle(28);
	gr2[1]->SetMarkerSize(1.4);
	gr2[1]->SetMarkerColor(kBlue);
	gr2[1]->SetLineColor(kBlue);
	gr2[1]->Draw("Psame");

	gr3[0]->SetMarkerStyle(24);
	gr3[0]->SetMarkerSize(1.4);
	gr3[0]->SetMarkerColor(kRed);
	gr3[0]->SetLineColor(kRed);
	gr3[0]->Draw("Psame");

	gr3[1]->SetMarkerStyle(25);
	gr3[1]->SetMarkerSize(1.4);
	gr3[1]->SetMarkerColor(kBlue);
	gr3[1]->SetLineColor(kBlue);
	gr3[1]->Draw("Psame");


    TBox *box1[50];
    TBox *box2[50];
    TBox *box3[50];
    TBox *box4[50];

    for(int mult = 3; mult < 8; mult++){

    	double xe = 5;
    	double ye = total_systematics_pPb;

    	double x1;
    	double value1;
    	gr1[0]->GetPoint(mult, x1, value1);

    	double x2;
    	double value2;
    	gr1[1]->GetPoint(mult, x2, value2);


    	box1[mult] = new TBox(x1-xe,value1-ye,x1+xe,value1+ye);
		box1[mult]->SetFillColor(kRed);
        box1[mult]->SetFillStyle(0);
    	box1[mult]->SetLineWidth(1);
    	box1[mult]->SetLineColor(kRed);
        box1[mult]->Draw("SAME");

		box2[mult] = new TBox(x2-xe,value2-ye,x2+xe,value2+ye);
		box2[mult]->SetFillColor(kBlue);
        box2[mult]->SetFillStyle(0);
    	box2[mult]->SetLineWidth(1);
    	box2[mult]->SetLineColor(kBlue);
        box2[mult]->Draw("SAME");

  //   	box3[mult] = new TBox(pPb_ntrkBinCenter[mult]-xe,value3[mult]-ye,pPb_ntrkBinCenter[mult]+xe,value3[mult]+ye);
		// box3[mult]->SetFillColor(kRed);
  //       box3[mult]->SetFillStyle(0);
  //   	box3[mult]->SetLineWidth(1);
  //   	box3[mult]->SetLineColor(kRed);
  //       box3[mult]->Draw("SAME");

		// box4[mult] = new TBox(pPb_ntrkBinCenter[mult]-xe,value4[mult]-ye,pPb_ntrkBinCenter[mult]+xe,value4[mult]+ye);
		// box4[mult]->SetFillColor(kBlue);
  //       box4[mult]->SetFillStyle(0);
  //   	box4[mult]->SetLineWidth(1);
  //   	box4[mult]->SetLineColor(kBlue);
  //       box4[mult]->Draw("SAME");
    }

    //VERZO:
	double alice_centrality_oppo1[] = {1.08E-4, -2.6E-5, -6.4E-5, -2.5E-5, -1.4E-5, -1.0E-6, 0.0 };
	double alice_centrality_oppo1_error[] = {8.0E-5, 4.1E-5, 2.3E-5, 1.4E-5, 9.0E-6, 6.0E-6, 4.0E-6};

	double alice_centrality_same1[] = {-3.66E-4, -2.48E-4, -1.47E-4, -9.8E-5, -4.6E-5, -1.0E-5, 0.0};
	double alice_centrality_same1_error[] = {8.1E-5, 4.1E-5, 2.4E-5, 1.4E-5, 9.0E-6, 6.0E-6, 4.0E-6};

	
	double alice_centrality_oppo2[7];
	double alice_centrality_oppo2_error[7];
	
	double alice_centrality_same2[7];
	double alice_centrality_same2_error[7];

	//TPC cumulant
	double p8370_d2x1y1_yval[] = { -6.0E-6, 2.0E-6, -5.0E-6, -1.1E-5, -1.0E-5, -3.0E-5, -1.0E-5 };
	double p8370_d2x1y1_yerrminus[] = { 3.0E-6, 4.0E-6, 3.0E-6, 5.0E-6, 7.0E-6, 1.2E-5, 2.2E-5 };

	double p8370_d2x1y2_yval[] = { -1.0E-5, -2.3E-5, -5.4E-5, -9.7E-5, -1.6E-4, -2.69E-4, -3.72E-4 };
	double p8370_d2x1y2_yerrminus[] = { 3.0E-6, 4.0E-6, 3.0E-6, 5.0E-6, 7.0E-6, 0.0001, 0.0001 };

	for(int mult = 0; mult < 7; mult++){

		alice_centrality_oppo2[mult] = p8370_d2x1y1_yval[6-mult];
		alice_centrality_oppo2_error[mult] = p8370_d2x1y1_yerrminus[6-mult];

		alice_centrality_same2[mult] = p8370_d2x1y2_yval[6-mult];
		alice_centrality_same2_error[mult] = p8370_d2x1y2_yerrminus[6-mult];
	}


	TGraphErrors* gr4 = new TGraphErrors(7, PbPb_ntrkCentralityBinCenter, alice_centrality_same2, xbinwidth, alice_centrality_same2_error);
	TGraphErrors* gr5 = new TGraphErrors(7, PbPb_ntrkCentralityBinCenter, alice_centrality_oppo2, xbinwidth, alice_centrality_oppo2_error);


	gr4->SetMarkerStyle(27);
	gr4->SetMarkerSize(1.4);
	gr4->SetMarkerColor(kBlack);
	gr4->SetLineColor(kBlack);
	gr4->Draw("Psame");


	gr5->SetMarkerStyle(28);
	gr5->SetMarkerSize(1.4);
	gr5->SetMarkerColor(kBlack);
	gr5->SetLineColor(kBlack);
	gr5->Draw("Psame");
	//c2->Print("../results/IntegratedResults.pdf");

	TLegend *w1 = new TLegend(0.40,0.15,0.8,0.4);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(20);
    w1->SetTextFont(43);
    w1->AddEntry(gr1[0], "pPb Pb-going, like sign", "P");
    w1->AddEntry(gr1[1], "pPb Pb-going, unlike sign", "P");
    w1->AddEntry(gr3[0], "CMS PbPb, like sign", "P");
    w1->AddEntry(gr3[1], "CMS PbPb, unlike sign", "P");
    w1->AddEntry(gr4, "ALICE PbPb, like sign", "P");
    w1->AddEntry(gr5, "ALICE PbPb, unlike sign", "P");
    w1->Draw("same");

   	TLatex* r11 = new TLatex(0.62,0.91, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);
    r11->Draw("same");

    TLatex* r22 = new TLatex(0.71,0.91, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    r22->Draw("same");

	return;

}