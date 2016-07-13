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
double PbPb_5TeV_ntrkBinCenter[] = {103.921, 134.061, 166.539, 201.615, 239.096, 279.13, 324.044, 374.081, 448.17};

double PbPb_ntrkCentralityBinCenter[] = {625.275, 811.118, 1025.79, 1257.64};

double PbPb_centralityBinCenter[] = {57.5, 52.5, 47.5, 42.5, 37.5, 32.5};

const int Nmults = 5;

double total_systematics_pPb = 0.0000;
double total_systematics_PbPb = 0.0000;

double average_v2_PbPb = 9.82293235240882662e-02;
double average_v2_pPb = 6.87095335452075351e-02;

void plotIntegrated_v2(){

	gStyle->SetErrorX(0);

	TFile* file1[7];

    file1[0] = new TFile("../dataPoints/pPb_data.root");
	file1[1] = new TFile("../dataPoints/PbPb5TeV_data.root");

	TGraphErrors* gr1[6];
    TGraphErrors* gr2[5];

	for(int i = 0; i < 6; i++){

		gr1[i] = (TGraphErrors*) file1[0]->Get(Form("Graph;%d", i+1));

	}
    for(int i = 0; i < 5; i++){

        gr2[i] = (TGraphErrors*) file1[1]->Get(Form("Graph;%d", i+1));

    }

    for(int i = 0; i < 6; i++){
        for(int mult = 0; mult < 3; mult++){
            double x, y;
            gr1[i]->GetPoint(mult, x, y);
            gr1[i]->SetPoint(mult, x, 100);
        }

        gr1[i]->SetPoint(9,x,100);
    }

//start plotting

    TGaxis::SetMaxDigits(3);

    TH1D* base1 = makeHist("base1", "Pb-going", "N^{offline}_{trk}", "#sqrt{#LTQ_{2,HF#pm}Q^{*}_{2,HF#mp}#GT#LTQ_{2,HF#pm}Q^{*}_{2,trk}#GT/#LTQ_{2,HF#mp}Q^{*}_{2,trk}#GT}", 10000,0.1,600,kBlack);
    base1->GetYaxis()->SetRangeUser(0.01, 0.15);
    //base1->GetYaxis()->SetRangeUser(-0.0003, 0.0002);
    //base1->GetXaxis()->SetRangeUser(0.8, 1.4);
    base1->GetXaxis()->SetTitleColor(kBlack);
    
    fixedFontHist1D(base1,1.1,1.25);

    base1->GetYaxis()->SetTitleOffset(1.5);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.6);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.6);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.6);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.6);
    base1->GetXaxis()->SetNdivisions(4,6,0);
    base1->GetYaxis()->SetNdivisions(4,6,0);

    //TCanvas* c1 = new TCanvas("c1","c1",1,1,1300,650);
    //c1->Divide(2,1,0.01,0.01);
    //c1->cd(1);
    TCanvas* c1 = new TCanvas("c1","c1",1,1,700,700);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.20);
    gPad->SetBottomMargin(0.13);
    gStyle->SetPadBorderMode(0.1);
    gStyle->SetOptTitle(0);
    //gPad->SetLogx(1);

    base1->Draw();

    gr1[4]->SetMarkerStyle(20);
    gr1[4]->SetMarkerSize(1.4);
    gr1[4]->SetMarkerColor(kRed);
    gr1[4]->SetLineColor(kRed);
    gr1[4]->Draw("Psame");

    gr1[5]->SetMarkerStyle(24);
    gr1[5]->SetMarkerSize(1.4);
    gr1[5]->SetMarkerColor(kRed);
    gr1[5]->SetLineColor(kRed);
    gr1[5]->Draw("Psame");

    gr2[3]->SetMarkerStyle(21);
    gr2[3]->SetMarkerSize(1.4);
    gr2[3]->SetMarkerColor(kBlue);
    gr2[3]->SetLineColor(kRed);
    gr2[3]->Draw("Psame");

    gr2[4]->SetMarkerStyle(25);
    gr2[4]->SetMarkerSize(1.4);
    gr2[4]->SetMarkerColor(kBlue);
    gr2[4]->SetLineColor(kBlue);
    //gr2[4]->Draw("Psame");


    TLegend *w1 = new TLegend(0.5,0.6,0.7,0.8);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(18);
    w1->SetTextFont(43);
    w1->AddEntry(gr1[4], "#phi_{c}(Pb-going)  ", "P");
    w1->AddEntry(gr1[5], "#phi_{c}(p-going)", "P");
    w1->AddEntry(gr2[3], "PbPb", "P");
    w1->Draw("same");

    TLatex* r4 = new TLatex(0.23, 0.78, "#sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);

    TLatex* lmult = new TLatex(0.18, 0.72, "185 #leq N^{offline}_{trk} < 260");
    //TLatex* lmult = new TLatex(0.18, 0.72, "30-40%");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);

    TLatex* latex1 = new TLatex(0.29, 0.25, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    //latex1->Draw("same");
    TLatex* latex2 = new TLatex(0.36, 0.25, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);
    //latex2->Draw("same");


    TLatex* r11 = new TLatex(0.23,0.84, "CMS");
    r11->SetNDC();
    r11->SetTextSize(0.04);

    TLatex* r22 = new TLatex(0.32,0.84, "Preliminary");
    r22->SetNDC();
    r22->SetTextSize(24);
    r22->SetTextFont(53);
    
    r4->Draw("same");
    //lmult->Draw("same");
    r11->Draw("same");
    r22->Draw("same");


}