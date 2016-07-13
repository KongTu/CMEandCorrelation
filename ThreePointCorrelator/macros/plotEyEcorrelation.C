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

double xbinwidth[] = {0.0,0.0,0.0,0.0,0.0};
double pPb_v2[] = {-0.15, -0.05, 0.0, 0.05, 0.1, 0.15};


const int Nmults = 5;

double total_systematics_pPb = 0.00006;
double total_systematics_PbPb = 0.000051;

double binNumber(double i){

	return (i - (-1.0) )/0.001;
}

void plotEyEcorrelation(){

	gStyle->SetErrorX(0);
	TGaxis::SetMaxDigits(3);

	TFile* file[2];

	file[0] = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v35.root");
	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v18.root");

	TH2D* QvsV2_pPb[3][2];
	TH2D* QvsV2_PbPb[3][2];

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF++){

			QvsV2_pPb[sign][HF] = (TH2D*) file[0]->Get(Form("ana/QvsV2_%d_%d",sign,HF));
			QvsV2_PbPb[sign][HF] = (TH2D*) file[1]->Get(Form("ana/QvsV2_%d_%d",sign,HF));

		}
	}

	for(int HF = 0; HF < 2; HF++){

		QvsV2_pPb[0][HF]->Add(QvsV2_pPb[1][HF], +1);
		QvsV2_pPb[0][HF]->Scale(0.5);

		QvsV2_PbPb[0][HF]->Add(QvsV2_PbPb[1][HF], +1);
		QvsV2_PbPb[0][HF]->Scale(0.5);

	}

	for(int sign = 0; sign < 3; sign++){

		QvsV2_PbPb[sign][0]->Add(QvsV2_PbPb[sign][1],+1);
		QvsV2_PbPb[sign][0]->Scale(0.5);

	}
	
	double pPb_v2_bincenter[10];
	
	double pPb_3p_value[3][2][10];
	double pPb_3p_value_error[3][2][10];

	double pPb_3p_diff[2][10];
	double pPb_3p_diff_error[2][10];

	double PbPb_3p_value[3][2][10];
	double PbPb_3p_value_error[3][2][10];

	double PbPb_3p_diff[2][10];
	double PbPb_3p_diff_error[2][10];

	for(int i = 0; i < Nmults; i++){
		for(int HF = 0; HF < 2; HF++){
			for(int sign = 0; sign < 3; sign++){

				pPb_v2_bincenter[i] = (pPb_v2[i+1] - pPb_v2[i])/2.0 + pPb_v2[i];
				pPb_3p_value[sign][HF][i] = QvsV2_pPb[sign][HF]->ProjectionY(Form("test_%d",i), binNumber(pPb_v2[i]), binNumber(pPb_v2[i+1])-1 )->GetMean();
				pPb_3p_value_error[sign][HF][i] = QvsV2_pPb[sign][HF]->ProjectionY(Form("test_%d",i), binNumber(pPb_v2[i]), binNumber(pPb_v2[i+1])-1 )->GetMeanError();
			
				PbPb_3p_value[sign][HF][i] = QvsV2_PbPb[sign][HF]->ProjectionY(Form("test_%d",i), binNumber(pPb_v2[i]), binNumber(pPb_v2[i+1])-1 )->GetMean();
				PbPb_3p_value_error[sign][HF][i] = QvsV2_PbPb[sign][HF]->ProjectionY(Form("test_%d",i), binNumber(pPb_v2[i]), binNumber(pPb_v2[i+1])-1 )->GetMeanError();
			}
		}
	}

    TLatex* r41 = new TLatex(0.24, 0.87, "pPb #sqrt{s_{NN}} = 5.02 TeV");
    r41->SetNDC();
    r41->SetTextSize(23);
    r41->SetTextFont(43);
    r41->SetTextColor(kBlack);

    TLatex* r42 = new TLatex(0.05, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.57,0.95, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.052);
    
    TLatex* r44 = new TLatex(0.7,0.95, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(24);
    r44->SetTextFont(53);
    
    TLatex* r45 = new TLatex(0.05, 0.80, "185 #leq N^{offline}_{trk} < 260");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r46 = new TLatex(0.24, 0.80, "185 #leq N^{offline}_{trk} < 260");
    r46->SetNDC();
    r46->SetTextSize(23);
    r46->SetTextFont(43);
    r46->SetTextColor(kBlack);

	TLatex* latex1 = new TLatex(0.48, 0.37, "same");
    latex1->SetNDC();
    latex1->SetTextSize(20);
    latex1->SetTextFont(43);
    latex1->SetTextColor(kBlack);
    TLatex* latex2 = new TLatex(0.59, 0.37, "oppo");
    latex2->SetNDC();
    latex2->SetTextSize(20);
    latex2->SetTextFont(43);
    latex2->SetTextColor(kBlack);

	TH1D* base1 = makeHist("base1", "", "#LT cos(2#phi_{#alpha}-2#phi_{c})#GT", "#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT", 1000,-0.25,0.25,kBlack);
	base1->GetYaxis()->SetRangeUser(-0.001,0.001);
	base1->GetXaxis()->SetRangeUser(-0.2, 4.8);
	base1->GetXaxis()->SetTitleColor(kBlack);

	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.5);
	base1->GetXaxis()->SetTitleOffset(1.0.);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.3);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.4);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.4);
	base1->GetXaxis()->SetNdivisions(8,18,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);

	TGraphErrors* gr1 = new TGraphErrors(Nmults, pPb_v2_bincenter, pPb_3p_value[0][0], xbinwidth, pPb_3p_value_error[0][0]);
	TGraphErrors* gr2 = new TGraphErrors(Nmults, pPb_v2_bincenter, pPb_3p_value[2][0], xbinwidth, pPb_3p_value_error[2][0]);
	
	TGraphErrors* gr3 = new TGraphErrors(Nmults, pPb_v2_bincenter, pPb_3p_value[0][1], xbinwidth, pPb_3p_value_error[0][1]);
	TGraphErrors* gr4 = new TGraphErrors(Nmults, pPb_v2_bincenter, pPb_3p_value[2][1], xbinwidth, pPb_3p_value_error[2][1]);

	TGraphErrors* gr5 = new TGraphErrors(Nmults, pPb_v2_bincenter, PbPb_3p_value[0][0], xbinwidth, PbPb_3p_value_error[0][0]);
	TGraphErrors* gr6 = new TGraphErrors(Nmults, pPb_v2_bincenter, PbPb_3p_value[2][0], xbinwidth, PbPb_3p_value_error[2][0]);

	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(kBlue);
	gr2->SetLineColor(kBlue);

	gr3->SetMarkerStyle(24);
	gr3->SetMarkerColor(kRed);
	gr3->SetLineColor(kRed);

	gr4->SetMarkerStyle(25);
	gr4->SetMarkerColor(kBlue);
	gr4->SetLineColor(kBlue);

	gr5->SetMarkerStyle(20);
	gr5->SetMarkerColor(kRed);
	gr5->SetLineColor(kRed);

	gr6->SetMarkerStyle(21);
	gr6->SetMarkerColor(kBlue);
	gr6->SetLineColor(kBlue);

    TLegend *w4 = new TLegend(0.5,0.2,0.95,0.35);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(23);
    w4->SetTextFont(45);
    w4->SetNColumns(2);
    w4->AddEntry(gr1, "  ", "P");
    w4->AddEntry(gr2, "  #phi_{c}(Pb-going)", "P");
    
    w4->AddEntry(gr3, "  ", "P");
    w4->AddEntry(gr4, "  #phi_{c}(p-going)", "P");

    TLegend *w40 = new TLegend(0.5,0.22,0.73,0.4);
    w40->SetLineColor(kWhite);
    w40->SetFillColor(0);
    w40->SetTextSize(23);
    w40->SetTextFont(45);
    w40->SetNColumns(2);
    w40->SetTextAlign();
    w40->AddEntry(gr5, " ", "P");
    w40->AddEntry(gr6, " ", "P");

    TLatex* r434 = new TLatex(0.70,0.29, "PbPb");
    r434->SetNDC();
    r434->SetTextSize(21);
    r434->SetTextFont(43);

	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	c1->Divide(2,1,0,0);
	c1->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

	r41->Draw("Psame");
	r46->Draw("Psame");
	w4->Draw("Psame");
	latex1->Draw("Psame");
	latex2->Draw("Psame");

	gr1->Draw("Psame");
	gr2->Draw("Psame");
	gr3->Draw("Psame");
	gr4->Draw("Psame");

	c1->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base1->Draw();

	w40->Draw("Psame");
	r434->Draw("Psame");
	latex1->Draw("Psame");
	latex2->Draw("Psame");
	

	r42->Draw("Psame");
	r43->Draw("Psame");
	r44->Draw("Psame");
	r45->Draw("Psame");

	gr5->Draw("Psame");
	gr6->Draw("Psame");

	TGraphErrors* dif1 = new TGraphErrors(Nmults);
	TGraphErrors* dif2 = new TGraphErrors(Nmults);
	TGraphErrors* dif3 = new TGraphErrors(Nmults);

	for(int i = 0; i < Nmults; i++){

	//pPb, Pb-going
		double x1,y1,ey1;
		gr1->GetPoint(i, x1, y1);
		ey1 = gr1->GetErrorY(i);
		double x2, y2,ey2;
		gr2->GetPoint(i, x2, y2);
		ey2 = gr2->GetErrorY(i);

		double x, y, ey;
		x = x1;
		y = y2-y1;
		ey = sqrt(ey1*ey1 + ey2*ey2);

		dif1->SetPoint(i, x, y);
		dif1->SetPointError(i, 0, ey);
	
	//pPb, p-going
		double x1,y1,ey1;
		gr3->GetPoint(i, x1, y1);
		ey1 = gr3->GetErrorY(i);
		double x2, y2,ey2;
		gr4->GetPoint(i, x2, y2);
		ey2 = gr4->GetErrorY(i);

		double x, y, ey;
		x = x1;
		y = y2-y1;
		ey = sqrt(ey1*ey1 + ey2*ey2);

		dif2->SetPoint(i, x, y);
		dif2->SetPointError(i, 0, ey);

	//PbPb
		double x1,y1,ey1;
		gr5->GetPoint(i, x1, y1);
		ey1 = gr5->GetErrorY(i);
		double x2, y2,ey2;
		gr6->GetPoint(i, x2, y2);
		ey2 = gr6->GetErrorY(i);

		double x, y, ey;
		x = x1;
		y = y2-y1;
		ey = sqrt(ey1*ey1 + ey2*ey2);

		dif3->SetPoint(i, x, y);
		dif3->SetPointError(i, 0, ey);
	}


	TH1D* base2 = (TH1D*) base1->Clone("base2");
	base2->GetYaxis()->SetRangeUser(-0.0006,0.001);
	base2->GetYaxis()->SetTitle("#LTcos(#phi_{#alpha}+#phi_{#beta}-2#phi_{c})#GT(oppo-same)");


	TLegend *w5 = new TLegend(0.5,0.2,0.95,0.35);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(23);
    w5->SetTextFont(45);
    w5->AddEntry(dif1, "#phi_{c}(Pb-going)", "P");
    w5->AddEntry(dif2, "#phi_{c}(p-going)", "P");
    

	TCanvas* c2 = new TCanvas("c2","c2",1000,600);
	c2->Divide(2,1,0,0);
	c2->cd(1);
	gPad->SetLeftMargin(0.20);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetTicks();
	base2->Draw();
	
	dif1->SetMarkerStyle(20);
	dif1->SetMarkerColor(kRed);
	dif1->SetLineColor(kRed);

	dif2->SetMarkerStyle(24);
	dif2->SetMarkerColor(kBlue);
	dif2->SetLineColor(kBlue);

	dif3->SetMarkerStyle(20);
	dif3->SetMarkerColor(kRed);
	dif3->SetLineColor(kRed);

	dif1->Draw("Psame");
	dif2->Draw("Psame");

	r41->Draw("Psame");
	r46->Draw("Psame");

	w5->Draw("Psame");

	c2->cd(2);
	gPad->SetBottomMargin(0.13);
	gPad->SetTopMargin(0.06);
	gPad->SetRightMargin(0.06);
	gPad->SetTicks();
	base2->Draw();

	r43->Draw("Psame");
	r44->Draw("Psame");
	
	dif3->Draw("Psame");
	

	r42->Draw("Psame");
	r45->Draw("Psame");



}