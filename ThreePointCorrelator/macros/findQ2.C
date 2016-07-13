#include "RiceStyle.h"

using namespace std;

void findQ2(){


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_Q2_v7_1.root");
	TH2D* q2_tracker_HF = (TH2D*) file->Get("ana/q2_tracker_HF");
	TH1D* q2_mag = (TH1D*) q2_tracker_HF->ProjectionY("q2_mag", 1,1000);
	
	q2_mag->Draw();

	double total = q2_mag->Integral();

	cout << "total: " << total << endl;
	
	for(int i = 0; i < q2_mag->GetNbinsX(); i++){

		double num = q2_mag->Integral(i+1, q2_mag->GetNbinsX());

		double ratio = num/total;

		if( ratio < 0.01 && ratio > 0.0095){

			cout << "bin number: " << i+1 << endl;
			cout << "value : " << q2_mag->Integral(i+1, 1000) << endl;
		}
	}

	return;


	TCanvas* c1 = new TCanvas("c1","c1",1,1,700,700);
	q2_mag->GetXaxis()->SetRangeUser(0,0.5);
	q2_mag->SetStats(kFALSE);


	q2_mag->Draw();

	gPad->SetTicks();
	gPad->SetLogy(1);

	TLine* l1[10];
	double temp = q2_mag->GetBinContent(30);
    l1[0] = new TLine(0.030,0, 0.030, temp);

    double temp = q2_mag->GetBinContent(70);
    l1[1] = new TLine(0.070,0, 0.070, temp);

    double temp = q2_mag->GetBinContent(92);
    l1[2] = new TLine(0.092,0, 0.092, temp);

    double temp = q2_mag->GetBinContent(108);
    l1[3] = new TLine(0.108,0, 0.108, temp);

    double temp = q2_mag->GetBinContent(124);
    l1[4] = new TLine(0.124,0, 0.124, temp);

    double temp = q2_mag->GetBinContent(143);
    l1[5] = new TLine(0.143,0, 0.143, temp);

    double temp = q2_mag->GetBinContent(165);
    l1[6] = new TLine(0.165,0, 0.165, temp);

    double temp = q2_mag->GetBinContent(198);
    l1[7] = new TLine(0.198,0, 0.198, temp);

    double temp = q2_mag->GetBinContent(226);
    l1[8] = new TLine(0.226,0, 0.226, temp);

    double temp = q2_mag->GetBinContent(280);
    l1[9] = new TLine(0.280,0, 0.280, temp);

   	TLatex* r1 = new TLatex(0.55,0.84, "Q2 ESE classes:");
    r1->SetNDC();
    r1->SetTextSize(0.04);
    r1->Draw("same");

	TLatex* s1[11];
	s1[0] = new TLatex(0.12,0.2, "1");
	s1[1] = new TLatex(0.18,0.2, "2");
	s1[2] = new TLatex(0.22,0.2, "3");
	s1[3] = new TLatex(0.25,0.2, "4");
	s1[4] = new TLatex(0.28,0.2, "5");
	s1[5] = new TLatex(0.31,0.2, "6");
	s1[6] = new TLatex(0.34,0.2, "7");
	s1[7] = new TLatex(0.39,0.2, "8");
	s1[8] = new TLatex(0.44,0.2, "9");
	s1[9] = new TLatex(0.50,0.2, "10");
	s1[10] = new TLatex(0.60,0.2, "11");
    s1[10]->SetNDC();
	s1[10]->SetTextSize(0.04);
	s1[10]->Draw("same");

	TLatex* s2[11];
	s2[0] = new TLatex(0.70,0.78, "1, 95-100%");
	s2[1] = new TLatex(0.70,0.75, "2, 80-95%");
	s2[2] = new TLatex(0.70,0.72, "3, 60-80%");
	s2[3] = new TLatex(0.70,0.69, "4, 50-60%");
	s2[4] = new TLatex(0.70,0.66, "5, 40-50%");
	s2[5] = new TLatex(0.70,0.63, "6, 30-40%");
	s2[6] = new TLatex(0.70,0.60, "7, 20-30%");
	s2[7] = new TLatex(0.70,0.57, "8, 10-20%");
	s2[8] = new TLatex(0.70,0.54, "9, 5-10%");
	s2[9] = new TLatex(0.70,0.51, "10, 1-5%");
	s2[10] = new TLatex(0.70,0.48, "11, 0-1%");
	s2[10]->SetNDC();
	s2[10]->SetTextSize(0.03);
	s2[10]->Draw("same");

    for(int i = 0; i < 10; i++){
    	l1[i]->SetLineWidth(1.3);
    	l1[i]->SetLineStyle(2);
    	l1[i]->SetLineColor(kRed);
    	l1[i]->Draw("Lsame");

		s1[i]->SetNDC();
		s1[i]->SetTextSize(0.04);
		s1[i]->Draw("same");

		s2[i]->SetNDC();
		s2[i]->SetTextSize(0.03);
		s2[i]->Draw("same");
    }

	double total = q2_mag->Integral();
	
	for(int i = 0; i < q2_mag->GetNbinsX(); i++){

		double num = q2_mag->Integral(i+1, q2_mag->GetNbinsX());

		double ratio = num/total;

		if( ratio < 0.97 && ratio > 0.965){

			cout << "bin number: " << i+1 << endl;
		}
	}

	TLatex* r4 = new TLatex(0.18, 0.82, "PbPb #sqrt{s_{NN}} = 5.02 TeV");
    r4->SetNDC();
    r4->SetTextSize(23);
    r4->SetTextFont(43);
    r4->SetTextColor(kBlack);

    TLatex* lmult = new TLatex(0.18, 0.76, "185 #leq N^{offline}_{trk} < 260");
    lmult->SetNDC();
    lmult->SetTextSize(23);
    lmult->SetTextFont(43);
    lmult->SetTextColor(kBlack);

	TCanvas*c2 = new TCanvas("c2","c2",1,1,700,700);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	//double v2_alpha[] = {0.0615888, 0.0632844, 0.0645851, 0.0664797, 0.0683457, 0.069671, 0.0722178, 0.0748879, 0.0782873, 0.0823528, 0.085673};
	//double v2_c_Pb[] = {0.011399, 0.0174723, 0.0291536, 0.0399171, 0.0445492, 0.0542639, 0.0687932, 0.0833151, 0.102611, 0.129255, 0.148249};
	double v2_c_p[] = {0.0269153, 0.0269532, 0.0287159, 0.0308912, 0.0294902, 0.0304618, 0.0324812, 0.0330777, 0.0345458, 0.036932, 0.0385128};
	
	double v2_alpha[] = {0.0850408, 0.0873448, 0.0901072, 0.0924002, 0.0939832, 0.0964905, 0.0989382, 0.10316, 0.10682, 0.111734, 0.118806};
	double v2_c_Pb[] = {0.0162113, 0.0241807, 0.0351651, 0.0434707, 0.0518347, 0.0600376, 0.0713306, 0.0865037, 0.102095, 0.122752, 0.158914};
	
	double xwidth[] = {0,0,0,0,0,0,0,0,0,0,0};
	TGraphErrors* gr1 = new TGraphErrors(11, v2_alpha, v2_c_Pb, xwidth, xwidth);
	TGraphErrors* gr2 = new TGraphErrors(11, v2_alpha, v2_c_p, xwidth, xwidth);

	TH1D* base1 = makeHist("base1","","v_{2,tracker}", "v_{2,c}", 1000,0.0,0.2,kBlack);
	base1->GetYaxis()->SetRangeUser(0.0, 0.2);
	base1->GetXaxis()->SetRangeUser(0.08, 0.13);
	base1->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base1,1.1,1.25);

	base1->GetYaxis()->SetTitleOffset(1.3);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.6);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.6);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.6);
	base1->GetXaxis()->SetNdivisions(4,6,0);
	base1->GetYaxis()->SetNdivisions(4,6,0);
	base1->Draw();
	

	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(1.4);
	gr1->SetMarkerColor(kRed);
	gr1->SetLineColor(kRed);
	gr1->Draw("Psame");

	gr1->Fit("pol1");
    TF1 * myFunc1 = gr1->GetFunction("pol1");
    myFunc1->SetLineStyle(2);
    double intersect_1 = myFunc1->GetParameter(0);
    double intersect_1_error = myFunc1->GetParError(0);
    double slope_1 = myFunc1->GetParameter(1);
    double slope_1_error = myFunc1->GetParError(1);

    TLatex* latex3 = new TLatex(0.18, 0.66, Form("slope: %.4f +/- %.4f",slope_1, slope_1_error ));
    latex3->SetNDC();
    latex3->SetTextSize(20);
    latex3->SetTextFont(43);
    latex3->SetTextColor(kRed);
    latex3->Draw("same");

    TLatex* latex4 = new TLatex(0.18, 0.63, Form("intersect: %.4f +/- %.4f",intersect_1, intersect_1_error ));
    latex4->SetNDC();
    latex4->SetTextSize(20);
    latex4->SetTextFont(43);
    latex4->SetTextColor(kRed);
    latex4->Draw("same");

	// gr2->SetMarkerStyle(24);
	// gr2->SetMarkerSize(1.4);
	// gr2->SetMarkerColor(kBlue);
	// gr2->SetLineColor(kBlue);
	// gr2->Draw("Psame");


 //    gr2->Fit("pol1");

 //    TF1 * myFunc2 = gr2->GetFunction("pol1");
 //    myFunc2->SetLineColor(kBlue);
 //    myFunc2->SetLineStyle(2);
 //    double intersect_2 = myFunc2->GetParameter(0);
 //    double intersect_2_error = myFunc2->GetParError(0);
 //    double slope_2 = myFunc2->GetParameter(1);
 //    double slope_2_error = myFunc2->GetParError(1);

 //    TLatex* latex5 = new TLatex(0.18, 0.59, Form("slope: %.4f +/- %.4f",slope_2, slope_2_error ));
 //    latex5->SetNDC();
 //    latex5->SetTextSize(20);
 //    latex5->SetTextFont(43);
 //    latex5->SetTextColor(kBlue);
 //    latex5->Draw("same");

 //    TLatex* latex6 = new TLatex(0.18, 0.56, Form("intersect: %.4f +/- %.4f",intersect_2, intersect_2_error ));
 //    latex6->SetNDC();
 //    latex6->SetTextSize(20);
 //    latex6->SetTextFont(43);
 //    latex6->SetTextColor(kBlue);
 //    latex6->Draw("same");
	
	r4->Draw("same");
    lmult->Draw("same");

	// TFile* file = new TFile("../rootfiles/CME_QvsdEta_EPOS_PbPb_v2.root");
	// TH1D* Ntrk = (TH1D*) file->Get("ana/Ntrk");

	// double total = Ntrk->Integral();
	
	// for(int i = 0; i < Ntrk->GetNbinsX(); i++){

	// 	double num = Ntrk->Integral(i+1, Ntrk->GetNbinsX());

	// 	double ratio = num/total;

	// 	if( ratio < 0.8 && ratio > 0.799){

	// 		cout << "bin number: " << i+1 << endl;
	// 	}
	// }


}