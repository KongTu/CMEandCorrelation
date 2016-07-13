#include "RiceStyle.h"

using namespace std;

void plotDetaBinDist(){

	TFile* file[10];
	

	file[0] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v22_4.root");
	file[1] = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v27_7.root");

	TH1D* QvsdEta[30][48][3][2];

	TH1D* delEta3p[30][3][2];

	TH1D* Ntrk[10];

	for(int mult = 0; mult < 2; mult++){

		Ntrk[mult] = (TH1D*) file[mult]->Get(Form("ana/Ntrk"));

		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				delEta3p[mult][sign][HF] = (TH1D*) file[mult]->Get(Form("ana/delEta3p_%d_%d",sign,HF));
			}
		}
	}



	delEta3p[0][1][0]->Scale( 1.0/delEta3p[0][1][0]->Integral() );
	delEta3p[0][1][0]->SetMarkerStyle(20);
	delEta3p[0][1][0]->SetMarkerColor(kBlue);
	delEta3p[0][1][0]->GetYaxis()->SetRangeUser(0,0.2);
	delEta3p[0][1][0]->Draw("P");

	delEta3p[1][0][0]->Scale( 1.0/delEta3p[1][0][0]->Integral() );
	delEta3p[1][0][0]->SetMarkerStyle(20);
	delEta3p[1][0][0]->SetMarkerColor(kRed);
	delEta3p[1][0][0]->Draw("Psame");


	string etarange1 = "|#eta_{#alpha,#beta}| < 2.4";
	string etarange2 = "|#eta_{#alpha,#beta}| < 0.8";


	TLegend *w4 = new TLegend(0.4,0.6,0.8,0.75);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(18);
    w4->SetTextFont(45);
    w4->AddEntry(delEta3p[0][0][0], etarange1.c_str());
    w4->AddEntry(delEta3p[1][0][0], etarange2.c_str());
    w4->Draw("same");

}