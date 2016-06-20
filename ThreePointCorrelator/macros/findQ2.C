#include "RiceStyle.h"

using namespace std;

void findQ2(){


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_PbPb_5TeV_30_100_v9.root");
	TH2D* q2_tracker_HF = (TH2D*) file->Get("ana/q2_tracker_HF");
	TH1D* q2_mag = (TH1D*) q2_tracker_HF->ProjectionY("q2_mag", 1,1000);

	double total = q2_mag->Integral();
	for(int i = 0; i < q2_mag->GetNbinsX(); i++){

		double num = q2_mag->Integral(i, q2_mag->GetNbinsX());

		double ratio = num/total;

		if( ratio < 0.95 && ratio > 0.945){

			cout << "bin number: " << i << endl;
		}
	}


}