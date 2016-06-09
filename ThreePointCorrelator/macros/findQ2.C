#include "RiceStyle.h"

using namespace std;

void findQ2(){


	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_Q2_v1.root");
	TH1D* q2_mag = (TH1D*) file->Get("ana/q2_mag");


	double total = q2_mag->Integral();
	for(int i = 0; i < q2_mag->GetNbinsX(); i++){

		double num = q2_mag->Integral(i, q2_mag->GetNbinsX());

		double ratio = num/total;

		if( ratio < 1.0 && ratio > 0.99){

			cout << "bin number: " << i << endl;
		}
	}


}