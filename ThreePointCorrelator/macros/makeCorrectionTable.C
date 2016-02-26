#include "RiceStyle.h"
#include "function.h"

using namespace std;

void makeCorrectionTable(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_HM_v12.root");

	vector<TH1D*> TRKcosPlusSum = loadingHistogram(file, "ana/TRKcosPlusSum_", Nbins);
	vector<TH1D*> TRKsinPlusSum = loadingHistogram(file, "ana/TRKsinPlusSum_", Nbins);
	vector<TH1D*> TRKcosMinusSum = loadingHistogram(file, "ana/TRKcosMinusSum_", Nbins);
	vector<TH1D*> TRKsinMinusSum = loadingHistogram(file, "ana/TRKsinMinusSum_", Nbins);

	ofstream myfile("v12_table.txt");

	for(int eta = 0; eta < Nbins; eta++){

		double cosPlusSum = getNormalizedSum( TRKcosPlusSum[eta] );
		double sinPlusSum = getNormalizedSum( TRKsinPlusSum[eta] );
		double cosMinusSum = getNormalizedSum( TRKcosMinusSum[eta] );
		double sinMinusSum = getNormalizedSum( TRKsinMinusSum[eta] );

		if( myfile.is_open() ){

			myfile << cosPlusSum << " " << sinPlusSum << " " << cosMinusSum << " " << sinMinusSum << endl;
		}	

	}

	myfile.close();

}