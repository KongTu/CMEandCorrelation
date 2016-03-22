#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

void makeCorrectionTable(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_PbPb_30_40_v10.root");

	vector<TH1D*> TRKcosPlusSum = loadingHistogram(file, "ana/TRKcosPlusSum_", Nbins);
	vector<TH1D*> TRKsinPlusSum = loadingHistogram(file, "ana/TRKsinPlusSum_", Nbins);
	vector<TH1D*> TRKcosMinusSum = loadingHistogram(file, "ana/TRKcosMinusSum_", Nbins);
	vector<TH1D*> TRKsinMinusSum = loadingHistogram(file, "ana/TRKsinMinusSum_", Nbins);

	vector<TH1D*> HFcosSum = loadingHistogram(file, "ana/HFcosSum_", 2);
	vector<TH1D*> HFsinSum = loadingHistogram(file, "ana/HFsinSum_", 2);
	vector<TH1D*> weightSum = loadingHistogram(file, "ana/weightSum_", 2);

	ofstream myfile("PbPb_v10.txt");

    double HFcosSumPlus = getWeightedSum( HFcosSum[0], weightSum[0] );
    double HFsinSumPlus = getWeightedSum( HFsinSum[0], weightSum[0] );
    double HFcosSumMinus = getWeightedSum( HFcosSum[1], weightSum[1] );
    double HFsinSumMinus = getWeightedSum( HFsinSum[1], weightSum[1] );

	for(int eta = 0; eta < Nbins; eta++){

		double cosPlusSum = getNormalizedSum( TRKcosPlusSum[eta] );
		double sinPlusSum = getNormalizedSum( TRKsinPlusSum[eta] );
		double cosMinusSum = getNormalizedSum( TRKcosMinusSum[eta] );
		double sinMinusSum = getNormalizedSum( TRKsinMinusSum[eta] );

		if( myfile.is_open() ){

			myfile << cosPlusSum << " " << sinPlusSum << " " << cosMinusSum << " " << sinMinusSum << 
			" " << HFcosSumPlus << " " << HFsinSumPlus << " " << HFcosSumMinus << " " << HFsinSumMinus << endl;
		}	

	}

	myfile.close();

}