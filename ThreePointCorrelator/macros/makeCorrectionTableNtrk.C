#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

void makeCorrectionTableNtrk(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_Qvsntrk_pPb_HM_185_220_v3.root");

	TH1D* TRKcosPlusSum = (TH1D*) file->Get("ana/TRKcosPlusSum");
	TH1D* TRKsinPlusSum = (TH1D*) file->Get("ana/TRKsinPlusSum");
	TH1D* TRK2cosPlusSum = (TH1D*) file->Get("ana/TRK2cosPlusSum");
	TH1D* TRK2sinPlusSum = (TH1D*) file->Get("ana/TRK2sinPlusSum");
		
	TH1D* TRKcosMinusSum = (TH1D*) file->Get("ana/TRKcosMinusSum");
	TH1D* TRKsinMinusSum = (TH1D*) file->Get("ana/TRKsinMinusSum");
	TH1D* TRK2cosMinusSum = (TH1D*) file->Get("ana/TRK2cosMinusSum");
	TH1D* TRK2sinMinusSum = (TH1D*) file->Get("ana/TRK2sinMinusSum");

	vector<TH1D*> HFcosSum = loadingHistogram(file, "ana/HFcosSum_", 2);
	vector<TH1D*> HFsinSum = loadingHistogram(file, "ana/HFsinSum_", 2);
	vector<TH1D*> weightSum = loadingHistogram(file, "ana/weightSum_", 2);

	ofstream myfile("ntrk_185_220.txt");

    double HFcosSumPlus = getWeightedSum( HFcosSum[0], weightSum[0] );
    double HFsinSumPlus = getWeightedSum( HFsinSum[0], weightSum[0] );
    double HFcosSumMinus = getWeightedSum( HFcosSum[1], weightSum[1] );
    double HFsinSumMinus = getWeightedSum( HFsinSum[1], weightSum[1] );

    double cosPlusSum = getNormalizedSumOnly( TRKcosPlusSum );
    double sinPlusSum = getNormalizedSumOnly( TRKsinPlusSum );
    double cos2PlusSum = getNormalizedSumOnly( TRK2cosPlusSum );
    double sin2PlusSum = getNormalizedSumOnly( TRK2sinPlusSum );

    double cosMinusSum = getNormalizedSumOnly( TRKcosMinusSum );
    double sinMinusSum = getNormalizedSumOnly( TRKsinMinusSum );
    double cos2MinusSum = getNormalizedSumOnly( TRK2cosMinusSum );
    double sin2MinusSum = getNormalizedSumOnly( TRK2sinMinusSum );

    double plusSumEntries = TRKcosPlusSum->GetEntries();
    double minusSumEntries = TRKcosMinusSum->GetEntries();



	if( myfile.is_open() ){

			myfile << cosPlusSum << " "
			       << sinPlusSum << " " 
			       << cos2PlusSum << " " 
			       << sin2PlusSum << " " 
			       << cosMinusSum << " " 
			       << sinMinusSum << " " 
			       << cos2MinusSum << " " 
			       << sin2MinusSum << " "
			       << plusSumEntries << " "			   
			       << minusSumEntries << " "
			       << HFcosSumPlus << " "
			       << HFsinSumPlus << " "
			       << HFcosSumMinus << " "
			       << HFsinSumMinus << endl;
	}	

	myfile.close();

}