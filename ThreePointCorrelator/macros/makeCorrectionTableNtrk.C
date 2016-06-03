#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

//CASE 1 and 2, cos1 and sin1 denotes sum of cos(phi) and sin(phi), and cos2 and sin2 denotes sum of cos(2phi) and sin(2phi)
//CASE 3, cos1 and sin1 denotes sum of cos(phi+) and sin(phi+), and cos2 and sin2 denote sum of cos(phi-) and sin(phi-)

double correctionFactor(double cos1, double sin1, double cos2, double sin2, double HFcos, double HFsin, int CASE){

	double realPart1 = 0.;
	double imagPart1 = 0.;

	double realPart2 = 0.;
	double imagPart2 = 0.;

	if( CASE == 1 || CASE == 2){

		realPart1 = cos1*cos1 - sin1*sin1 - cos2;
		imagPart1 = 2*cos1*sin1 - sin2;

		realPart2 = realPart1 * HFcos + imagPart1 * HFsin;

		return realPart2;

	}
	else if( CASE == 3 ){

		realPart1 = cos1*cos2 - sin1*sin2;
		imagPart1 = sin1*cos2 + cos1*sin2;

		realPart2 = realPart1 * HFcos + imagPart1 * HFsin;

		return realPart2;
	}
}

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

    //case1 (++) with HF+
    double factor1plus = correctionFactor(cosPlusSum, sinPlusSum, cos2PlusSum, sin2PlusSum, HFcosSumPlus, HFsinSumPlus, 1);
    factor1plus = factor1plus/(plusSumEntries*(plusSumEntries-1));
    //case1 (++) with HF-
    double factor1minus = correctionFactor(cosPlusSum, sinPlusSum, cos2PlusSum, sin2PlusSum, HFcosSumMinus, HFsinSumMinus, 1);
    factor1minus = factor1minus/(plusSumEntries*(plusSumEntries-1));

    //case2 (--) with HF+
    double factor2plus = correctionFactor(cosMinusSum, sinMinusSum, cos2MinusSum, sin2MinusSum, HFcosSumPlus, HFsinSumPlus, 2);
    factor2plus = factor2plus/(minusSumEntries*(minusSumEntries-1));
    //case2 (--) with HF-
    double factor2minus = correctionFactor(cosMinusSum, sinMinusSum, cos2MinusSum, sin2MinusSum, HFcosSumMinus, HFsinSumMinus, 2);
    factor2minus = factor2minus/(minusSumEntries*(minusSumEntries-1));

    //case3 (+/-) with HF+
    double factor3plus = correctionFactor(cosPlusSum, sinPlusSum, cosMinusSum, sinMinusSum, HFcosSumPlus, HFsinSumPlus, 3);
    factor3plus = factor3plus/(plusSumEntries*minusSumEntries);
    //case3 (+/-) with HF-
    double factor3minus = correctionFactor(cosPlusSum, sinPlusSum, cosMinusSum, sinMinusSum, HFcosSumMinus, HFsinSumMinus, 3);
    factor3minus = factor3minus/(plusSumEntries*minusSumEntries);

    ofstream myfile("ntrk_185_220.txt");
	if( myfile.is_open() ){

			myfile << factor1plus << " "
			       << factor1minus << " " 
			       << factor2plus << " " 
			       << factor2minus << " " 
			       << factor3plus << " " 
			       << factor3minus << endl; 
			        
	}	

	myfile.close();

}