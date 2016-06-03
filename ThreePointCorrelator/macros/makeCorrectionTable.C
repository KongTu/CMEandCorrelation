#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

void makeCorrectionTable(){

	int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;

	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_v15_2nd.root");

	//QaQb2Dhistograms:
	TH2D* QaQbvsdEta[5][5];
	TH2D* NaNbvsdEta[5];

	//QaQc2Dhistograms:
	TH2D* QaQcvsdEta[5][5][5];
	TH2D* QbvsdEta[5][5][5];
	TH2D* NaNcvsdEta[5][5];

	//QaQbsingle2Dhistogram
	TH2D* QaSinglevsdEta[5][5];
	TH2D* QbSinglevsdEta[5][5];

	for(int type = 0; type < 3; type++){
		for(int sign = 0; sign < 3; sign++){

	 		 NaNcvsdEta[type][sign] = (TH2D*) file->Get( Form("ana/NaNcvsdEta_%d_%d", type, sign) );
		}
	}

	for(int sign = 0; sign < 3; sign++){

		NaNbvsdEta[sign] = (TH2D*)file->Get( Form("ana/NaNbvsdEta_%d", sign)  ); 

		for(int j = 0; j < 2; j++){

			QaQbvsdEta[sign][j] = (TH2D*)file->Get( Form("ana/QaQbvsdEta_%d_%d", sign, j) );
			QaSinglevsdEta[sign][j] = (TH2D*)file->Get( Form("ana/QaSinglevsdEta_%d_%d", sign, j) );
			QbSinglevsdEta[sign][j] = (TH2D*)file->Get( Form("ana/QbSinglevsdEta_%d_%d", sign, j) );

			for(int type = 0; type < 3; type++){

				QaQcvsdEta[type][sign][j] = (TH2D*)file->Get( Form("ana/QaQcvsdEta_%d_%d_%d", type, sign, j) );
				QbvsdEta[type][sign][j] = (TH2D*)file->Get( Form("ana/QbvsdEta_%d_%d_%d", type, sign, j) );
			}

		}
	}

	//HF Qvector, with Pb and p going side separated
	vector<TH1D*> HFcosSum = loadingHistogram(file, "ana/HFcosSum_", 2);
	vector<TH1D*> HFsinSum = loadingHistogram(file, "ana/HFsinSum_", 2);
	vector<TH1D*> weightSum = loadingHistogram(file, "ana/weightSum_", 2);

    double HFcosSumPlus = getWeightedSum( HFcosSum[0], weightSum[0] );
    double HFsinSumPlus = getWeightedSum( HFsinSum[0], weightSum[0] );
    double HFcosSumMinus = getWeightedSum( HFcosSum[1], weightSum[1] );
    double HFsinSumMinus = getWeightedSum( HFsinSum[1], weightSum[1] );
    //end HFQvector

    //Using HF plus:
    
    TH1D* temp1Dhist;
    TH1D* temp2Dhist;
    TH1D* temp3Dhist;
    TH1D* temp4Dhist;
    TH1D* temp5Dhist;

    double term1_HFplus, term2, term3, term4_HFplus;
    double term1_HFminus, term2, term3, term4_HFminus;

    vector< vector<double> > VdEtaTerm1, VdEtaTerm2, VdEtaTerm4;
    vector<double> signVector;

	for(int deta = 0; deta < NdEtaBins; deta++){
    	 //term1
    	 for(int sign = 0; sign < 3; sign++){
			
			temp1Dhist= (TH1D*)QaQbvsdEta[sign][0]->ProjectionY(Form("QaQbvsdEta_%d_0_dEta_%d",sign,deta), deta+1, deta+1 );
			double real = getSumOnly(temp1Dhist);

		 	temp2Dhist= (TH1D*)QaQbvsdEta[sign][1]->ProjectionY(Form("QaQbvsdEta_%d_1_dEta_%d",sign,deta), deta+1, deta+1 );
			double imag = getSumOnly(temp2Dhist);

			temp3Dhist = (TH1D*)NaNbvsdEta[sign]->ProjectionY( Form("NaNbvsdEta_%d_dEta_%d",sign,deta), deta+1, deta+1 );
			double sum = getSumOnly( temp3Dhist );

			real /= sum; 
			imag /= sum;

			term1_HFplus = real*HFcosSumPlus + imag*HFsinSumPlus;
			term1_HFminus = real*HFcosSumMinus + imag*HFsinSumMinus;

			signVector.push_back( term1_HFplus ); signVector.push_back( term1_HFminus );

    	}

    	VdEtaTerm1.push_back( signVector );
    	signVector.clear();

    	//term2 and term3
		for(int sign = 0; sign < 3; sign++){
			for(int type = 0; type < 2; type++){

    			temp1Dhist = (TH1D*)QaQcvsdEta[type][sign][0]->ProjectionY(Form("QaQcvsdEta_%d_%d_0_dEta_%d",type,sign,deta), deta+1, deta+1);
    			double real = getSumOnly(temp1Dhist);

    			temp2Dhist = (TH1D*)QaQcvsdEta[type][sign][1]->ProjectionY(Form("QaQcvsdEta_%d_%d_1_dEta_%d",type,sign,deta), deta+1, deta+1);
				double imag = getSumOnly(temp2Dhist);

				temp3Dhist = (TH1D*)NaNcvsdEta[type][sign]->ProjectionY( Form("NaNcvsdEta_%d_%d_dEta_%d",type,sign,deta), deta+1, deta+1 );
				double sum = getSumOnly( temp3Dhist );

				temp4Dhist = (TH1D*)QbvsdEta[type][sign][0]->ProjectionY( Form("QbvsdEta_%d_%d_0_dEta_%d",type,sign,deta), deta+1, deta+1);
    			double real_b = getNormalizedSum(temp4Dhist);

				temp5Dhist = (TH1D*)QbvsdEta[type][sign][1]->ProjectionY( Form("QbvsdEta_%d_%d_1_dEta_%d",type,sign,deta), deta+1, deta+1);
				double imag_b = getNormalizedSum(temp5Dhist);

				real /= sum;
				imag /= sum;

				term2 = real*real_b - imag*imag_b;

				signVector.push_back( 2*term2 );

    		}
    	}

		VdEtaTerm2.push_back( signVector );
		signVector.clear();

		//term3
		
		for(int sign = 0; sign < 3; sign++){

			temp1Dhist= (TH1D*)QaSinglevsdEta[sign][0]->ProjectionY(Form("QaSinglevsdEta_%d_0_dEta_%d",sign,deta), deta+1, deta+1 );
			double real_a = getNormalizedSum(temp1Dhist);

		 	temp2Dhist= (TH1D*)QaSinglevsdEta[sign][1]->ProjectionY(Form("QaSinglevsdEta_%d_1_dEta_%d",sign,deta), deta+1, deta+1 );
			double imag_a = getNormalizedSum(temp2Dhist);

			temp3Dhist= (TH1D*)QbSinglevsdEta[sign][0]->ProjectionY(Form("QbSinglevsdEta_%d_0_dEta_%d",sign,deta), deta+1, deta+1 );
			double real_b = getNormalizedSum(temp3Dhist);

		 	temp4Dhist= (TH1D*)QbSinglevsdEta[sign][1]->ProjectionY(Form("QbSinglevsdEta_%d_1_dEta_%d",sign,deta), deta+1, deta+1 );
			double imag_b = getNormalizedSum(temp4Dhist);

			double real = real_a*real_b - imag_a*imag_b;
			double imag = real_a*imag_b + real_b*imag_a;

			term4_HFplus = real*HFcosSumPlus + imag*HFsinSumPlus;
			term4_HFminus = real*HFcosSumMinus + imag*HFsinSumMinus;

			signVector.push_back(2*term4_HFplus); signVector.push_back(2*term4_HFminus);
		}

		VdEtaTerm4.push_back( signVector );
		signVector.clear();
    }

	ofstream myfile("EPOS_v15_2nd.txt");
	if( myfile.is_open() ){
		for(int deta = 0; deta < NdEtaBins; deta++){

			double v[6];
			for(int i = 0; i < 6; i++){

				v[i] = VdEtaTerm1[deta][i] + VdEtaTerm2[deta][i] - VdEtaTerm4[deta][i];
			}

			if(deta == 0 ) myfile << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
			else{
			myfile << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " " << v[4] << " " << v[5] << endl;}
			//       HF+, ++        HF-,++         HF+,--         HF-,--         HF+,+-         HF-,+-
		}
	}

	myfile.close();	

	

}