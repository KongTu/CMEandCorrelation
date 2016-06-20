#include "RiceStyle.h"
#include <vector>

using namespace std;

const int MAXTRACKS = 40000;
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;


double get3Real(double R1, double R2, double R3, double I1, double I2, double I3) 
{

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;
}
double get3RealOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*R3;
      double t2 = (2*R1*I1-I2)*I3;
      double N = N1*(N1-1)*N3;

      return (t1-t2)/N;
}
double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}
double get3ImagOverlap(double R1, double R2, double R3, double I1, double I2, double I3, double N1, double N3){

      double t1 = (R1*R1 - I1*I1 - R2)*I3;
      double t2 = (2*R1*I1-I2)*R3;
      double N = N1*(N1-1)*N3;

      return (t1+t2)/N;
}

double get2Real( double R1, double R2, double I1, double I2){

      double real = R1*R2 - I1*I2;
      return real;

}
double get2RealOverlap( double R1, double R2, double I1, double I2){

      double real = R1*R1 - I1*I1 - R2;
      return real;

}
double get2Imag( double R1, double R2, double I1, double I2){

      double imag = R1*I2 + R2*I1;
      return imag;
}
double get2ImagOverlap( double R1, double R2, double I1, double I2){

      double imag = (2*R1*I1-I2);
      return imag;
} 

double getAverage( vector<double> v ){

	int size = v.size();
	double sum = 0.;
	for(int i = 0; i < size; i ++){

		sum = sum + v[i];
	}

	return sum/size;
}

double getAverageHist( TH1D* hist ){

	double integral2 = 0.;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		integral2 += temp*center;
	
	}
	double ave = integral2/( hist->GetEntries() );
	return ave;

}

void plotCorrectionTreeNestedLoopTest(){

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/test_tree.root");
	TTree* tree = (TTree*)file->Get("ana_tree/trackTree");

	int nTrk;
	int nTower;

	double vtxZ;
	double trkEta[MAXTRACKS];
	double trkPhi[MAXTRACKS];
	int    trkCharge[MAXTRACKS];
	double trkPt[MAXTRACKS];
	double trkPtError[MAXTRACKS];
	double trkDCAxy[MAXTRACKS];
	double trkDCAz[MAXTRACKS];

	double towEta[MAXTRACKS];
	double towPhi[MAXTRACKS];
	double towEt[MAXTRACKS];
	double towEnergy[MAXTRACKS];

	tree->SetBranchAddress("nTrk", &nTrk);
	tree->SetBranchAddress("nTower", &nTower);
	tree->SetBranchAddress("vtxZ", &vtxZ);

	tree->SetBranchAddress("trkEta", &trkEta);
	tree->SetBranchAddress("trkPhi", &trkPhi);
	tree->SetBranchAddress("trkCharge", &trkCharge);
	tree->SetBranchAddress("trkPt", &trkPt);
	tree->SetBranchAddress("trkPtError", &trkPtError);
	tree->SetBranchAddress("trkDCAxy", &trkDCAxy);
	tree->SetBranchAddress("trkDCAz", &trkDCAz);

	tree->SetBranchAddress("towEta", &towEta);
	tree->SetBranchAddress("towPhi", &towPhi);
	tree->SetBranchAddress("towEt", &towEt);
	tree->SetBranchAddress("towEnergy", &towEnergy);

	int numberOfEvents = tree->GetEntries();

	TH1D* Ntrk = new TH1D("Ntrk", "", 400,0,400);

	double Q1[48][2][2];
	double Q2[48][2][2];
	double Q3[2][2];
	double Q1_count[48][2]; double ETT[2] = 0.;

	TH1D* QvsdEta[48][3][2];

	for(int deta = 0; deta < NdEtaBins; deta++){
		for( int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				QvsdEta[deta][sign][HF] = new TH1D(Form("QvsdEta_%d_%d_%d",deta,sign,HF), "test", 20000,-1.0,1.0);
		
			}
		}
	}
	for(int iEvent = 0; iEvent < 1000; iEvent++){

		tree->GetEntry( iEvent );

		if( fabs(vtxZ) > 15.0 ) continue;
		
		double real_term = 0.;
		double real_count = 0.;
		for(int itrk = 0; itrk < nTrk; itrk++ ){

			if( trkPtError[itrk] > 0.1 ) continue;
			if( trkDCAz[itrk] > 3.0 ) continue;
			if( trkDCAxy[itrk] > 3.0 ) continue;
			if( trkPt[itrk] < 0.4 || fabs(trkEta[itrk]) > 2.4 ) continue;
			if( trkCharge[itrk] != 1 ) continue;

			for(int jtrk = 0; jtrk < nTrk; jtrk++){

				if( trkPtError[jtrk] > 0.1 ) continue;
				if( trkDCAz[jtrk] > 3.0 ) continue;
				if( trkDCAxy[jtrk] > 3.0 ) continue;
				if( trkPt[jtrk] < 0.4 || fabs(trkEta[jtrk]) > 2.4 ) continue;
				if( trkCharge[jtrk] != 1 ) continue;

				if( itrk == jtrk ) continue;
				if( fabs(trkEta[itrk] - trkEta[jtrk]) >= 0.1 ) continue;
			
				for(int itow = 0; itow < nTower; itow++){
					//double w = towEt[itow];

					if( towEta[itow] > 4.4 && towEta[itow] < 5.0 ){
						
						real_term += cos(trkPhi[itrk]+trkPhi[jtrk] - 2*towPhi[itow]);
						real_count ++;
					}
				}

			}
		}

        QvsdEta[0][0][0]->Fill( real_term/real_count , real_count);

		cout << "ievent: " << iEvent << endl;
		//cout << "number of pairs: " << real_count << endl;

	}

	cout << "the integrated should be: " << endl;
	cout << "value: " << QvsdEta[0][0][0]->GetMean() << endl;





}