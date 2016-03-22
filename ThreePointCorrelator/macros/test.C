#include "RiceStyle.h"
#include "function.C"
#include "inputHistogram.h"

using namespace std;

void test(){


	TFile* file = new TFile("~/Desktop/test.root");

	TGraph* gr = (TGraph*) file->Get("ana/Graph;1");

	int size = gr->GetN();

	vector<double> xV,yV;

	for(int i = 0; i < size; i++ ){
		double x, y;
		gr->GetPoint(i, x, y);
		
		if( x == 0.1 ){
			yV.push_back( y );
		}		
	}

	double sum = 0.;
	for(int j = 0; j < yV.size(); j++ ){

		sum += yV[j];
	}

	cout << "average: " << sum/yV.size() << endl;


}