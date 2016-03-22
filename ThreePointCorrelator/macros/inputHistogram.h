#include "TH1.h"
#include <sstream>
#include <vector>
#include "TString.h"
#include "TFile.h"

using namespace std;

vector<TH1D*> loadingHistogram( TFile * file, TString hName, int Nbins ){

	vector<TH1D *> spectra;

  	stringstream HistName;

  	for (int mult = 0; mult < Nbins; mult++){

	    HistName.str("");

	    HistName << hName;
	    HistName << mult;

	    TH1D* temp = (TH1D*)file->Get( HistName.str().c_str() );
	    temp->SetMarkerSize(1.3);
		temp->SetMarkerStyle(20);
		temp->SetStats(kFALSE);
		temp->SetTitle("");
		temp->SetXTitle("");
		temp->SetYTitle("");
	    
	    spectra.push_back( temp );
  	}

  	return spectra;

}

class FitDataPoint
	{
	public:
		double x, y, ey, bw;

		FitDataPoint(){
			x = 0;
			y = 0;
			ey = 0;
			bw = 0;
		}
		FitDataPoint( double _x, double _y, double _ey, double _bw = 1.0 ){
			x = _x;
			y = _y;
			ey = _ey;
			bw = _bw;
		}
		~FitDataPoint(){}
		
	};