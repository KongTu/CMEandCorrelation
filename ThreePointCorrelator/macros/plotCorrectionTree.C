#include "RiceStyle.h"
#include <vector>

using namespace std;

const int MAXTRACKS = 40000;
double etabins[] = {-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4};
const int Nbins = sizeof(etabins) / sizeof(etabins[0]) - 1;
double dEtaBins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8};
const int NdEtaBins = sizeof(dEtaBins) / sizeof(dEtaBins[0]) - 1;

double get3Real(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*R3;
  double t2 = R1*I2*I3;
  double t3 = R2*I1*I3;
  double t4 = I1*I2*R3;

  return t1-t2-t3-t4;

}

double get3Imag(double R1, double R2, double R3, double I1, double I2, double I3){

  double t1 = R1*R2*I3;
  double t2 = R1*R3*I2;
  double t3 = R2*R3*I1;
  double t4 = I1*I2*I3;

  return t1+t2+t3-t4;

}

double get2Real( double R1, double R2, double I1, double I2){

	double real = R1*R2 - I1*I2;
	return real;

}
double get2Imag( double R1, double R2, double I1, double I2){

	double imag = R1*I2 + R2*I1;
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

void plotCorrectionTree(){

	for(int deta = 0; deta < NdEtaBins+1; deta++){
		dEtaBins[deta] = dEtaBins[deta] - 0.00001;//fix bin boundary
	}


	TFile* file = new TFile("../rootfiles/test_tree_1.root");
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

	vector< vector<double>> Q_total_real;
	vector< vector<double>> Q_total_imag;

	double Q1[48][2][2];
	double Q3[2];
	int Q1_count[48][3]; int Q3_count = 0;
	for(int part = 0; part < 2; part++){
		Q3[part] = 0.;
			for(int eta = 0; eta < Nbins; eta++){
				for(int sign = 0; sign < 2; sign++ ){
					Q1[eta][sign][part] = 0.;
					Q1_count[eta][sign] = 0;  
				}
			}
	}

	TH1D* QvsdEta[48][3];
	TH1D* XY_real[48][3];TH1D* XY_imag[48][3];
	TH1D* XZ_real[48][3];TH1D* XZ_imag[48][3];
	TH1D* YZ_real[48][3];TH1D* YZ_imag[48][3];
	TH1D* X_real[48][3]; TH1D* X_imag[48][3];
	TH1D* Y_real[48][3]; TH1D* Y_imag[48][3];
	TH1D* Z_real[48][3]; TH1D* Z_imag[48][3];

	for(int deta = 0; deta < NdEtaBins; deta++){
		for( int sign = 0; sign < 3; sign++){
			QvsdEta[deta][sign] = new TH1D(Form("QvsdEta_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			
			XY_real[deta][sign] = new TH1D(Form("XY_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			XZ_real[deta][sign] = new TH1D(Form("XZ_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			YZ_real[deta][sign] = new TH1D(Form("YZ_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			
			XY_imag[deta][sign] = new TH1D(Form("XY_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			XZ_imag[deta][sign] = new TH1D(Form("XZ_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			YZ_imag[deta][sign] = new TH1D(Form("YZ_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			
			X_real[deta][sign] = new TH1D(Form("X_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			Y_real[deta][sign] = new TH1D(Form("Y_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			Z_real[deta][sign] = new TH1D(Form("Z_real_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			
			X_imag[deta][sign] = new TH1D(Form("X_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			Y_imag[deta][sign] = new TH1D(Form("Y_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
			Z_imag[deta][sign] = new TH1D(Form("Z_imag_%d_%d",deta,sign), "test", 20000,-1.0,1.0);
		}
	}

	for(int iEvent = 0; iEvent < 100; iEvent++){

		tree->GetEntry( iEvent );
		if( fabs(vtxZ) > 15.0 ) continue;
		int nTracks = 0;
		for(int itrk = 0; itrk < nTrk; itrk++ ){

			if( trkPtError[itrk] > 0.1 ) continue;
			if( trkDCAz[itrk] > 3.0 ) continue;
			if( trkDCAxy[itrk] > 3.0 ) continue;
			nTracks++;
			//if( trkPhi[itrk] < -1.5 ) continue;

			for(int eta = 0; eta < Nbins; eta++){
				if( trkEta[itrk] > etabins[eta] && trkEta[itrk] < etabins[eta+1] ){

					if( trkCharge[itrk] == 1){
						Q1[eta][0][0] += cos( trkPhi[itrk] );
						Q1[eta][0][1] += sin( trkPhi[itrk] );
						Q1_count[eta][0]++;
					}
					else if( trkCharge[itrk] == -1){
						Q1[eta][1][0] += cos( trkPhi[itrk] );
						Q1[eta][1][1] += sin( trkPhi[itrk] );
						Q1_count[eta][1]++;
					}	
				}
			}
		}
		Ntrk->Fill( nTracks );
		for(int itow = 0; itow < nTower; itow++){
			//if( towPhi[itow] < -1.5 ) continue;
			if( fabs(towEta[itow]) > 4.4 && fabs(towEta[itow]) < 5.0 ){
				Q3[0] += cos( -2*towPhi[itow] );
				Q3[1] += sin( -2*towPhi[itow] );
				Q3_count++;

			}
		}

		for(int ieta = 0; ieta < Nbins; ieta++ ){
			for(int jeta = 0; jeta < Nbins; jeta++){

				if(ieta == jeta ) continue;
				double deltaEta = fabs( etabins[ieta] - etabins[jeta] );

				for(int deta = 0; deta < NdEtaBins; deta++){

					if( deltaEta > dEtaBins[deta] && deltaEta < dEtaBins[deta+1] ){

						for(int sign = 0; sign < 2; sign++ ){

						//fill 3p correlators:
						if( Q1_count[ieta][sign] == 0 || Q1_count[jeta][sign] == 0 || Q3_count == 0 ) continue;
						
						double Q_real = get3Real(Q1[ieta][sign][0]/Q1_count[ieta][sign],Q1[jeta][sign][0]/Q1_count[jeta][sign],Q3[0]/Q3_count, Q1[ieta][sign][1]/Q1_count[ieta][sign], Q1[jeta][sign][1]/Q1_count[jeta][sign], Q3[1]/Q3_count);
						QvsdEta[deta][sign]->Fill( Q_real );	

						double XY_real_temp = get2Real(Q1[ieta][sign][0], Q1[jeta][sign][0], Q1[ieta][sign][1], Q1[jeta][sign][1]);
						double XY_imag_temp = get2Imag(Q1[ieta][sign][0], Q1[jeta][sign][0], Q1[ieta][sign][1], Q1[jeta][sign][1]);
						
						double XZ_real_temp = get2Real(Q1[ieta][sign][0], Q3[0], Q1[ieta][sign][1], Q3[1]);
						double XZ_imag_temp = get2Imag(Q1[ieta][sign][0], Q3[0], Q1[ieta][sign][1], Q3[1]);

						double YZ_real_temp = get2Real(Q1[jeta][sign][0], Q3[0], Q1[jeta][sign][1], Q3[1]);
						double YZ_imag_temp = get2Imag(Q1[jeta][sign][0], Q3[0], Q1[jeta][sign][1], Q3[1]);

						XY_real[deta][sign]->Fill( XY_real_temp/(Q1_count[ieta][sign]*Q1_count[jeta][sign]), Q1_count[ieta][sign]*Q1_count[jeta][sign] );
						XY_imag[deta][sign]->Fill( XY_imag_temp/(Q1_count[ieta][sign]*Q1_count[jeta][sign]), Q1_count[ieta][sign]*Q1_count[jeta][sign] );
						
						XZ_real[deta][sign]->Fill( XZ_real_temp/(Q1_count[ieta][sign]*Q3_count), Q1_count[ieta][sign]*Q3_count );
						XZ_imag[deta][sign]->Fill( XZ_imag_temp/(Q1_count[ieta][sign]*Q3_count), Q1_count[ieta][sign]*Q3_count );

						YZ_real[deta][sign]->Fill( YZ_real_temp/(Q1_count[jeta][sign]*Q3_count), Q1_count[jeta][sign]*Q3_count );
						YZ_imag[deta][sign]->Fill( YZ_imag_temp/(Q1_count[jeta][sign]*Q3_count), Q1_count[jeta][sign]*Q3_count );

						double X_real_temp = Q1[ieta][sign][0]; double X_imag_temp = Q1[ieta][sign][1]; 
						double Y_real_temp = Q1[jeta][sign][0]; double Y_imag_temp = Q1[jeta][sign][1]; 
						double Z_real_temp = Q3[0];       double Z_imag_temp = Q3[1]; 

						X_real[deta][sign]->Fill( X_real_temp/Q1_count[ieta][sign], Q1_count[ieta][sign]); 	  
						Y_real[deta][sign]->Fill( Y_real_temp/Q1_count[jeta][sign], Q1_count[jeta][sign]); 	  
						Z_real[deta][sign]->Fill( Z_real_temp/Q3_count, Q3_count); 	
					
						X_imag[deta][sign]->Fill( X_imag_temp/Q1_count[ieta][sign], Q1_count[ieta][sign]); 	  
						Y_imag[deta][sign]->Fill( Y_imag_temp/Q1_count[jeta][sign], Q1_count[jeta][sign]); 	  
						Z_imag[deta][sign]->Fill( Z_imag_temp/Q3_count, Q3_count);
						}

					}
				}
			}
		}
		cout << "ievent: " << iEvent << endl;
	}

	TH1D* hist1[2];
	TH1D* hist2[2];
	for(int sign = 0; sign < 2; sign++){
		hist1[sign] = new TH1D(Form("hist1_%d",sign),"test", NdEtaBins, dEtaBins);
		hist2[sign] = new TH1D(Form("hist2_%d",sign),"test", NdEtaBins, dEtaBins);
	}

	for(int deta = 0; deta < NdEtaBins; deta++){

		if(deta == 0) continue;

		for(int sign = 0; sign < 2; sign++){

			double Q_total_real_dEta = QvsdEta[deta][sign]->GetMean();
			double Q_total_real_dEta_error = QvsdEta[deta][sign]->GetMeanError();
			
			double XY_real_temp = 0.; double XZ_real_temp = 0.; double YZ_real_temp = 0.; double X_real_temp = 0.; double Y_real_temp = 0.; double Z_real_temp = 0.;
			double XY_imag_temp = 0.; double XZ_imag_temp = 0.; double YZ_imag_temp = 0.; double X_imag_temp = 0.; double Y_imag_temp = 0.; double Z_imag_temp = 0.;

			XY_real_temp = XY_real[deta][sign]->GetMean();
			XZ_real_temp = XZ_real[deta][sign]->GetMean();
			YZ_real_temp = YZ_real[deta][sign]->GetMean();

			XY_imag_temp = XY_imag[deta][sign]->GetMean();
			XZ_imag_temp = XZ_imag[deta][sign]->GetMean();
			YZ_imag_temp = YZ_imag[deta][sign]->GetMean();

			X_real_temp = X_real[deta][sign]->GetMean();
			Y_real_temp = Y_real[deta][sign]->GetMean();
			Z_real_temp = Z_real[deta][sign]->GetMean();

			X_imag_temp = X_imag[deta][sign]->GetMean();
			Y_imag_temp = Y_imag[deta][sign]->GetMean();
			Z_imag_temp = Z_imag[deta][sign]->GetMean();

			double term1_real = get2Real(XY_real_temp, Z_real_temp, XY_imag_temp, Z_imag_temp);
			double term2_real = get2Real(XZ_real_temp, Y_real_temp, XZ_imag_temp, Y_imag_temp);
			double term3_real = get2Real(YZ_real_temp, X_real_temp, YZ_imag_temp, X_imag_temp);
			double term4_real = get3Real(X_real_temp, Y_real_temp, Z_real_temp, X_imag_temp, Y_imag_temp, Z_imag_temp);
			double all_real = term1_real+term2_real+term3_real-2*term4_real;

			if( Q_total_real_dEta == 0.00000 ) continue;

			cout << "before correction is " << Q_total_real_dEta << endl;
			cout << "correction factor is: " << all_real << endl;
			cout << "after correction is " << Q_total_real_dEta - all_real << endl;
			cout << "-------- " << endl;
			hist1[sign]->SetBinContent(deta+1, Q_total_real_dEta );
			hist1[sign]->SetBinError(deta+1,  Q_total_real_dEta_error);
			hist2[sign]->SetBinContent(deta+1, Q_total_real_dEta - all_real);
			hist2[sign]->SetBinError(deta+1, Q_total_real_dEta_error);
		}
	}

	TCanvas* c1 = new TCanvas();
	hist1[0]->SetMarkerColor(kBlack);
	hist1[0]->SetMarkerStyle(20);
	hist1[0]->GetYaxis()->SetRangeUser(-0.5, 0.5);
	hist1[0]->GetYaxis()->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3})");
	hist1[0]->GetXaxis()->SetTitle("#Delta#eta");
	hist1[0]->Draw("P");

	hist2[0]->SetMarkerColor(kRed);
	hist2[0]->SetMarkerStyle(20);
	hist2[0]->Draw("Psame");

	TCanvas* c2 = new TCanvas();
	hist1[1]->SetMarkerColor(kBlack);
	hist1[1]->SetMarkerStyle(20);
	hist1[1]->GetYaxis()->SetRangeUser(-0.5, 0.5);
	hist1[1]->GetYaxis()->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3})");
	hist1[1]->GetXaxis()->SetTitle("#Delta#eta");
	hist1[1]->Draw("P");

	hist2[1]->SetMarkerColor(kRed);
	hist2[1]->SetMarkerStyle(20);
	hist2[1]->Draw("Psame");


}