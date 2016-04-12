#include "RiceStyle.h"
#include <vector>

using namespace std;

const int MAXTRACKS = 40000;

int Ntrkbins[] = {10,20,30,40,50,60,80,100};
double Ntrkbins_fill[] = {10,20,30,40,50,60,80,100};
const int nNtrkBins = sizeof(Ntrkbins)/ sizeof(Ntrkbins[0]) - 1; 

double get3Real(double R1, double R2, double R3, double I1, double I2, double I3){

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

void plotCorrectionTreeNtrk(){

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

	double Q1[2][2];
	double Q2[2][2];
	double Q3[2][2];
	double Q1_count[2]; double Q3_count[2] = 0.;

	TH1D* QvsNtrk[20][3][2];
	TH1D* XY_real[20][3][2];TH1D* XY_imag[20][3][2];
	TH1D* XZ_real[20][3][2];TH1D* XZ_imag[20][3][2];
	TH1D* YZ_real[20][3][2];TH1D* YZ_imag[20][3][2];
	TH1D* X_real[20][3][2]; TH1D* X_imag[20][3][2];
	TH1D* Y_real[20][3][2]; TH1D* Y_imag[20][3][2];
	TH1D* Z_real[20][3][2]; TH1D* Z_imag[20][3][2];

	for(int trk = 0; trk < nNtrkBins; trk++){
		for( int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){

				QvsNtrk[trk][sign][HF] = new TH1D(Form("QvsNtrk_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				
				XY_real[trk][sign][HF] = new TH1D(Form("XY_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				XZ_real[trk][sign][HF] = new TH1D(Form("XZ_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				YZ_real[trk][sign][HF] = new TH1D(Form("YZ_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				
				XY_imag[trk][sign][HF] = new TH1D(Form("XY_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				XZ_imag[trk][sign][HF] = new TH1D(Form("XZ_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				YZ_imag[trk][sign][HF] = new TH1D(Form("YZ_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				
				X_real[trk][sign][HF] = new TH1D(Form("X_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				Y_real[trk][sign][HF] = new TH1D(Form("Y_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				Z_real[trk][sign][HF] = new TH1D(Form("Z_real_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				
				X_imag[trk][sign][HF] = new TH1D(Form("X_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				Y_imag[trk][sign][HF] = new TH1D(Form("Y_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
				Z_imag[trk][sign][HF] = new TH1D(Form("Z_imag_%d_%d_%d",trk,sign,HF), "test", 20000,-1.0,1.0);
			}
		}
	}
	
	for(int iEvent = 0; iEvent < numberOfEvents; iEvent++){

		tree->GetEntry( iEvent );

		for(int HF = 0; HF < 2; HF++){
			for(int part = 0; part < 2; part++){
				Q3[part][HF] = 0.;
			}
			Q3_count[HF] = 0.;
		}
		for(int sign = 0; sign < 2; sign++ ){
			Q1_count[sign] = 0.;
			for(int part = 0; part < 2; part++){
				Q1[sign][part] = 0.;
				Q2[sign][part] = 0.;	  
			}
		}
		if( fabs(vtxZ) > 15.0 ) continue;
		int nTracks = 0;
		for(int itrk = 0; itrk < nTrk; itrk++ ){

			if( trkPtError[itrk] > 0.1 ) continue;
			if( trkDCAz[itrk] > 3.0 ) continue;
			if( trkDCAxy[itrk] > 3.0 ) continue;
			if( trkPt[itrk] < 0.4 || fabs(trkEta[itrk]) > 2.4 ) continue;
			nTracks++;
			if( trkPhi[itrk] < -1.5 ) continue;

			if( trkCharge[itrk] == 1 ){
				Q1[0][0] += cos( trkPhi[itrk] );
				Q1[0][1] += sin( trkPhi[itrk] );
				Q1_count[0]++;

				Q2[0][0] += cos( 2*trkPhi[itrk] );
				Q2[0][1] += sin( 2*trkPhi[itrk] );

			}
			else if( trkCharge[itrk] == -1 ){
				Q1[1][0] += cos( trkPhi[itrk] );
				Q1[1][1] += sin( trkPhi[itrk] );
				Q1_count[1]++;

				Q2[1][0] += cos( 2*trkPhi[itrk] );
				Q2[1][1] += sin( 2*trkPhi[itrk] );
			}

		}

		Ntrk->Fill( nTracks );
		if(nTracks < 10 ) continue;

		for(int itow = 0; itow < nTower; itow++){
			if( towPhi[itow] < -1.5 ) continue;
			double w = towEt[itow];

			if( towEta[itow] > 4.4 && towEta[itow] < 5.0 ){
				Q3[0][0] += w*cos( -2*towPhi[itow] );
				Q3[0][1] += w*sin( -2*towPhi[itow] );
				Q3_count[0] += w;

			}
			else if (towEta[itow] < -4.4 && towEta[itow] > -5.0){
				Q3[1][0] += w*cos( -2*towPhi[itow] );
				Q3[1][1] += w*sin( -2*towPhi[itow] );
				Q3_count[1] += w;

			}
		}

		for(int trk = 0; trk < nNtrkBins; trk++){
			if(  nTracks >= Ntrkbins[trk] &&  nTracks < Ntrkbins[trk+1] ){
				
				for(int HF = 0; HF < 2; HF++){
					for(int sign = 0; sign < 2; sign++){
						if( Q1_count[sign] == 0 || Q3_count[HF] == 0 ) continue;
						double Q_real = 0.;
						Q_real = get3RealOverlap(Q1[sign][0], Q2[sign][0], Q3[HF][0], Q1[sign][1], Q2[sign][1], Q3[HF][1], Q1_count[sign],Q3_count[HF]);
						QvsNtrk[trk][sign][HF]->Fill( Q_real );

						double XY_real_temp = get2RealOverlap(Q1[sign][0],Q2[sign][0],Q1[sign][1],Q2[sign][1]);
						double XZ_real_temp = get2Real(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);
						double YZ_real_temp = get2Real(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);

						double XY_imag_temp = get2ImagOverlap(Q1[sign][0],Q2[sign][0],Q1[sign][1],Q2[sign][1]);
						double XZ_imag_temp = get2Imag(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);
						double YZ_imag_temp = get2Imag(Q1[sign][0],Q3[HF][0],Q1[sign][1],Q3[HF][1]);

						double X_real_temp = Q1[sign][0];
						double Y_real_temp = Q1[sign][0];
						double Z_real_temp = Q3[HF][0];

						double X_imag_temp = Q1[sign][1];
						double Y_imag_temp = Q1[sign][1];
						double Z_imag_temp = Q3[HF][1];

						XY_real[trk][sign][HF]->Fill( XY_real_temp/(Q1_count[sign]*(Q1_count[sign]-1) ), (Q1_count[sign]*(Q1_count[sign]-1) ) );
						XZ_real[trk][sign][HF]->Fill( XZ_real_temp/(Q1_count[sign]*Q3_count[HF]), (Q1_count[sign]*Q3_count[HF]));
						YZ_real[trk][sign][HF]->Fill( YZ_real_temp/(Q1_count[sign]*Q3_count[HF]), (Q1_count[sign]*Q3_count[HF]));

						XY_imag[trk][sign][HF]->Fill( XY_imag_temp/(Q1_count[sign]*(Q1_count[sign]-1) ), (Q1_count[sign]*(Q1_count[sign]-1) ) );
						XZ_imag[trk][sign][HF]->Fill( XZ_imag_temp/(Q1_count[sign]*Q3_count[HF]), (Q1_count[sign]*Q3_count[HF]));
						YZ_imag[trk][sign][HF]->Fill( YZ_imag_temp/(Q1_count[sign]*Q3_count[HF]), (Q1_count[sign]*Q3_count[HF]));

						X_real[trk][sign][HF]->Fill( X_real_temp/Q1_count[sign], Q1_count[sign]);
						Y_real[trk][sign][HF]->Fill( Y_real_temp/Q1_count[sign], Q1_count[sign]);
						Z_real[trk][sign][HF]->Fill( Z_real_temp/Q3_count[HF], Q3_count[HF]);
						
						X_imag[trk][sign][HF]->Fill( X_imag_temp/Q1_count[sign], Q1_count[sign]);
						Y_imag[trk][sign][HF]->Fill( Y_imag_temp/Q1_count[sign], Q1_count[sign]);
						Z_imag[trk][sign][HF]->Fill( Z_imag_temp/Q3_count[HF], Q3_count[HF]);
					}	
					
					if( Q1_count[0] == 0 || Q1_count[1] == 0 || Q3_count[HF] == 0 ) continue;
					double Q_real = 0.;
					Q_real = get3Real(Q1[0][0]/Q1_count[0], Q1[1][0]/Q1_count[1], Q3[HF][0]/Q3_count[HF], Q1[0][1]/Q1_count[0], Q2[1][1]/Q1_count[1], Q3[HF][1]/Q3_count[HF]);
					QvsNtrk[trk][2][HF]->Fill( Q_real );

					double XY_real_temp = get2Real(Q1[0][0],Q1[1][0],Q1[0][1],Q1[1][1]);
					double XZ_real_temp = get2Real(Q1[0][0],Q3[HF][0],Q1[0][1],Q3[HF][1]);
					double YZ_real_temp = get2Real(Q1[1][0],Q3[HF][0],Q1[1][1],Q3[HF][1]);

					double XY_imag_temp = get2Imag(Q1[0][0],Q1[1][0],Q1[0][1],Q1[1][1]);
					double XZ_imag_temp = get2Imag(Q1[0][0],Q3[HF][0],Q1[0][1],Q3[HF][1]);
					double YZ_imag_temp = get2Imag(Q1[1][0],Q3[HF][0],Q1[1][1],Q3[HF][1]);

					double X_real_temp = Q1[0][0];
					double Y_real_temp = Q1[1][0];
					double Z_real_temp = Q3[HF][0];

					double X_imag_temp = Q1[0][1];
					double Y_imag_temp = Q1[1][1];
					double Z_imag_temp = Q3[HF][1];

					XY_real[trk][2][HF]->Fill( XY_real_temp/(Q1_count[0]*Q1_count[1]), (Q1_count[0]*Q1_count[1])  );
					XZ_real[trk][2][HF]->Fill( XZ_real_temp/(Q1_count[0]*Q3_count[HF]), (Q1_count[0]*Q3_count[HF]));
					YZ_real[trk][2][HF]->Fill( YZ_real_temp/(Q1_count[1]*Q3_count[HF]), (Q1_count[1]*Q3_count[HF]));

					XY_imag[trk][2][HF]->Fill( XY_imag_temp/(Q1_count[0]*Q1_count[1]), (Q1_count[0]*Q1_count[1])  );
					XZ_imag[trk][2][HF]->Fill( XZ_imag_temp/(Q1_count[0]*Q3_count[HF]), (Q1_count[0]*Q3_count[HF]));
					YZ_imag[trk][2][HF]->Fill( YZ_imag_temp/(Q1_count[1]*Q3_count[HF]), (Q1_count[1]*Q3_count[HF]));

					X_real[trk][2][HF]->Fill( X_real_temp/Q1_count[0], Q1_count[0]);
					Y_real[trk][2][HF]->Fill( Y_real_temp/Q1_count[1], Q1_count[1]);
					Z_real[trk][2][HF]->Fill( Z_real_temp/Q3_count[HF], Q3_count[HF]);
					
					X_imag[trk][2][HF]->Fill( X_imag_temp/Q1_count[0], Q1_count[0]);
					Y_imag[trk][2][HF]->Fill( Y_imag_temp/Q1_count[1], Q1_count[1]);
					Z_imag[trk][2][HF]->Fill( Z_imag_temp/Q3_count[HF], Q3_count[HF]);

				}				
			}

		}
		cout << "ievent: " << iEvent << endl;
	}

	TH1D* hist1[3][2];
	TH1D* hist2[3][2];

	for(int sign = 0; sign < 3; sign++){
		for(int HF = 0; HF < 2; HF ++ ){
			hist1[sign][HF] = new TH1D(Form("hist1_%d_%d",sign,HF),"test", nNtrkBins, Ntrkbins_fill);
			hist2[sign][HF] = new TH1D(Form("hist2_%d_%d",sign,HF),"test", nNtrkBins, Ntrkbins_fill);
		}
	}

	for(int trk = 0; trk < nNtrkBins; trk++){

		for(int sign = 0; sign < 3; sign++){
			for(int HF = 0; HF < 2; HF++){
				double Q = QvsNtrk[trk][sign][HF]->GetMean();
				double Q_error = QvsNtrk[trk][sign][HF]->GetMeanError();

				double XY_real_temp = 0.; double XZ_real_temp = 0.; double YZ_real_temp = 0.; double X_real_temp = 0.; double Y_real_temp = 0.; double Z_real_temp = 0.;
				double XY_imag_temp = 0.; double XZ_imag_temp = 0.; double YZ_imag_temp = 0.; double X_imag_temp = 0.; double Y_imag_temp = 0.; double Z_imag_temp = 0.;

				XY_real_temp = XY_real[trk][sign][HF]->GetMean();
				XZ_real_temp = XZ_real[trk][sign][HF]->GetMean();
				YZ_real_temp = YZ_real[trk][sign][HF]->GetMean();

				XY_imag_temp = XY_imag[trk][sign][HF]->GetMean();
				XZ_imag_temp = XZ_imag[trk][sign][HF]->GetMean();
				YZ_imag_temp = YZ_imag[trk][sign][HF]->GetMean();

				X_real_temp = X_real[trk][sign][HF]->GetMean();
				Y_real_temp = Y_real[trk][sign][HF]->GetMean();
				Z_real_temp = Z_real[trk][sign][HF]->GetMean();

				X_imag_temp = X_imag[trk][sign][HF]->GetMean();
				Y_imag_temp = Y_imag[trk][sign][HF]->GetMean();
				Z_imag_temp = Z_imag[trk][sign][HF]->GetMean();

				double term1_real = get2Real(XY_real_temp, Z_real_temp, XY_imag_temp, Z_imag_temp);
				double term2_real = get2Real(XZ_real_temp, Y_real_temp, XZ_imag_temp, Y_imag_temp);
				double term3_real = get2Real(YZ_real_temp, X_real_temp, YZ_imag_temp, X_imag_temp);
				double term4_real = get3Real(X_real_temp, Y_real_temp, Z_real_temp, X_imag_temp, Y_imag_temp, Z_imag_temp);
				double all_real = term1_real+term2_real+term3_real-2*term4_real;

				cout << "before correction is " << Q << endl;
				cout << "correction factor is: " << all_real << endl;
				cout << "after correction is " << Q - all_real << endl;
				cout << "-------- " << endl;

				hist1[sign][HF]->SetBinContent(trk+1, Q);
				hist1[sign][HF]->SetBinError(trk+1, Q_error);

				hist2[sign][HF]->SetBinContent(trk+1, Q - all_real );
				hist2[sign][HF]->SetBinError(trk+1, Q_error);

			}
		}

	}

	TCanvas* c1 = new TCanvas();
	c1->Divide(2,1);
	for(int i = 0; i < 2; i++){
		c1->cd(i+1);
		hist1[0][i]->SetMarkerColor(kBlack);
		hist1[0][i]->SetMarkerStyle(20);
		hist1[0][i]->GetYaxis()->SetRangeUser(-0.5, 0.5);
		hist1[0][i]->GetYaxis()->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3})");
		hist1[0][i]->GetXaxis()->SetTitle("N^{offline}_{trk}");
		hist1[0][i]->Draw("P");

		hist2[0][i]->SetMarkerColor(kRed);
		hist2[0][i]->SetMarkerStyle(20);
		hist2[0][i]->Draw("Psame");
	}


	TCanvas* c2 = new TCanvas();
	c2->Divide(2,1);
	for(int i = 0; i < 2; i++){
		c2->cd(i+1);
		hist1[1][i]->SetMarkerColor(kBlack);
		hist1[1][i]->SetMarkerStyle(20);
		hist1[1][i]->GetYaxis()->SetRangeUser(-0.5, 0.5);
		hist1[1][i]->GetYaxis()->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3})");
		hist1[1][i]->GetXaxis()->SetTitle("N^{offline}_{trk}");
		hist1[1][i]->Draw("P");

		hist2[1][i]->SetMarkerColor(kRed);
		hist2[1][i]->SetMarkerStyle(20);
		hist2[1][i]->Draw("Psame");
	}

	TCanvas* c3 = new TCanvas();
	c3->Divide(2,1);
	for(int i = 0; i < 2; i++){
		c3->cd(i+1);
		hist1[2][i]->SetMarkerColor(kBlack);
		hist1[2][i]->SetMarkerStyle(20);
		hist1[2][i]->GetYaxis()->SetRangeUser(-0.5, 0.5);
		hist1[2][i]->GetYaxis()->SetTitle("cos(#phi_{1}+#phi_{2}-2#phi_{3})");
		hist1[2][i]->GetXaxis()->SetTitle("N^{offline}_{trk}");
		hist1[2][i]->Draw("P");

		hist2[2][i]->SetMarkerColor(kRed);
		hist2[2][i]->SetMarkerStyle(20);
		hist2[2][i]->Draw("Psame");
	}


	return;



}