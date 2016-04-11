#include "RiceStyle.h"

using namespace std;

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

double getAverage( TH1D* hist ){

	double integral2 = 0.;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		integral2 += temp*center;
	
	}
	double ave = integral2/( hist->GetEntries() );
	return ave;

}

double getSum( TH1D* hist ){

	double integral2 = 0.;
	for(int i = 0; i < hist->GetNbinsX(); i++ ){

		double temp = hist->GetBinContent(i+1);
		double center = hist->GetBinCenter(i+1);
		if(temp != 0 ){
			integral2 = integral2 + temp*center;
		}
	}

	double sum = integral2;
	return sum;
}

void plotNewTest(){
	
	TFile* file = new TFile("../rootfiles/CME_QvsdEta_pPb_EPOS_correction_test_v1.root");

	TH1D* QvsdEta[2];

	for(int i = 0; i < 2; i++){
		QvsdEta[i] = (TH1D*)file->Get(Form("ana/QvsdEta_%d",i) );
	}
	

	TH1D* XY[2];
	TH1D* XZ[2];
	TH1D* YZ[2];

	TH1D* XYcount;
	TH1D* XZcount;
	TH1D* YZcount;

	TH1D* X[2];
	TH1D* Y[2];
	TH1D* Z[2];

	TH1D* Xcount;
	TH1D* Ycount;
	TH1D* Zcount;

	for(int real = 0; real < 2; real++){

		XY[real] = (TH1D*) file->Get(Form("ana/XY_%d", real));
		XZ[real] = (TH1D*) file->Get(Form("ana/XZ_%d", real));
		YZ[real] = (TH1D*) file->Get(Form("ana/YZ_%d", real));

		X[real] = (TH1D*) file->Get(Form("ana/X_%d", real));
		Y[real] = (TH1D*) file->Get(Form("ana/Y_%d", real));
		Z[real] = (TH1D*) file->Get(Form("ana/Z_%d", real));

	}

	// XYcount = (TH1D*) file->Get("ana/XYcount");
	// XZcount = (TH1D*) file->Get("ana/XZcount");
	// YZcount = (TH1D*) file->Get("ana/YZcount");

	// Xcount = (TH1D*) file->Get("ana/Xcount");
	// Ycount = (TH1D*) file->Get("ana/Ycount");
	// Zcount = (TH1D*) file->Get("ana/Zcount");


	double XY_real = 0.; double XZ_real = 0.; double YZ_real = 0.; double X_real = 0.; double Y_real = 0.; double Z_real = 0.;
	double XY_imag = 0.; double XZ_imag = 0.; double YZ_imag = 0.; double X_imag = 0.; double Y_imag = 0.; double Z_imag = 0.;
	double XY_count = 0.; double XZ_count = 0.; double YZ_count = 0.; double X_count = 0.; double Y_count = 0.; double Z_count = 0.;

	cout << "getMean: " << XY[0]->GetMean() << endl;
	cout << "getAverage: " << getAverage(XY[0]) << endl;

	XY_real = XY[0]->GetMean(); XY_imag = XY[1]->GetMean();
	XZ_real = XZ[0]->GetMean(); XZ_imag = XZ[1]->GetMean();
	YZ_real = YZ[0]->GetMean(); YZ_imag = YZ[1]->GetMean();

	//XY_count = getSum( XYcount ); XZ_count = getSum( XZcount ); YZ_count = getSum( YZcount );

	X_real = X[0]->GetMean(); X_imag = X[1]->GetMean();
	Y_real = Y[0]->GetMean(); Y_imag = Y[1]->GetMean();
	Z_real = Z[0]->GetMean(); Z_imag = Z[1]->GetMean();

	//X_count = getSum( Xcount ); Y_count = getSum( Ycount ); Z_count = getSum( Zcount );

	// double term1_real = get2Real(XY_real/XY_count, Z_real/Z_count, XY_imag/XY_count, Z_imag/Z_count);
	// double term2_real = get2Real(XZ_real/XZ_count, Y_real/Y_count, XZ_imag/XZ_count, Y_imag/Y_count);
	// double term3_real = get2Real(YZ_real/YZ_count, X_real/X_count, YZ_imag/YZ_count, X_imag/X_count);
	// double term4_real = get3Real(X_real/X_count, Y_real/Y_count, Z_real/Z_count, X_imag/X_count, Y_imag/Y_count, Z_imag/Z_count);
	// double all_real = term1_real+term2_real+term3_real-2*term4_real;

	// double term1_imag = get2Imag(XY_real/XY_count, Z_real/Z_count, XY_imag/XY_count, Z_imag/Z_count);
	// double term2_imag = get2Imag(XZ_real/XZ_count, Y_real/Y_count, XZ_imag/XZ_count, Y_imag/Y_count);
	// double term3_imag = get2Imag(YZ_real/YZ_count, X_real/X_count, YZ_imag/YZ_count, X_imag/X_count);
	// double term4_imag = get3Imag(X_real/X_count, Y_real/Y_count, Z_real/Z_count, X_imag/X_count, Y_imag/Y_count, Z_imag/Z_count);
	// double all_imag = term1_imag+term2_imag+term3_imag-2*term4_imag;

	double term1_real = get2Real(XY_real, Z_real, XY_imag, Z_imag);
	double term2_real = get2Real(XZ_real, Y_real, XZ_imag, Y_imag);
	double term3_real = get2Real(YZ_real, X_real, YZ_imag, X_imag);
	double term4_real = get3Real(X_real, Y_real, Z_real, X_imag, Y_imag, Z_imag);
	double all_real = term1_real+term2_real+term3_real-2*term4_real;

	double term1_imag = get2Imag(XY_real, Z_real, XY_imag, Z_imag);
	double term2_imag = get2Imag(XZ_real, Y_real, XZ_imag, Y_imag);
	double term3_imag = get2Imag(YZ_real, X_real, YZ_imag, X_imag);
	double term4_imag = get3Imag(X_real, Y_real, Z_real, X_imag, Y_imag, Z_imag);
	double all_imag = term1_imag+term2_imag+term3_imag-2*term4_imag;

	double real_Qplusplus = QvsdEta[0]->GetMean();
	double imag_Qplusplus = QvsdEta[1]->GetMean();

	cout << "before correction real: " << real_Qplusplus << endl;
	cout << "after correction real: " << real_Qplusplus - all_real << endl;
	cout << "all_real: " << all_real << endl;
	cout << "before correction imag: " << imag_Qplusplus << endl;
	cout << "after correction imag: " << imag_Qplusplus - all_imag << endl;
	cout << "all_imag: " << all_imag << endl;

}