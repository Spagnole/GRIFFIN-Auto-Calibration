#include <vector>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include "TChannel.h"
#include <tuple>
using namespace std;

TRandom3 rand_par;

const int NDet = 64;
const int MaxSources=10;
int NSource;
string Source[MaxSources];// = {"60Co","152Eu","133Ba","56Co"};	
vector<double> LitEnergy[MaxSources];

TFile *fSpectra[MaxSources]; //read in spectra
TFile *SpectraFile; //writing spectra to a file
TFile *GraphsFile; //write graphs and fits to file
TH2D *mat[MaxSources];	TH1D *hist[MaxSources][NDet]; //read in matrices and make spectra	
TH2D *calmat[MaxSources];	TH1D *calhist[MaxSources][NDet];	
TChannel *mychannel; //calibration information
TGraphErrors *gData[MaxSources][NDet];	TGraphErrors *gRes[MaxSources][NDet];
TGraphErrors *gFit[NDet];
TF1 *fit[NDet];

//fitting information
double cent, err, chi2, width;
double low, upp, new_low, new_upp;
//fit information for calibration
vector<int> source_num[NDet]; vector<int> FitNumber[NDet];
//calibration parameters to read from initial cal file
vector<Float_t> cal_pars[NDet];
//calibration parameters to write to new file
double cal_a1[NDet], cal_a2[NDet], cal_a3[NDet];
//fitting functions
TPeakFitter *pf[MaxSources][NDet][100];
TRWPeak *P1[MaxSources][NDet][100];

void AddSource( string sourcename, string rootspectrum ){	
	Source[NSource] = sourcename;
	fSpectra[NSource] = TFile::Open(rootspectrum.c_str());
	NSource++;
}

void ReadLitEnergies(){	
	
	//LitEnergy[i].clear();
	double read_in;		
	ifstream file_energies[NSource];
	for(int i = 0; i < NSource; i++){
		LitEnergy[i].clear();
		file_energies[i].open( Form("LitEn_%s.dat",Source[i].c_str() ) );
		file_energies[i] >> read_in;
		while( !file_energies[i].eof() ){
			LitEnergy[i].push_back(read_in);
			file_energies[i] >> read_in;
		}
		file_energies[i].close();
	}
}

void PrintLitEnergies(int source_number){

	cout << "Printing list of energies for source " << Source[source_number].c_str() << endl;
	for(int i = 0; i < LitEnergy[source_number].size(); i++){
		cout << LitEnergy[source_number].at(i) << "\t";
	}
	cout << endl;
}

void GetSpectra(){	
	for(int i = 0; i < NSource; i++){
		mat[i] = (TH2D*)fSpectra[i]->Get("gCA");
		for(int j = 0; j < NDet; j++){
			hist[i][j] = mat[i]->ProjectionY( Form("h%s_det%d", Source[i].c_str(), j+1 ) ,j+1 ,j+1 );
		}
	}
}


void GetCalPars(){

	//TFile *AnalysisFile = TFile::Open("EnergyCalib_21601_000.root");	
	TChannel *chn = (TChannel*)fSpectra[0]->Get("Channel"); //This line is needed as without it the mychannel is a null ptr //I do not know why
	for(int i = 0; i < 16; i++){
		for(int j = 0; j < 4; j++){
	 		mychannel = TChannel::FindChannelByName(Form("GRG%02d%sN00A", i+1, TGriffin::GetColorFromNumber(j)));
	 		cal_pars[i*4+j] = mychannel->GetENGCoeff();
	 		cal_a1[i*4+j] = cal_pars[i*4+j].at(0);
	 		cal_a2[i*4+j] = cal_pars[i*4+j].at(1);
		}
	}
	TChannel::WriteCalFile("Initial.cal");
}

//here is the fitting procedure that is used in the MyFits() and ReFitPeak() finctions
tuple<double, double, double, double, TPeakFitter*, TRWPeak*> FitResult(TH1D *histogram, double cent_val, double low_bound, double up_bound){
	TPeakFitter *peak_fit = new TPeakFitter(low_bound,up_bound);
	TRWPeak *peak = new TRWPeak(cent_val);
	peak_fit->AddPeak(peak);
	peak_fit->Fit(histogram,"EQLS+");
	double new_cent_val = peak->GetFitFunction()->GetParameter(1);
	double error = peak->GetFitFunction()->GetParError(1);
	double width = peak->GetFitFunction()->GetParameter(2);
	double chi2 = peak->GetFitFunction()->GetParameter(0) / peak_fit->GetFitFunction()->GetNDF();
	return {new_cent_val, error, width, chi2, peak_fit, peak};
}



void MyFits( int DetNum ){

	SpectraFile->cd();
	bool good_fit;
	int fitnum_counter=0;
	for(int i = 0; i < NSource; i++){
		int j = DetNum-1;	
		GraphsFile->cd();	
		gData[i][j] = new TGraphErrors();
		gData[i][j]->SetName( Form("gDet%d_%s",DetNum,Source[i].c_str()) );
		gData[i][j]->SetTitle( Form("gDet%d_%s",DetNum,Source[i].c_str()) );
		gData[i][j]->SetMarkerStyle(20);
		gData[i][j]->SetMarkerColor(i+1);
		SpectraFile->cd();			
		TCanvas *C = new TCanvas();
		if( hist[i][j]->Integral() < 10000 ) continue; //skipping empty spectra
		for(int k = 0; k < LitEnergy[i].size(); k++){
			cout << "Fiting peak " << LitEnergy[i].at(k) << " From SOURCE " << Source[i] << " Detector " << DetNum << endl;
			int fit_counter = 0;			
			good_fit = true;
			cent = LitEnergy[i].at(k) / cal_pars[j].at(1);
			low = cent - 10;
			upp = cent + 10;
			new_low = low;
			new_upp = upp;
			auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[i][j],cent,new_low,new_upp);
			pf[i][j][k] = pf_curr;
			P1[i][j][k] = P1_curr;
			cent = centr;
			while( err == 0 || width > 3. || width < 0.5 || chi2 > 1E5 || cent < new_low || cent > new_upp ){
				cent = LitEnergy[i].at(k) / cal_pars[j].at(1);
				new_low = low + 8*(rand_par.Rndm()-0.5);
				new_upp = upp + 8*(rand_par.Rndm()-0.5);
				auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[i][j],cent,new_low,new_upp);
				pf[i][j][k] = pf_curr;
				P1[i][j][k] = P1_curr;
				cent = centr;
				fit_counter++;
				if( fit_counter == 10 ){
					cout << "No good fit found...skipping data point" << endl;
					good_fit = false;
					break;
				}
			}			
			if( good_fit == true ){
				cent = centr;
				P1[i][j][k]->GetFitFunction()->SetLineColor(kRed);
				source_num[DetNum-1].push_back(i);
				gData[i][j]->SetPoint(gData[i][j]->GetN(),cent, LitEnergy[i].at(k) );
				gData[i][j]->SetPointError(gData[i][j]->GetN()-1,err, 0.1 );
				fitnum_counter++;
			}
		}
		for(int k = 0; k < LitEnergy[i].size(); k++){ 
			pf[i][j][k]->Draw("same");
			P1[i][j][k]->Draw("same");
		}
		hist[i][j]->Write(hist[i][j]->GetName(),TObject::kOverwrite);
		GraphsFile->cd(); gData[i][j]->Write(gData[i][j]->GetName(),TObject::kOverwrite);
		SpectraFile->cd();		
		//C->Close();
	}
}

void GetGraph(int DetNum, int SourceNumber){

	int i = DetNum-1;
	int n = SourceNumber;
	gData[n][i] = (TGraphErrors*)GraphsFile->Get( Form("gDet%d_%s",DetNum,Source[n].c_str()) );
	gRes[n][i] = (TGraphErrors*)GraphsFile->Get( Form("gRes%d_%s",DetNum,Source[n].c_str()) );
}

void PrintFitResults(int DetNum){

	GraphsFile->cd();	
	int i = DetNum-1;
	for(int n = 0; n < NSource; n++){
		GetGraph(DetNum,n);
		for(int j = 0; j < gData[n][i]->GetN(); j++){		
			cout << n << " " << Source[n]<< " " << " " << j << " " << gData[n][i]->GetPointY(j) << "\t" << gData[n][i]->GetPointX(j) << " +/- " << gData[n][i]->GetErrorX(j) << endl;
		} 	
	}
}

void CombineGraphs(int DetNum){

	int i = DetNum-1;
	gFit[i] = new TGraphErrors();
	gFit[i]->SetName( Form("gDet%d_AllSources",DetNum) );
	gFit[i]->SetTitle( Form("gDet%d All Sources",DetNum) );
	gFit[i]->SetMarkerStyle(20);
	gFit[i]->SetMarkerColor(kRed);	
	for(int n = 0; n < NSource; n++){
		GetGraph(DetNum,n);
		for(int j = 0; j < gData[n][i]->GetN(); j++){
			gFit[i]->SetPoint(gFit[i]->GetN(), gData[n][i]->GetPointX(j), gData[n][i]->GetPointY(j) );
			gFit[i]->SetPointError(gFit[i]->GetN()-1, gData[n][i]->GetErrorX(j), gData[n][i]->GetErrorY(j) );
		} 	
	}

}

void FitCalibration(int DetNum){

	GraphsFile->cd();
	cout << "PERFORMING FIT OF CALIBRATION PARAMETERS\n";
	int i = DetNum-1;
	cal_a1[i] = 0.;	cal_a2[i] = 1.;	cal_a3[i] = 0.;	
	if( DetNum == 49 || DetNum == 50 || DetNum == 51 || DetNum == 52 ){
		cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
	}
	
	TCanvas *C = new TCanvas();
	C->Divide(1,2);
	C->cd(1);
	
	CombineGraphs(DetNum);
	gFit[i]->Draw("AP");
	fit[i] = new TF1( Form("Fit_Det%d", DetNum), "[0]+[1]*x+[2]*x*x", 0, 4000);
	fit[i]->SetParameters(cal_a1[i], cal_a2[i], cal_a3[i] );
	fit[i]->FixParameter(2,0.0);
	gFit[i]->Fit( Form("Fit_Det%d", DetNum), "Q" );
	cal_a1[i] = fit[i]->GetParameter(0);
	cal_a2[i] = fit[i]->GetParameter(1);
	
	C->cd(2);
	TH1F *hFrame = C->DrawFrame(0,-10,5000,10);
	for(int n = 0; n < NSource; n++){
		gRes[n][i] = new TGraphErrors();
		gRes[n][i]->SetName( Form("gRes%d_%s",DetNum,Source[n].c_str()) );
		gRes[n][i]->SetMarkerStyle(20);
		gRes[n][i]->SetMarkerColor(n+1);
		for(int j = 0; j < gData[n][i]->GetN(); j++){
			double channel =  gData[n][i]->GetPointX(j);
			double LitEnergy = gData[n][i]->GetPointY(j);
			double calc_energy = cal_a1[i] + cal_a2[i]*channel + cal_a3[i]*channel*channel;
			gRes[n][i]->SetPoint(j, LitEnergy, LitEnergy-calc_energy );
			gRes[n][i]->SetPointError(j, 0.05, gData[n][i]->GetErrorX(j));		
		}
		gRes[n][i]->Draw("P");
		gRes[n][i]->Write(gRes[n][i]->GetName(),TObject::kOverwrite);
	}		
}

void WriteFit(int DetNum){
	
	if( DetNum == 49 || DetNum == 50 || DetNum == 51 || DetNum == 52 ){
		cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
		return;
	}
	else{
		GraphsFile->cd();
		gFit[DetNum-1]->Write(gFit[DetNum-1]->GetName(),TObject::kOverwrite);
		fit[DetNum-1]->Write(fit[DetNum-1]->GetName(),TObject::kOverwrite);
	}
}

//this will write all calibration parameters to a new calibration files
/*
void WriteCalToFile(string calname){

	TF1 *fTemp; //
	for(int i = 0; i < 16; i++){
		for(int j = 0; j < 4; j++){
			cout << i*4+j+1 << endl;
			if( GraphsFile->GetListOfKeys()->Contains(Form("Fit_Det%d", i*4+j+1) ) ){
				fTemp = (TF1*)GraphsFile->Get(Form("Fit_Det%d", i*4+j+1));
		 		mychannel = TChannel::FindChannelByName(Form("GRG%02d%sN00A", i+1, TGriffin::GetColorFromNumber(j)));
		 		mychannel->DestroyENGCal();
		 		mychannel->AddENGCoefficient(fTemp->GetParameter(0));
		 		mychannel->AddENGCoefficient(fTemp->GetParameter(1));
		 		mychannel->Print();
	 		}
	 		else if( !GraphsFile->GetListOfKeys()->Contains(Form("Fit_Det%d", i*4+j+1) ) ){
	 			cout << "WARNING!!!!!!\n\n\nNo calibration exists for detector " << i*4+j+1 << endl;
	 			cout << "USING ORIGINAL CALIBRATION??? Double check not empty\n";	 		
	 		}
		}
	}
	TChannel::WriteCalFile(calname.c_str());
}
*/ //Should this be removed?

//this will write new calibration parameters for a specific detector to a new file
void WriteChannelToCalFile(int DetNum){

	TF1 *fTemp; //
	int i = (DetNum-1)/4;
	int j = (DetNum-1)%4;
	cout << DetNum << "\t" << i << "\t" << j << endl;	
	if( GraphsFile->GetListOfKeys()->Contains(Form("Fit_Det%d", DetNum) ) ){
		fTemp = (TF1*)GraphsFile->Get(Form("Fit_Det%d", DetNum));
		mychannel = TChannel::FindChannelByName(Form("GRG%02d%sN00A", i+1, TGriffin::GetColorFromNumber(j)));
		mychannel->DestroyENGCal();
 		mychannel->AddENGCoefficient(fTemp->GetParameter(0));
 		mychannel->AddENGCoefficient(fTemp->GetParameter(1));
		mychannel->Print();
	}
}

void BuildCalFile(string calname){
	for(int i = 0; i < NDet; i++){
		int DetNum = i+1;
		WriteChannelToCalFile(DetNum);	
	}
	TChannel::WriteCalFile(calname.c_str());
}

//This allows the user to re-fit any bad peaks to improve the calibration
void ReFitPeak(int DetNum){
		
	int j = DetNum-1;
	bool good_fit = false;
	int source_number;
	int fitnum;
	PrintFitResults(DetNum);
	SpectraFile->cd();
	cout << "Would you like to re-fit a peak? Energy (y/Y) to refit a peak or any other key to quit ReFitPeak()\n" << endl;
	string answer;
	cin >> answer;
	if( answer == "y" ||  answer ==  "Y" ){
		TCanvas *C = new TCanvas();	
		C->Update();	
		cout << "What is the source number of the peak you wish to fit?" << endl;
		cin >> source_number;
		if( source_number > NSource){
			cout << "No such source number";
			return;
		}
		hist[source_number][j]->Draw();	
		hist[source_number][j]->Draw("histsame");
		C->Update();
		cout << "What is the fit number of the peak to be refit?" << endl;
		cin >> fitnum;
		int k = fitnum;
		if ( fitnum > gData[source_number][j]->GetN()-1){
			cout << "Calibration point does not exist. Do you wish to add a new data point?" << endl;
			cin >> answer;
			if( answer == "y" ||  answer ==  "Y" ){
				cout << "What is the literature energy of the transition?" << endl;
				PrintLitEnergies(source_number);
				double litEn;
				cin >> litEn;
				gData[source_number][j]->SetPoint(fitnum , litEn, litEn );		
			}
			else{
				 cout<<"Data point not added"<<endl;
			}
		}
		cout << "Performing refit of data\ngive new range";
		cout << "give lower range = "<< endl;
		cin >> new_low; 
		cout << "give upper range = "<< endl;
		cin >> new_upp;
		C->Update();
		///
		cent = (new_low+new_upp)/2;
		auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[source_number][j],cent,new_low,new_upp);
		pf[source_number][j][k] = pf_curr;
		P1[source_number][j][k] = P1_curr;
		cent = centr;
		P1[source_number][j][k]->GetFitFunction()->SetLineColor(kGreen+2);		///
		//FitSinglePeak_pegs( (TH1F*)hist[source_number][j], new_low, new_upp, "E", cent, err, chi2, width );
		C->Update();
		cout << "ARE YOU HAPPY WITH THIS FIT? (answer y or Y for yes or any other key if you are not satisfied with the fit" << endl;
		cin >> answer;
		if( answer == "y" ||  answer ==  "Y" ) good_fit = true;
		while( good_fit == false ){
			cout << "Performing refit of data\ngive new range";
			cout << "give lower range = "<< endl;
			cin >> new_low; 
			cout << "give upper range = "<< endl;
			cin >> new_upp;
			cent = (new_low+new_upp)/2;
			auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[source_number][j],cent,new_low,new_upp);
			pf[source_number][j][k] = pf_curr;
			P1[source_number][j][k] = P1_curr;
			cent = centr;	
			P1[source_number][j][k]->GetFitFunction()->SetLineColor(kGreen+2);	///
			//FitSinglePeak_pegs( (TH1F*)hist[source_number][j], new_low, new_upp, "E", cent, err, chi2, width );
			hist[source_number][j]->GetXaxis()->SetRangeUser(new_low-50,new_upp+50 );
			C->Update();
			cout << "ARE YOU HAPPY WITH THIS FIT? (answer y or Y for yes or skip to discard refit or erase to remove data point feom callbration" << endl;
			cin >> answer;
			if( answer == "y" ||  answer ==  "Y" ){
				good_fit = true;
				break;
			}
			else if(answer == "skip"){
				cout << "Leaving original fit as is...No action performed\n";
				break;
			
			}
			else if(answer == "erase"){
				cout << "Erasing data point from graph...The data point will not be used\n";
				GraphsFile->cd();
				gData[source_number][j]->RemovePoint(fitnum);
				break;
			}
		}
		if( good_fit == true ){
			GraphsFile->cd();
			cout << source_number << "\t" << fitnum << " " << gData[source_number][j]->GetPointX(fitnum) << " " << gData[source_number][j]->GetPointY(fitnum) << endl;
			gData[source_number][j]->SetPoint(fitnum , cent, gData[source_number][j]->GetPointY(fitnum) );
			gData[source_number][j]->SetPointError(fitnum , err, 0.1 );
			cout << source_number << "\t" << fitnum << " " << gData[source_number][j]->GetPointX(fitnum) << " " << gData[source_number][j]->GetPointY(fitnum) << endl;
			gData[source_number][j]->Write(gData[source_number][j]->GetName(),TObject::kOverwrite);
		}
	}
}

void WriteFitsToText(int DetNum){

	ofstream fitdata_out(Form("FitData/FitDet%d.dat",DetNum));
	int i = DetNum-1;
	for(int n = 0; n < NSource; n++){
		GetGraph(DetNum,n);
		for(int j = 0; j < gData[n][i]->GetN(); j++){		
			fitdata_out << n << " " << Source[n]<< " " << " " << j << " " << gData[n][i]->GetPointY(j) << "\t" << gData[n][i]->GetPointX(j) << " +/- " << gData[n][i]->GetErrorX(j) << endl;
		} 	
	}
}




void MyAutoCal(){

	cout << "Welcome to my energy calibration script which is used to produce an energy calibration\n";
	cout << "This script will read in lists of energies for different sources to be fitted\n";
	cout << "The energies are read in using the ReadLitEnergies() function\n";
	cout << "Make sure your order of sources matches the order of the input spectra so you find the correct peaks\n";
	cout << "The input root spectra were produced using a grsiproof selector using GetCharge() versus GetDetector()\nThe matrices are labelled gCA";
	cout << "After reading in the energies lists the script will then open the input spectrum using the OpenRootFiles() function\n";
	cout << "Then the GetSpectra() will convert the matrices into individual spectra for eash detector\n";
	cout << "Initial calibration parameters are obtained from the first input root file. These parameters are used to estimate the location of peaks in the Charge spectrum\n";
	cout << "The peaks are fit using the MyFits(Detector Number) function. This will try to fit all the peaks from the list and is needed to try to perform a calibration\n";
	cout << "You can use the PrintFitResults(Detector Number) to see the results of the fits\n";
	cout << "The fits are written to a file currently called Fits.root you can refit the same detector and it will overwrite the old fits\n";
	cout << "Perform an energy calibration using the FitCalibration(Detector Number) function\n";
	cout << "The graphs and energy calibration are written to a root file which allows the same detector to be recalibrated any number of times\n";
	cout << "Refit bad peaks using the  ReFitPeak(Detector Number) function\nThis allows you to refit any bad peak\n";
	cout << "After refitting perform a new energy calibration using the FitCalibration(Detector Number) function\n";
	cout << "Use WriteCalToFile(your cal file name) function to write your calibration to a file\n";



	//ReadLitEnergies(); //read in energy lists	
	//OpenRootFiles(); //open root files  make spectra	
	//GetSpectra(); //makes spectra from root file
	//GetCalPars(); //gets the calibration parameters
	//SpectraFile = new TFile("PietroFits.root","UPDATE");
	//GraphsFile = new TFile("PietroGraphs.root","UPDATE");
/*	for(int i = 0; i < 64; i++){
		MyFits(i+1);
		FitCalibration(i+1);
		WriteGraph(i+1);
	}
	WriteCalToFile("test.cal");*/
	//MyFits(10);
	//PrintFitResults(10);
	//FitCalibration(10);
	//FitCalibration(10);
	//	MyFits(i+1);
	//	FitCalibration(i+1);
	//}
	//WriteCalToFile("test.cal");

}

void ReFitCalibration(int DetNum, vector<int> SourceNum){

	GraphsFile->cd();
	cout << "PERFORMING FIT OF CALIBRATION PARAMETERS\n";
	int i = DetNum-1;
	if( DetNum == 49 || DetNum == 50 || DetNum == 51 || DetNum == 52 ){
		cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
	}
	
	string graphname;
	for(int n=0; n < SourceNum.size();n++){
		graphname += "_" + Source[SourceNum.at(n)];
	}
	
	TCanvas *C = new TCanvas();
	C->Divide(1,2);
	C->cd(1);
	CombineGraphs(DetNum);
	gFit[i]->Draw("AP");
	for(int n = 0; n < SourceNum.size(); n++){
		gData[SourceNum.at(n)][i]->Draw("P");
	}
	fit[i] = new TF1( Form("Fit_Det%d", DetNum), "[0]+[1]*x+[2]*x*x", 0, 4000);
	fit[i]->SetParameters(cal_a1[i], cal_a2[i], cal_a3[i] );
	fit[i]->FixParameter(2,0.0);
	gFit[i]->Fit( Form("Fit_Det%d", DetNum), "Q" );
	cal_a1[i] = fit[i]->GetParameter(0);
	cal_a2[i] = fit[i]->GetParameter(1);	
	C->cd(2);
	TH1F *hFrame = C->DrawFrame(0,-10,5000,10);
	for(int n = 0; n < SourceNum.size(); n++){
		gRes[SourceNum.at(n)][i] = new TGraphErrors();
		gRes[SourceNum.at(n)][i]->SetName( Form("gRes%d%s",DetNum,Source[i].c_str()) );
		gRes[SourceNum.at(n)][i]->SetMarkerStyle(20);
		gRes[SourceNum.at(n)][i]->SetMarkerColor(n+1);
		for(int j = 0; j < gData[SourceNum.at(n)][i]->GetN(); j++){
			double channel =  gData[SourceNum.at(n)][i]->GetPointX(j);
			double LitEnergy = gData[SourceNum.at(n)][i]->GetPointY(j);
			double calc_energy = cal_a1[i] + cal_a2[i]*channel + cal_a3[i]*channel*channel;
			gRes[SourceNum.at(n)][i]->SetPoint(j, LitEnergy, LitEnergy-calc_energy );
			gRes[SourceNum.at(n)][i]->SetPointError(j, 0.05, gData[SourceNum.at(n)][i]->GetErrorX(j));		
		}
		gRes[SourceNum.at(n)][i]->Draw("P");
	}	
}

void FirstAttempt(){
	for(int i = 48; i < 64; i++){
		int DetNum = i+1;
		if( DetNum != 49 && DetNum != 50 && DetNum != 51 && DetNum != 52 ){			
			MyFits(DetNum);
			FitCalibration(DetNum);
			WriteFit(DetNum);
		}
		else cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
	}	
}

void OpenFits(string FileName){
	SpectraFile = TFile::Open(FileName.c_str(),"UPDATE");
}
void OpenGraphs(string FileName){
	GraphsFile = TFile::Open(FileName.c_str(),"UPDATE");
}
void NewFits(string FileName, string NewFileName){
	gSystem->CopyFile(FileName.c_str(), NewFileName.c_str() );
	SpectraFile= TFile::Open(NewFileName.c_str(),"UPDATE");
}
void NewGraphs(string FileName, string NewFileName){
	gSystem->CopyFile(FileName.c_str(), NewFileName.c_str() );
	GraphsFile= TFile::Open(NewFileName.c_str(),"UPDATE");
}

void DrawGraphs(int DetNum){

	TCanvas *C = new TCanvas();
	C->Divide(1,2);
	C->cd(1);
	CombineGraphs(DetNum);
	gFit[DetNum-1]->Draw("AP");
	for(int i = 0; i < NSource; i++){
		gData[i][DetNum-1]->Draw("P");	
	}
	C->cd(2);
	TH1F *hFrame = C->DrawFrame(0,-10,5000,10);
	for(int i = 0; i < NSource; i++){		
		gRes[i][DetNum-1]->Draw("P");	
	}
}

