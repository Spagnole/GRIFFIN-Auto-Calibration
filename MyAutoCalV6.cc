#include <vector>
#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include "TChannel.h"
#include <tuple>

/*
PURPOSE:
Calibrate data in energy, refit if necessary.

REQUIRES:
Calibration sources .root file (Analysis trees) that's produced from selector file

// WORKFLOW:
// - in GRSISort, load this file .cc
// - MyAutoCal() initializes everything
// FirstAttempt() will make a rough guess at fitting the peaks with MyFits(), associate the fits with calibrations through FitCalibration()
// and lastly write the results to .root files.
// - Go detector by detector with FitCalibration() to check quality of fits, where residuals and fits between sources are printed onto
// a canvas.
// - ReFitPeak() can be used to fit any peaks there were not initially fit well by MyFit()
// - ReFitCalibration() can be used to separate the cal files written incase of drift between calibration runs. Format being feeding a 
// vector in the form of {0,1,...} where the numbers index the sources you want grouped together.
   - when ready BuildCalFile() to create .cal files for analysis tree
*/
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
string name_SpectraFile;
string name_GraphsFile;


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

//Below are a series of more legible ways to open files or start files
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
//void NewGraphs(string FileName, string NewFileName){
void NewGraphs(string NewFileName){
	//gSystem->CopyFile(FileName.c_str(), NewFileName.c_str() );
	GraphsFile->Cp(NewFileName.c_str());
	GraphsFile= TFile::Open(NewFileName.c_str(),"UPDATE");
}

void AddSource( string sourcename, string rootspectrum ){	
	Source[NSource] = sourcename;
	fSpectra[NSource] = TFile::Open(rootspectrum.c_str());
	NSource++;
}

void ReadLitEnergies(){	// Finds list of literature energies necessary for calibrations
	
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

void GetSpectra(){	// From input .root file made by selector, draw spectra for each detector
	for(int i = 0; i < NSource; i++){
		mat[i] = (TH2D*)fSpectra[i]->Get("gCA");
		for(int j = 0; j < NDet; j++){
			hist[i][j] = mat[i]->ProjectionY( Form("h%s_det%d", Source[i].c_str(), j+1 ) ,j+1 ,j+1 );
		}
	}
}


void GetCalPars(){ //PURPOSE: initial guess parameters used in selctor file is written to initial.cal  REQUIRES: input .root file from selector

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



void MyFits( int DetNum ){  // Fit peaks chosen detector

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
				if( fit_counter == 2 ){// if no good fit is found after 10 cycles, move on
					cout << "No good fit found...skipping data point" << endl;
					good_fit = false;
					break;
				}
			}			
			if( good_fit == true ){ //set fit if good fit is found
				cent = centr;
				P1[i][j][k]->GetFitFunction()->SetLineColor(kRed);
				//source_num[DetNum-1].push_back(i); is this line doing nothing?
				gData[i][j]->SetPoint(gData[i][j]->GetN(),cent, LitEnergy[i].at(k) );
				gData[i][j]->SetPointError(gData[i][j]->GetN()-1,err, 0.1 );
				fitnum_counter++;
			}
		}
		for(int k = 0; k < LitEnergy[i].size(); k++){ 
			pf[i][j][k]->Draw("same"); //Draws all fits onto histogram
			P1[i][j][k]->Draw("same");
		}
		hist[i][j]->Write(hist[i][j]->GetName(),TObject::kOverwrite);
		GraphsFile->cd(); gData[i][j]->Write(gData[i][j]->GetName(),TObject::kOverwrite);
		SpectraFile->cd();		
		//C->Close();
	}
}

void GetGraph(int DetNum, int SourceNumber){// helper function to find graphs

	int i = DetNum-1;
	int n = SourceNumber;
	gData[n][i] = (TGraphErrors*)GraphsFile->Get( Form("gDet%d_%s",DetNum,Source[n].c_str()) );
	gRes[n][i] = (TGraphErrors*)GraphsFile->Get( Form("gRes%d_%s",DetNum,Source[n].c_str()) );
}

void PrintFitResults(int DetNum){// helper function prints all peaks fitted so far

	GraphsFile->cd();	
	int i = DetNum-1;
	for(int n = 0; n < NSource; n++){
		GetGraph(DetNum,n);
		for(int j = 0; j < gData[n][i]->GetN(); j++){		
			cout << n << " " << Source[n]<< " " << " " << j << " " << gData[n][i]->GetPointY(j) << "\t" << gData[n][i]->GetPointX(j) << " +/- " << gData[n][i]->GetErrorX(j) << endl;
		} 	
	}
}

void CombineGraphs(int DetNum){// Merge all source graphs together

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
	gFit[i]->Write(gFit[i]->GetName(),TObject::kOverwrite);
}

void FitCalibration(int DetNum){ // Compares fitted peaks with known transition energies, find parameters to calibrate energies
/*
REQUIRES:
MyFits() has to be run on this detector at least once, so that literature energies can be compared against fitted energies
*/

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
		gRes[n][i]->Write(gRes[n][i]->GetName(),TObject::kOverwrite); // Writes residuals between calibrated energy and line found in calibration
	}
			
}

void WriteFit(int DetNum){// Fits are written to the graph file
	
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

void BuildCalFile(string calname){ // PURPOSE: Write .cal files for Analysis trees/selector scripts
// REQUIRES: User to be finished with all calibrations and other functions
//           User needs to input name for the resulting cal file
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
				gData[source_number][j]->Write(gData[source_number][j]->GetName(),TObject::kOverwrite);
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

void WriteFitsToText(int DetNum){// Creates .dat files to be read in elsewhere

	ofstream fitdata_out(Form("FitData/FitDet%d.dat",DetNum));
	int i = DetNum-1;
	for(int n = 0; n < NSource; n++){
		GetGraph(DetNum,n);
		for(int j = 0; j < gData[n][i]->GetN(); j++){		
			fitdata_out << n << " " << Source[n]<< " " << " " << j << " " << gData[n][i]->GetPointY(j) << "\t" << gData[n][i]->GetPointX(j) << " +/- " << gData[n][i]->GetErrorX(j) << endl;
		} 	
	}
}




void MyAutoCal(){// run this first to initialize everything

	name_SpectraFile = "Fits_92Rb";  //Name of the fit files, will be made if non-existent 
	name_GraphsFile = "Graphs_92Rb"; //Name of the graph files, will be made if non-existent 

	cout << "WELCOME to my energy calibration script which is used to produce an energy calibration\n";
	cout << "This script will read in lists of energies for different sources to be fitted\n";
	cout << "The energies are read in using the ReadLitEnergies() function\n";
	cout << "Make sure your order of sources matches the order of the input spectra so you find the correct peaks\n";
	cout << "The input root spectra were produced using a grsiproof selector using GetCharge() versus GetDetector()\nThe matrices are labelled gCA";
	cout << "After reading in the energies lists the script will then open the input spectrum using the OpenRootFiles() function\n";
	cout << "Then the GetSpectra() will convert the matrices into individual spectra for eash detector\n\n";
	cout << "Initial calibration parameters are obtained from the first input root file. These parameters are used to estimate the location of peaks in the Charge spectrum\n";
	cout << "- The peaks are fit using the MyFits(Detector Number) function. This will try to fit all the peaks from the list and is needed to try to perform a calibration\n";
	cout << "-  You can use the PrintFitResults(Detector Number) to see the results of the fits\n";
	cout << "- It is possible to use FirstAttempt() to use MyFits() over all detectors\n";
	cout << Form(  "The fits are written to a file currently called %s you can refit the same detector and it will overwrite the old fits\n",name_SpectraFile.c_str());
	cout << "-  Perform an energy calibration using the FitCalibration(Detector Number) function\n";
	cout << Form("  The graphs and energy calibration are written to %s which allows the same detector to be recalibrated any number of times\n",name_GraphsFile.c_str());
	cout << "-  Refit bad peaks using the  ReFitPeak(Detector Number) function\nThis allows you to refit any bad peak\n";
	cout << "- After refitting perform a new energy calibration using the FitCalibration(Detector Number) function\n";
	cout << "-  Use WriteCalToFile(your cal file name) function to write your calibration to a file\n\n";
	cout << "For a higher level overview, please read README_energy_calibration_and_nonlinearities.txt for overarching instructions"<< endl;


	//AddSource("152Eu","SumRuns/EnergyCalib_21662_sum.root"); //Selector created file that we need to read in
	//AddSource("60Co","SumRuns/EnergyCalib_21601_sum.root");
	//AddSource("133Ba","SumRuns/EnergyCalib_21658_sum.root");
	//AddSource("56Co","SumRuns/EnergyCalib_21598_21600_sum.root");
	AddSource("92Rb","EnergyCalib_Rb92_Sum.root");
	ReadLitEnergies(); //read in energy lists	
	//OpenRootFiles(); //open root files  make spectra	
	GetSpectra(); //makes spectra from root file
	GetCalPars(); //gets the calibration parameters
	SpectraFile = new TFile(Form("%s.root",name_SpectraFile.c_str()),"UPDATE");
	GraphsFile = new TFile(Form("%s.root",name_GraphsFile.c_str()),"UPDATE");
	cout<<name_SpectraFile<<endl;
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

inline bool exists_test (const std::string& name) {// Test if file exists
    ifstream f(name.c_str());
    return f.good();
}

void MakeGraphCombination(int DetNum, vector<int> SourceNum){
// Given a detector and combination of sources, create a new possible GraphsFile
// REQUIRES: vector of possible combination of sources ex. {1,2} for a combination
// of sources 1 and 2 to be calibrated together, whilst ignoring other sources for this detector
	string graphname;
	for(int n=0; n < SourceNum.size();n++){
		graphname += "_" + Source[SourceNum.at(n)];
	}

	string file_to_open = Form("%s_%s.root",name_GraphsFile.c_str(),graphname.c_str()); 
	if(exists_test (file_to_open) == 0){
		//NewGraphs(name_GraphsFile,file_to_open); // If no file exists for the source combination, make it based off the full file.
		NewGraphs(file_to_open);
	}

	OpenGraphs(file_to_open); // GraphsFile now points or has changed to the new file
}

void ReFitCalibration(int DetNum, vector<int> SourceNum){// Allows partitioning of calibration sources to individual graph files instead of all being combind into one.

	int i = DetNum-1;
	if( DetNum == 49 || DetNum == 50 || DetNum == 51 || DetNum == 52 ){
		cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
	}

	else{
		MakeGraphCombination(DetNum, SourceNum);
		GraphsFile->cd();
		cout << "PERFORMING FIT OF CALIBRATION PARAMETERS\n";
		
		TCanvas *C = new TCanvas();
		C->Divide(1,2);
		C->cd(1);
		//CombineGraphs(DetNum);
		//int i = DetNum-1;
		gFit[i] = new TGraphErrors();
		gFit[i]->SetName( Form("gDet%d_AllSources",DetNum) );
		gFit[i]->SetTitle( Form("gDet%d All Sources",DetNum) );
		gFit[i]->SetMarkerStyle(20);
		gFit[i]->SetMarkerColor(kRed);	
		for(int n = 0; n < SourceNum.size(); n++){
			GetGraph(DetNum,SourceNum.at(n));
			for(int j = 0; j < gData[SourceNum.at(n)][i]->GetN(); j++){
				gFit[i]->SetPoint(gFit[i]->GetN(), gData[SourceNum.at(n)][i]->GetPointX(j), gData[SourceNum.at(n)][i]->GetPointY(j) );
				gFit[i]->SetPointError(gFit[i]->GetN()-1, gData[SourceNum.at(n)][i]->GetErrorX(j), gData[SourceNum.at(n)][i]->GetErrorY(j) );
			} 	
		}
		gFit[i]->Write(gFit[i]->GetName(),TObject::kOverwrite);
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
		fit[DetNum-1]->Write(fit[DetNum-1]->GetName(),TObject::kOverwrite);
	}			
}

void FirstAttempt(){// Runs MyFits() for all detectors.
	for(int i = 0; i < 64; i++){
		int DetNum = i+1;
		if( DetNum != 49 && DetNum != 50 && DetNum != 51 && DetNum != 52 ){			
			MyFits(DetNum);
			FitCalibration(DetNum);
			WriteFit(DetNum);
		}
		else cout << "Do not perform fit for detector " << DetNum << "... Detector is empty\n";	
	}	
}

void AddPeak(int DetNum, int source_number, double NewLitEnergy){ 
/*
PURPOSE: Adds a new peak as defined by user to be fit

REQUIRES: literature energy to be added, source's index/number, detector number
*/

	GetGraph(DetNum, source_number);
	SpectraFile->cd();
	bool good_fit;
	int fit_counter=0;
	int i = source_number;
	int j = DetNum-1;
	int k = gData[i][j]->GetN();
	
	hist[source_number][j]->Draw();
	good_fit = true;
	cent = NewLitEnergy / cal_pars[j].at(1);
	low = cent - 10;
	upp = cent + 10;
	new_low = low;
	new_upp = upp;
	auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[source_number][j],cent,new_low,new_upp);
	pf[source_number][j][k] = pf_curr;
	P1[source_number][j][k] = P1_curr;
	cent = centr;
	while( err == 0 || width > 3. || width < 0.5 || chi2 > 1E5 || cent < new_low || cent > new_upp ){
		cent = NewLitEnergy / cal_pars[j].at(1);
		new_low = low + 8*(rand_par.Rndm()-0.5);
		new_upp = upp + 8*(rand_par.Rndm()-0.5);
		auto [centr, err, width, chi2, pf_curr, P1_curr] = FitResult(hist[source_number][j],cent,new_low,new_upp);
		pf[source_number][j][k] = pf_curr;
		P1[source_number][j][k] = P1_curr;
		cent = centr;
		fit_counter++;
		if( fit_counter == 10 ){
			cout << "No good fit found...skipping data point" << endl;
			good_fit = false;
			break;
		}
	}			
	if( good_fit == true ){
		cout << "GOOD FIT ADDING DATA POINT TO GRAPH" << endl;
		cent = centr;
		P1[source_number][j][k]->GetFitFunction()->SetLineColor(kRed);
		//source_num[DetNum-1].push_back(i);	
		GraphsFile->cd();	
		gData[source_number][j]->SetPoint(gData[source_number][j]->GetN(),cent, NewLitEnergy );
		gData[source_number][j]->SetPointError(gData[source_number][j]->GetN()-1,err, 0.1 );
		gData[source_number][j]->Write(gData[source_number][j]->GetName(),TObject::kOverwrite);
	}
		
	SpectraFile->cd();	hist[source_number][j]->Write(hist[source_number][j]->GetName(),TObject::kOverwrite);

}

void AddPeakToAllDetectors(){

	for(int i = 56; i < 64; i++){
		int DetNum = i+1;
		AddPeak(DetNum,3,122.06);
	}

}


void DrawGraphs(int DetNum){// Draw graphs created. Nothing special

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

void CalculateEnergy(int DetNum, double channel){

	double energy = channel * cal_pars[DetNum-1].at(1);
	cout << "Detector\tchannel\tenergy\n";
	cout << DetNum << "\t\t" << channel << "\t" << energy << endl;

}

void CalculateCharge(int DetNum, double energy){

	double charge = energy / cal_pars[DetNum-1].at(1);
	cout << "Detector\tchannel\tenergy\n";
	cout << DetNum << "\t\t" << charge << "\t" << energy << endl;

}

