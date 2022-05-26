#include "MyAutoCalV6.cc"


TGraphErrors *gNonLin[MaxSources][NDet];	
TGraphErrors *gNonLinFull[NDet];

vector<tuple<double, double> > NonLinearityPoint[NDet];
TFile *CalFits;
TFile *NonLinGraphs;

TPeakFitter *pf_NonLin[MaxSources][NDet][50];
TRWPeak *Peak_NonLin[MaxSources][NDet][50];

void GetCalSpectra(){	
	for(int i = 0; i < NSource; i++){
		calmat[i] = (TH2D*)fSpectra[i]->Get("gEA_NoNonLin");
		for(int j = 0; j < NDet; j++){
			calhist[i][j] = calmat[i]->ProjectionY( Form("h%s_det%d", Source[i].c_str(), j+1 ) ,j+1 ,j+1 );
			//cout << Source[i] << "\tDetector " << j+1 << "\t counts = " << hist[i][j]->Integral() << endl;
		}
	}
}


void FitNonLinearities( int DetNum ){

	SpectraFile->cd();
	bool good_fit;
	int fitnum_counter=0;
	for(int i = 0; i < NSource; i++){
		int j = DetNum-1;	
		GraphsFile->cd();	
		gNonLin[i][j] = new TGraphErrors();
		gNonLin[i][j]->SetName( Form("gNonLin%d_%s",DetNum,Source[i].c_str()) );
		gNonLin[i][j]->SetTitle( Form("gNonLin%d_%s",DetNum,Source[i].c_str()) );
		gNonLin[i][j]->SetMarkerStyle(20);
		gNonLin[i][j]->SetMarkerColor(i+1);
		SpectraFile->cd();			
		TCanvas *C = new TCanvas();
		if( calhist[i][j]->Integral() < 10000 ) continue; //skipping empty spectra
		for(int k = 0; k < LitEnergy[i].size(); k++){
			cout << "Fiting peak " << LitEnergy[i].at(k) << " From SOURCE " << Source[i] << " Detector " << DetNum << endl;
			int fit_counter = 0;			
			good_fit = true;
			cent = LitEnergy[i].at(k);
			low = cent - 10;
			upp = cent + 10;
			new_low = low;
			new_upp = upp;
			auto [centr, err, width, chi2, pf_NonLin_curr, Peak_NonLin_curr] = FitResult(calhist[i][j],cent,new_low,new_upp);
			pf_NonLin[i][j][k] = pf_NonLin_curr;
			Peak_NonLin[i][j][k] = Peak_NonLin_curr;
			cent = centr;
			if( err == 0 || width > 3. || width < 0.5 || chi2 > 1E5 || cent < new_low || cent > new_upp ) good_fit = false;
			/*while( err == 0 || width > 3. || width < 0.5 || chi2 > 1E5 || cent < new_low || cent > new_upp ){
				cent = LitEnergy[i].at(k);
				new_low = low + 8*(rand_par.Rndm()-0.5);
				new_upp = upp + 8*(rand_par.Rndm()-0.5);
				auto [centr, err, width, chi2, pf_NonLin_curr, Peak_NonLin_curr] = FitResult(calhist[i][j],cent,new_low,new_upp);
				pf_NonLin[i][j][k] = pf_NonLin_curr;
				Peak_NonLin[i][j][k] = Peak_NonLin_curr;
				cent = centr;
				fit_counter++;
				if( fit_counter == 10 ){
					cout << "No good fit found...skipping data point" << endl;
					good_fit = false;
					break;
				}
			}	*/		
			if( good_fit == true ){
				cent = centr;
				Peak_NonLin[i][j][k]->GetFitFunction()->SetLineColor(kRed);
				source_num[DetNum-1].push_back(i);
				gNonLin[i][j]->SetPoint(gNonLin[i][j]->GetN(),LitEnergy[i].at(k), cent - LitEnergy[i].at(k) );
				gNonLin[i][j]->SetPointError(gNonLin[i][j]->GetN()-1,0.1, err );
				fitnum_counter++;
			}
		}	
		for(int k = 0; k < LitEnergy[i].size(); k++){ 
			pf_NonLin[i][j][k]->Draw("same");
			Peak_NonLin[i][j][k]->Draw("same");
		}
		calhist[i][j]->Write(calhist[i][j]->GetName(),TObject::kOverwrite);
		GraphsFile->cd(); gNonLin[i][j]->Write(gNonLin[i][j]->GetName(),TObject::kOverwrite);
		SpectraFile->cd();	
	}
}




void SortNonLinearities( int DetNum ){

	//read non linearities from each source
	for(int i = 0; i < NSource; i++){
		gNonLin[i][DetNum-1] = (TGraphErrors*)GraphsFile->Get( Form("gNonLin%d_%s",DetNum,Source[i].c_str()) );
		for(int j = 0; j < gNonLin[i][DetNum-1]->GetN(); j++){
			cout << i << " " << Source[i] << " " << gNonLin[i][DetNum-1]->GetPointX(j) << " " << gNonLin[i][DetNum-1]->GetPointY(j) << endl;
			NonLinearityPoint[DetNum-1].push_back(make_tuple(gNonLin[i][DetNum-1]->GetPointX(j), gNonLin[i][DetNum-1]->GetPointY(j)));
		}
	}
	
	//sort the non linearities
	sort(NonLinearityPoint[DetNum-1].begin(), NonLinearityPoint[DetNum-1].end());
	
	//make the graph
	GraphsFile->cd();	
	gNonLinFull[NDet-1] = new TGraphErrors();
	gNonLinFull[NDet-1]->SetName( Form("gNonLinFull%d",DetNum ) );
	gNonLinFull[NDet-1]->SetTitle( Form("gNonLinFull%d",DetNum ) );
	gNonLinFull[NDet-1]->SetMarkerStyle(20);
	//now that they are sorted
	gNonLinFull[NDet-1]->SetPoint(gNonLinFull[NDet-1]->GetN(),  10., 0.0 );
	for(int k = 0; k <  NonLinearityPoint[DetNum-1].size(); k++){
		//cout << NonLinearityPoint[DetNum-1].at(k) << endl;	
		cout << get<0>(NonLinearityPoint[DetNum-1][k]) << " " << get<1>(NonLinearityPoint[DetNum-1][k]) << endl;
		gNonLinFull[NDet-1]->SetPoint(gNonLinFull[NDet-1]->GetN(),  get<0>(NonLinearityPoint[DetNum-1][k]), get<1>(NonLinearityPoint[DetNum-1][k]) );
	}
	gNonLinFull[NDet-1]->SetPoint(gNonLinFull[NDet-1]->GetN(),  10000., 0.0 );
	gNonLinFull[NDet-1]->Draw("ALP");
	gNonLinFull[NDet-1]->Write(gNonLinFull[NDet-1]->GetName(), TObject::kOverwrite);
}

void MakeCalChannel(int DetNum){

	TF1 *fTemp;
	int i = (DetNum-1)/4;
	int j = (DetNum-1)%4;
	cout << DetNum << "\t" << i << "\t" << j << endl;	
	if( GraphsFile->GetListOfKeys()->Contains(Form("gNonLinFull%d", DetNum) ) ){
		fTemp = (TF1*)GraphsFile->Get(Form("gNonLinFull%d", DetNum));
		mychannel = TChannel::FindChannelByName(Form("GRG%02d%sN00A", i+1, TGriffin::GetColorFromNumber(j)));
		//cal_pars[i*4+j] = mychannel->GetENGCoeff();
		//mychannel->DestroyENGCal();
		//mychannel->DestroyEnergyNonlinearity();
		//mychannel->DestroyCTCal();
		//mychannel->AddEnergyNonlinearityPoint(10, 0.0 );	
		for(int k = 0; k < gNonLinFull[NDet-1]->GetN(); k++){
			mychannel->AddEnergyNonlinearityPoint(gNonLinFull[NDet-1]->GetPointX(k), gNonLinFull[NDet-1]->GetPointY(k) );		
		}
		//mychannel->AddEnergyNonlinearityPoint(10000., 0.0 );	
		mychannel->Print();
	}
	//cout <<  cal_pars[i*4+j].at(0) << "\t" <<  cal_pars[i*4+j].at(1) << endl;
	//TChannel::WriteCalFile(calname.c_str());

}

void MakeCalFile(string calfile){

	for(int i = 0; i < 2; i++){
		MakeCalChannel(i+1);	
	}
	TChannel::WriteCalFile(calfile.c_str());
}

void GetNonLinearitiesSourceGraph(int DetNum, int SourceNumber){

	int i = DetNum-1;
	int n = SourceNumber;
	gNonLin[n][i] = (TGraphErrors*)GraphsFile->Get( Form("gNonLin%d_%s",DetNum,Source[n].c_str()) );
}

void PrintNonlinearities(int DetNum){

	GraphsFile->cd();	
	int i = DetNum-1;
	for(int n = 0; n < NSource; n++){
		GetNonLinearitiesSourceGraph(DetNum,n);
		for(int j = 0; j < gNonLin[n][i]->GetN(); j++){		
			cout << n << " " << Source[n]<< " " << " " << j << " " << gNonLin[n][i]->GetPointX(j) << "\t" << gNonLin[n][i]->GetPointY(j) << endl;
		} 	
	}
}

void SetUp(){

	//AddSource("152Eu","SumRuns/EnergyCalib_21662_sum.root");
	//AddSource("60Co","SumRuns/EnergyCalib_21601_sum.root");
	//AddSource("133Ba","SumRuns/EnergyCalib_21658_sum.root");
	AddSource("56Co","SumRuns/EnergyCalib_21598_21600_sum.root");
	ReadLitEnergies(); //read in energy lists	
	//OpenRootFiles(); //open root files  make spectra
	GetCalSpectra();	
	GetCalPars();
	OpenFits("CalSpectra-pietro.root");
	OpenGraphs("NonLinearGraphs-pietro.root");
	for(int i = 49; i < 64; i++){
		int DetNum=i+1;
		if( DetNum == 49 || DetNum == 50 || DetNum == 51 || DetNum == 52 ) continue;
		FitNonLinearities(i+1);
		SortNonLinearities(i+1);
	}	
}

void ReFitNonLinearities(int DetNum){

	int j = DetNum-1;
	bool good_fit = false;
	int source_number;
	int fitnum;
	PrintNonlinearities(DetNum);
	SpectraFile->cd();
	cout << "Would you like to re-fit a peak?" << endl;
	string answer;
	cin >> answer;
	if( answer == "y" ||  answer ==  "Y" ){
		TCanvas *C = new TCanvas();	
		cout << "What is the source number of the peak you wish to fit?" << endl;
		cin >> source_number;
		if( source_number > NSource){
			cout << "No such source number";
			return;
		}
		calhist[source_number][j]->Draw();	
		calhist[source_number][j]->Draw("histsame");
		C->Update();
		cout << "What is the fit number of the peak to be refit?" << endl;
		cin >> fitnum;
		int k = fitnum;
		if ( fitnum > gNonLin[source_number][j]->GetN()-1){
			cout << "Calibration point does not exist. Do you wish to add a new data point?" << endl;
			cin >> answer;
			if( answer == "y" ||  answer ==  "Y" ){
				cout << "What is the literature energy?" << endl;
				double litEn;
				cin >> litEn;
				gNonLin[source_number][j]->SetPoint(fitnum , litEn, 0.0 );		
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
		cent = (new_low+new_upp)/2;
		auto [centr, err, width, chi2, pf_NonLin_curr, Peak_NonLin_curr] = FitResult(calhist[source_number][j],cent,new_low,new_upp);
		pf_NonLin[source_number][j][k] = pf_NonLin_curr;
		Peak_NonLin[source_number][j][k] = Peak_NonLin_curr;
		cent = centr;
		Peak_NonLin[source_number][j][k]->GetFitFunction()->SetLineColor(kGreen+2);		///
		calhist[source_number][j]->GetXaxis()->SetRangeUser(new_low-50,new_upp+50 );
		C->Update();
		cout << "ARE YOU HAPPY WITH THIS FIT? (answer y or Y for yes or skip to discard refit" << endl;
		cin >> answer;
		if( answer == "y" ||  answer ==  "Y" ) good_fit = true;
		while( good_fit == false ){
			cout << "Performing refit of data\ngive new range";
			cout << "give lower range = "<< endl;
			cin >> new_low; 
			cout << "give upper range = "<< endl;
			cin >> new_upp;
			cent = (new_low+new_upp)/2;
			auto [centr, err, width, chi2, pf_NonLin_curr, Peak_NonLin_curr] = FitResult(calhist[source_number][j],cent,new_low,new_upp);
			pf_NonLin[source_number][j][k] = pf_NonLin_curr;
			Peak_NonLin[source_number][j][k] = Peak_NonLin_curr;
			cent = centr;
			Peak_NonLin[source_number][j][k]->GetFitFunction()->SetLineColor(kGreen+2);		///
			calhist[source_number][j]->GetXaxis()->SetRangeUser(new_low-50,new_upp+50 );
			C->Update();
			cout << "ARE YOU HAPPY WITH THIS FIT? (answer y or Y for yes)" << endl;
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
			cout << source_number << "\t" << fitnum << " " << gNonLin[source_number][j]->GetPointX(fitnum) << " " << gNonLin[source_number][j]->GetPointY(fitnum) << endl;			
			gNonLin[source_number][j]->SetPoint(fitnum , gNonLin[source_number][j]->GetPointX(fitnum), cent - gNonLin[source_number][j]->GetPointX(fitnum) );			
			cout << source_number << "\t" << fitnum << " " << gNonLin[source_number][j]->GetPointX(fitnum) << " " << gNonLin[source_number][j]->GetPointY(fitnum) << endl;			
			gNonLin[source_number][j]->Write(gNonLin[source_number][j]->GetName(),TObject::kOverwrite);
		}
	}
}

