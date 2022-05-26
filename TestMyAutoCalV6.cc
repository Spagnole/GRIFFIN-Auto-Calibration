#include "MyAutoCalV6.cc"

int NBins = 10000;
double first = 0.25;
double last = 5000.25;

void CalSpectra(int source_number,int DetNum){	

	TF1 *fTemp; //
	int i = source_number;
	if( GraphsFile->GetListOfKeys()->Contains(Form("Fit_Det%d", DetNum) ) ){
		fTemp = (TF1*)GraphsFile->Get(Form("Fit_Det%d", DetNum));
		double a0 = fTemp->GetParameter(0);
		double a1 = fTemp->GetParameter(1);		
		int j = DetNum-1;		
		calhist[i][j] = new TH1D( Form("hCal%s_det%d", Source[i].c_str(), DetNum ),Form("hCal%s_det%d", Source[i].c_str(), DetNum ), NBins, first, last);
		cout << hist[i][j]->GetName() << endl;
		cout << a0 << "\t" << a1 << endl;
		for(int k = 0; k < hist[i][j]->GetXaxis()->GetNbins(); k++){
			int counts = hist[i][j]->GetBinContent(k+1);
			double channel =  hist[i][j]->GetBinCenter(k+1);
			for(int l = 0; l < counts; l++){
				double temp_chn = channel+(rand_par.Rndm()-0.5)*hist[i][j]->GetBinWidth(k+1);
				double energy = a0 + a1*temp_chn;
				calhist[i][j]->Fill(energy);				
			}		
		}
		calhist[i][j]->Write();
	}
	else if( GraphsFile->GetListOfKeys()->Contains(Form("Fit_Det%d", DetNum) ) ){
		cout << "No calibration exists for detector " << DetNum << endl;
		cout << "Must write energy calibration using WriteGraph() function\n";	
	}
}




void BuildCalMatrices(int source_number){
	
	int i = source_number;
	calmat[i] = new TH2D(Form("gEA_%s",Source[i].c_str()), Form("Test calibration %s",Source[i].c_str()), 64, 0.5, 64.5, NBins, first, last);
	for(int j = 0; j < NDet; j++){
		if( j == 48 || j == 49 || j == 50 || j == 51 ) continue;
		for(int k = 0; k < calhist[i][j]->GetXaxis()->GetNbins(); k++){
			calmat[i]->SetBinContent(j+1 ,k+1, calhist[i][j]->GetBinContent(k+1) );
		}
	}
	calmat[i]->Write();
}


void TestCalibration(string outputfile, string graphfilename){
	//MyAutoCal();
	
	AddSource("56Co","EnergyCalib_21600_008.root");
	GetSpectra();
	OpenGraphs( graphfilename.c_str() );
	//GetCalSpectra();
	TFile *fTest = TFile::Open(outputfile.c_str(),"RECREATE");
	for(int i = 0; i < 64; i++) CalSpectra(0,i+1);
	BuildCalMatrices(0);
}
