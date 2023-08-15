#include "MyNonLinearities.cc"

/*
WORKFLOW:
- test() will write various histograms, hCorrectedCalHist[MaxSources][NDet], both with and without nonlinearity corrections. 
  matrix, mCorrectedCalMat[MaxSources], has nonlinearities applied. 

PURPOSE:
produce mCorrectedCalMat[MaxSources] with nonlinearity corrections applied and calmat[MaxSources] without nonlinearity corrections applied
for visual comparison

NOTE:
This is a debugging script, not really making anything for further analysis.
*/

TFile *OutputSpectra;

TH1D *hCorrectedCalHist[MaxSources][NDet];
TH2D *mCorrectedCalMat[MaxSources];

int NBins = 10000;
double first = 0.25;
double last = 5000.25;

void GetNonLinearitiesGraph(int DetNum){
/*
REQUIRES:
GraphsFiles .root file with nonlinearities applied.

PURPOSE:
Read in histograms to make matrices with
*/
	int i = DetNum-1;
	gNonLinFull[i] = (TGraphErrors*)GraphsFile->Get(Form("gNonLinFull%d",DetNum) );
	gNonLinFull[i]->Draw("AP");
}

void TestNonLinearities(int source_number,int DetNum){
/*
REQUIRES:
GraphsFiles .root file with nonlinearities applied.

PURPOSE:
Makes histograms with both linearities and nonlinearities applied
*/
	TGraph *gTemp = new TGraph();
	TF1 *fTemp; //
	int i = DetNum-1;
	if( GraphsFile->GetListOfKeys()->Contains(Form("gNonLinFull%d",DetNum) ) ){	
		//fTemp = (TF1*)GraphsFile->Get(Form("gNonLinFull%d",DetNum) );
		GetNonLinearitiesGraph(DetNum);
		hCorrectedCalHist[source_number][i] = new TH1D( Form("hCorrectedCal%s_det%d", Source[source_number].c_str(), DetNum ),
															Form("hCorrectedCal%s_det%d", Source[source_number].c_str(), DetNum ), NBins, first, last);
		for(int k = 0; k < calhist[source_number][i]->GetXaxis()->GetNbins(); k++){
			int counts = calhist[source_number][i]->GetBinContent(k+1);
			double channel =  calhist[source_number][i]->GetBinCenter(k+1);
			double energy;
			for(int l = 0; l < counts; l++){
				double temp_chn = channel+(rand_par.Rndm()-0.5)*calhist[source_number][i]->GetBinWidth(k+1);
				if(temp_chn > gNonLinFull[DetNum-1]->GetPointX(0) && temp_chn < gNonLinFull[DetNum-1]->GetPointX(gNonLinFull[DetNum-1]->GetN()-1) ){					
					energy = temp_chn - gNonLinFull[DetNum-1]->Eval(temp_chn);
					//if( l == 1 )cout << "Correcting energy of event " << temp_chn << " ----> " << energy << endl;
				}
				else{ 
					energy = temp_chn;
				}
				hCorrectedCalHist[source_number][i]->Fill(energy);				
			}		
		}
		OutputSpectra->cd();
		hCorrectedCalHist[source_number][i]->Write(hCorrectedCalHist[source_number][i]->GetName(),TObject::kOverwrite);
		calhist[source_number][i]->Write(calhist[source_number][i]->GetName(),TObject::kOverwrite);
		/*
		//cout << gNonLinFull[DetNum-1]->GetName() << endl;
		for(int x = 0; x < 5000; x++){
			if(x > gNonLinFull[DetNum-1]->GetPointX(0) && x < gNonLinFull[DetNum-1]->GetPointX(gNonLinFull[DetNum-1]->GetN()-1) ){
				gTemp->SetPoint(gTemp->GetN(), x, gNonLinFull[DetNum-1]->Eval(x) );
			}
		}
		gTemp->Draw("P");
		*/
	}
}

void BuildNonLinearCorrectedMatrices(int source_number){
/*
PURPOSE:
Makes matrices from histograms produced by TestNonLinearities() with both linearities and nonlinearities applied
*/	
	int i = source_number;
	mCorrectedCalMat[i] = new TH2D(Form("gEA_%s",Source[i].c_str()), Form("Test calibration with Non Linearity Correction%s",Source[i].c_str()), 64, 0.5, 64.5, NBins, first, last);
	for(int j = 0; j < NDet; j++){
		int DetNum = j+1;
		if( (PACES && DetNum == 49) || (PACES && DetNum == 50) || (PACES && DetNum == 51) || (PACES && DetNum == 52 ) ) continue;
		cout << "Filling detector " << j+1 << endl;
		for(int k = 0; k < hCorrectedCalHist[i][j]->GetXaxis()->GetNbins(); k++){
			mCorrectedCalMat[i]->SetBinContent(j+1 ,k+1, hCorrectedCalHist[i][j]->GetBinContent(k+1) );
		}
		cout << "finished filling detector\n";
	}
	mCorrectedCalMat[i]->Write(mCorrectedCalMat[i]->GetName(),TObject::kOverwrite);
	calmat[i]->Write(calmat[i]->GetName(),TObject::kOverwrite);
}


void MakeOutputSpectra(string filename){

	OutputSpectra = TFile::Open(filename.c_str(),"UPDATE");
}

//SetUp();

void TestNonLinearities(string outputfilename){
/*
REQUIRES:
Calibration sources .root files (analysis trees) . These should have been energy calibrated already.
*/
	//AddSource("152Eu","SumRuns/EnergyCalib_21662_sum.root");
	//AddSource("60Co","SumRuns/EnergyCalib_21601_sum.root");
	//AddSource("133Ba","SumRuns/EnergyCalib_21658_sum.root");
	AddSource("56Co","../SumRuns/EnergyCalib_21598_21600_sum.root");
	GetCalSpectra();
	OpenGraphs("NonLinearGraphs-pietro.root"); // Nonlinearity graphs
	MakeOutputSpectra( outputfilename.c_str() );
	for(int i = 0; i < 64; i++){
		int DetNum = i+1;
		if( (PACES && DetNum == 49) || (PACES && DetNum == 50) || (PACES && DetNum == 51) || (PACES && DetNum == 52 ) ) continue;
		TestNonLinearities(0,DetNum); // Every tag noted as 0 is just an index of the source we are checking.
	}
	BuildNonLinearCorrectedMatrices(0);
	//TestNonLinearities(0,1);

}
