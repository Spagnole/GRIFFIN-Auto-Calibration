#define EnergyCalSelector_cxx
// The class definition in EnergyCalSelector.h has been generated automatically
#include "EnergyCalSelector.h"

int kValue_grif = 379; // default k-value for not piled-up hits
int kValue_zds = 59; // default k-value for not piled-up hits
int kValue_paces = 200; // default k-value for not piled-up hits //paces has different values for different detectors it seems 
//flags
bool PileUpRej	 = true; // If true, pile-up is rejected
bool AddbOn		 = false; // If true, addback spectra are sorted
bool SuppOn		 = true; // If true, suppressed spectra are sorted
bool AddbSuppOn	 = false; // If true, addback-suppressed spectra are sorted
// Ancillaries
bool BgoOn		 = true; // If true, spectra and matrices for BGos are sorted
bool ZdsOn		 = false; // If true, spectra and matrices for Zds are sorted
bool SceptarOn	 = false; // If true, spectra and matrices for Sceptar are sorted (ONLY IF ZDS IS ALSO ON!)
bool LabrOn		 = false; // If true, spectra and matrices for LaBr are sorted
bool LabrSuppOn	 = false; // If true, spectra and matrices for LaBr with suppression are sorted
bool PacesOn	 = false; // If true, spectra and matrices for Paces are sorted


////Global variables for the time coincidence windows

double ggtime   = 250.;//Time coincidence windows for GRIFFIN-GRIFFIN, assumed symmetrical
double gztime_l = -200.;//Time coincidence windows for GRIFFIN-ZDS, lower limit
double gztime_u = 0.;//Time coincidence windows for GRIFFIN-ZDS, upper limit
double gptime_l = -150.;//Time coincidence windows for GRIFFIN-PACES, lower limit
double gptime_u = -50.;//Time coincidence windows for GRIFFIN-PACES, upper limit

double randtime_l=1000.;//Random time coincidence window, common for all coincidences, lower limit
double randtime_u=2000.;////Random time coincidence window, common for all coincidences, upper limit

double gz_randtime_l=800.;//Random time coincidence window, common for all coincidences, lower limit
double gz_randtime_u=1800.;////Random time coincidence window, common for all coincidences, upper limit

double gp_randtime_l=200.;//Random time coincidence window, common for all coincidences, lower limit
double gp_randtime_u=1200.;////Random time coincidence window, common for all coincidences, upper limit

double MyCycleLength = 45+13+1.5+0.5;

// prompt coincidences
bool PromptCoincidence(TGriffinHit *h1, TGriffinHit *h2){  //Check if GRIFFIN hits are less then 120 ns apart.
   return std::fabs(h1->GetTime() - h2->GetTime()) < ggtime;
}
bool PromptCoincidence(TGriffinHit *h1, TZeroDegreeHit *h2){  // coincidence window GRIFFIN-ZDS: (-450,-325) ns
   return gztime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gztime_u;
}
bool PromptCoincidence(TGriffinHit *h1, TPacesHit *h2){  // coincidence window GRIFFIN-PACES: (-450,-325) ns
   return gptime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gptime_u;
}

// time random
bool TimeRandom(TGriffinHit *h1, TGriffinHit *h2){  //Check if hits are less then 1000 ns apart.
   return randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TZeroDegreeHit *h2){  //Check if hits are less then 1000 ns apart.
   return gz_randtime_l < std::fabs(h1->GetTime() - h2->GetTime()) && std::fabs(h1->GetTime() - h2->GetTime()) < gz_randtime_u;
}
bool TimeRandom(TGriffinHit *h1, TPacesHit *h2){  //Check if hits are less then 1000 ns apart.
   return gp_randtime_l < h1->GetTime() - h2->GetTime() && h1->GetTime() - h2->GetTime() < gp_randtime_u;
}


void EnergyCalSelector::CreateHistograms() {
   // get cycle length from ODB - this is in nanoseconds! (Actually it's not in ns, I got it in mu s!)
   //fCycleLength = fPpg->OdbCycleLength();
   //std::cout<<"got ODB cycle length "<<fCycleLength<<" mu s, "<<fCycleLength/1e6<<" s"<<std::endl;
	
////////////////////////////////////////////////////////////////////////////////////////////////////
	fCycleLength = fPpg->OdbCycleLength();
	std::cout<<"got ODB cycle length "<<fCycleLength<<" ns, "<<fCycleLength/1e6<<" s"<<std::endl;
	
	
		
	//energy spectra
	fH2["gEA"] = new TH2D("gEA","#gamma singles vs. crystal #", 64, 0.5, 64.5, 10000, 0.25, 5000.25); //griffin energy (all channels)
	fH2["gCA"] = new TH2D("gCA","#gamma charge vs. crystal #", 64, 0.5, 64.5, 10000, 0.25, 5000.25); //griffin energy (all channels)
	fH2["gEA_NoNonLin"] = new TH2D("gEA_NoNonLin","#gamma singles (Ignore NonLinearity correction) vs. crystal #", 64, 0.5, 64.5, 10000, 0.25, 5000.25); //griffin energy (all channels)
	

/////////////////////////////////////////////////////////////////////////// Multiplicity Histograms
///////////////////////////////////////////////////////////////////////////  GRIFFIN ENERGY Spectra
   
   //Send histograms to Output list to be added and written.
   for(auto it : fH1) {
      GetOutputList()->Add(it.second);
	  // it->Sumw2();//If this works, it will propagate errors correctly
   }
   for(auto it : fH2) {
      GetOutputList()->Add(it.second);
	  // it.Sumw2();//If this works, it will propagate errors correctly
   }
   for(auto it : fHSparse) {
      GetOutputList()->Add(it.second);
	  // it.Sumw2();//If this works, it will propagate errors correctly
   }
}


void EnergyCalSelector::FillHistograms() {

///////////////////////////////////////   //GRIFFIN No addback with suppression

	for(auto g = 0; g < fGrif->GetSuppressedMultiplicity(fBgo); ++g) {
		auto grif = fGrif->GetSuppressedHit(g);
		if(PileUpRej) if(grif->GetKValue() != kValue_grif) continue;
		fH2.at("gEA")->Fill(grif->GetArrayNumber(), grif->GetEnergy());
		fH2.at("gCA")->Fill(grif->GetArrayNumber(), grif->GetCharge());	
		fH2.at("gEA_NoNonLin")->Fill(grif->GetArrayNumber(), grif->GetChannel()->CalibrateENG(grif->GetCharge()));	
	}
 }



