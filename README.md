# GRIFFIN-Auto-Calibration
Tool for calibrating GRIFFIN detectors and writing output calibration file
Please note you will require an initial energy calibration file to run this tool
The online calibration from your experiment should be suitable for this

This tool has a few steps
1. Run the EnergyCalSelector.C selector code
2. Run the MyAutoCalV6.cc code
3. Load the output cal file from step 2 to your analysis file
4. Run the MyNonLinearities.cc code
5. Load the nonlinearities cal file to your analysis file
6. Run the EnergyCalSelector.C selector code to see the qualitiy of your now calibrated spectra


################
The fisrt step is to run the EnergyCalSelector.C selector on the analysis files to produce spectra
This will create three matrices

gCA - This is the uncalibrated energy GetCharge() versus detector array number, Used by the MyAutoCalV6.cc for the energy calibration
gEA - This is the calibrated energy GetEnergy() versus detector array number, final output
gEA_NoNonLin - This is a calibrated energy spectrum without the nonlinearities applied, this matrix is needed for the MyNonLinearities.cc code
#################

The next step is to use the MyAutoCalV6.cc code. This code has a lot of functions and a number of step
1. AddSource("60Co","SumRuns/EnergyCalib_21601_sum.root");
This step can be run several times to add the different calibration sources
2. ReadLitEnergies();
This will read in your literature values
The files should be named as LitEn_XXX.dat where XXX is your source name e.g. 60Co
3. GetSpectra()
This will get the spectra from the root file in step loaded in step 1
4. GetCalPars()
This will extract the initial calibration parameters from the analysis tree
5. OpenFits(string FileName)
This will open a root file to write to. It will save the individaul spectra and the fits for a quality check
You can open a file that already exists. This allows you to come back to your earlier calibration so you don't have to do it all in one go.
6. OpenGraphs(string FileName)
This will open a root file for your graphs for the calibration, your peak fits are added to graphs to then get a linear calibration 
You can open a file that already exists. This allows you to come back to your earlier calibration so you don't have to do it all in one go.
7. MyFits(Det Number)
This will fit all the peaks in your list (given in step 2) for the given detector number
8. FitCalibration(Det Number)
This will perform a calibration with the data obtained from step 7.
9. WriteFit(Det Number)
If you are happy with the calibration from step 8, write the graph to the root file.
10. BuildEnergyCalFile("yourfile.cal")
If you are happy with your calibration, you can now write the calibration file. 

If you are not happy with some of your calibrations you can use the ReFitPeak() function.
This will allow you to refit a given peak for a given detector
run the function for the detector you want to fit a new peak for and follow the instructions printed to the terminal
You can also erase the peak from the graph if you cannot get a good fit
You will have to run FitCalibration() [step 8] again if you refit peaks

You can also add a peak to your graph that is not already in your list of gammas with the AddPeak(int DetNum, int source_number, double NewLitEnergy) function

###############
It is now recommended to test the quality of your calibration
This can be performed with the TestMyAutoCalV6.cc code
Load the code in grsisort and run the function
TestCalibration(string outputfile, string graphfilename)
Give the name of your root file containing your calibration graphs and give a name to write the spectra to a file
See the TestCalibration() function for more info

###############
Load your cal file to the analysis files and run the EnergyCalSelector.C selector again
###############
Now it is time to perform the non linearity correction
1. open grsisort and load the code
2. AddSource("56Co","../SumRuns/EnergyCalib_21598_21600_sum.root");
3. ReadLitEnergies()
4. GetCalSpectra()
5. GetCalPars()
6. OpenFits("CalSpectra.root");
7. OpenGraphs("NonLinearGraphs.root");
8. FitNonLinearities(Det Number)
9. SortNonLinearities(Det Number)
10. MakeNonLinCalFile("nonlin.cal")

As with the calibration file, you can re-fit peaks that you are not happy with using the ReFitNonLinearities() function

You can test the quaility of your non-linearity corrections using the TestMyNonLinearities.cc code
Load the code in grsisort and use the TestNonLinearities(string outputfilename) function
###################
If you are happy with your non-linearity corrections
save the non-linearity cal file to your analysis files
Run the EnergyCalSelector.C selector one final time to see the quality of your calibrated spectra

####################
You are finished


