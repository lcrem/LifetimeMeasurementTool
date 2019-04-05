#include "RawDigitiser.h"
#include "LifetimeConventions.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;


string findFilename(string dirname, int itdc);

int main(int argc, char *argv[]){

  Int_t run;
  string baseDir;
  
  if((argc!=2)&&(argc!=3)){
    std::cerr << "Usage 1: " << argv[0] << " [irun] (basedir)" << std::endl;
    return 1;
  } else {
    run = atoi(argv[1]);
  if (argc==3) baseDir += argv[2];
  else  baseDir +=  whereIsMyPrMdata;
  }

  string dirname = Form("%s/run%d/", baseDir.c_str(), run);

  // string filename[2];
  // for (int itdc=0; itdc<2; itdc++){
  //   filename[itdc] = findFilename(dirname, itdc+1);
  // }

  // string line;

   // Waveform .dat file to read in 
  string inputFile = "/Users/linda/DUNE/exampleDataFromDigitiser/ST_950V_200ns_007_wf_0_test.dat";
  string outputFile = "output.root";
  
  ifstream digitiserFile;
  digitiserFile.open(inputFile.c_str());
    
  if(digitiserFile.is_open()){
    
    std::cout << "File selected to analyse: " << inputFile << std::endl;
    
    string data;
    
    // Find the full window of the waveform (at what point does the left column go back to 0?)
    // Counter to get number of total lines in file:
    int pointsCounter = 0;

    // Counter to fill trees and arrays
    int counter = 0;

    // Counter to count number of waveforms
    int waveform_counter = 0;

    // If the timebase doesn't start at 0 we want to find out what it does start at:
    Double_t timebase = 0.;

    // Variables to find bin corresponding to the start of the long gate
    int longgate_previous = 0;
    int longgate_start = 0;

    // Open file:       
    while (getline(digitiserFile, data)) {
     
      stringstream dataReadin;
      dataReadin << data;

      RawDigitiser raw;
      
      dataReadin >> raw.timeStamp >> raw.ADC >> raw.baseline >> raw.trigger >> raw.longGate >> raw.shortGate >> raw.zero;
      
      // If the first line starts at zero, then set that to be the start of the timebase. When the next timebase is reached (i.e. the start of another waveform) break out of the loop.
      if (pointsCounter == 0) {
	timebase = raw.timeStamp;
      }
      else if (timebase == raw.timeStamp) {
	break;
      }
      
      // Increment counter to find out how many points per waveform
      ++pointsCounter;
      
    }

    cout << "Number of points per wavelength: " << pointsCounter << endl;
    
    // Define an integer to have the number of points per waveform to use to define arrays
    const unsigned int nLines = pointsCounter;
    
    // Go back to the beginning of the file to read the first waveform again
    digitiserFile.clear();
    digitiserFile.seekg(0, ios::beg);

    // Define graph/hist and its arrays
    double* timebase_value;
    timebase_value = new double[nLines];

    double* amplitudeVolts;
    amplitudeVolts = new double[nLines];

    Int_t longgate_value;
    Int_t longgate_start_bin;

    Int_t countGraphs=0;

    TFile *fout = new TFile(outputFile.c_str(), "recreate");
    
    // Re-read in file variables
    while (getline(digitiserFile, data)) {

      stringstream dataReadin;
      dataReadin << data;
      
      RawDigitiser raw;
      
      dataReadin >> raw.timeStamp >> raw.ADC >> raw.baseline >> raw.trigger >> raw.longGate >> raw.shortGate >> raw.zero;
      
      // Fill arrays and tree
 
      timebase_value[counter] = raw.timeStamp;
      amplitudeVolts[counter] = raw.ADC*1e-3; // convert ADC (mV) to V
      
      longgate_previous = longgate_value;
      
      longgate_value = raw.longGate;

      if (longgate_value > longgate_previous) {
	longgate_start = longgate_value;
	longgate_start_bin = counter + 1;
	//cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	//cout << "Long gate value at change: " << longgate_start << " for bin number " << longgate_start_bin << endl;
      }

      if (counter == pointsCounter - 1) {
	
	TGraph *g = new TGraph(counter, timebase_value, amplitudeVolts);
	g->Write(Form("graph%d",countGraphs));
	counter = 0;
	countGraphs++;
	
      }
      else {++counter;}
    } //close while loop to read in file
    
    cout << "I've made it outside of the loop!!!!!!" << endl;
    cout << "Number of waveforms read in: " << countGraphs << endl;
    fout->Close();

    digitiserFile.close();

  }//close if to see if file is open

  cout << "And now I've made it outside the if file open loop" << endl;
  
    
  return 0;
}


string findFilename(string dirname, int itdc){
  
  DIR *dp;
  dirent *d;
  
  
  if((dp = opendir(dirname.c_str())) == NULL)
    perror("opendir");

  while((d = readdir(dp)) != NULL){
    if(!strcmp(d->d_name,".") || !strcmp(d->d_name,".."))
      continue;
    
    string temp = d->d_name;
    if (temp.find(Form("parsed_%d_", itdc)) != std::string::npos){
      if (temp.find(".root") != std::string::npos) continue;
      else return dirname+"/"+temp;
      
    }
  }  
  return "nothing";
}
