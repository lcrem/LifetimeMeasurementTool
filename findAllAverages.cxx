#include "RawDigitiser.h"
#include "LifetimeConventions.h"
#include "UsefulFunctions.h"

#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "FFTtools.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;

TGraph *graphs[1000];

int fillGraphs(string fname);

int main(int argc, char *argv[]){

  bool recreate;
  string filename="";
  
  if((argc!=2)&&(argc!=3)){
    std::cerr << "Usage 1: " << argv[0] << " [filename] (recreate)" << std::endl;
    return 1;
  } else {
    filename += argv[1];
    if (argc==2){
      recreate=false;
    } else{
      recreate = atoi(argv[2]);
    }
    if (recreate==true) cout << "Recreating averages ... " << endl;
  }

  string finput = filename + ".root";
  string foutput = filename + "_averages.root";

  TGraph *gavg20[100];
  TGraph *gavg50[20];
  TGraph *gavg100[10];
  TGraph *gavg200[10];
  
  
  TFile *fintemp = new TFile(foutput.c_str(), "read");
   
  if (fintemp->IsZombie() || recreate){
     
    int ngraphs = fillGraphs(finput);

    TGraph *justAvg      = UsefulFunctions::justAverage( ngraphs, graphs);
    std::cout << " Done just average " << std::endl;

    TFile *fout = new TFile(foutput.c_str(), "recreate");
    justAvg       ->Write("justAvg");

#ifdef FFTW_UTIL_EXISTS
    double deltat = justAvg->GetX()[1]-justAvg->GetX()[0];
    // Interpolate to 8 times the sampling size
    // Correlate and average
    TGraph *upCorAvg    = FFTtools::interpolateCorrelateAndAverage(deltat,  ngraphs, graphs);
    // Downsample again to initial sampling size
    TGraph *intCorAvg   = FFTtools::getInterpolatedGraph(upCorAvg, deltat);
    std::cout << " Done interpolate, correlate and average " << std::endl;
    intCorAvg->Write("intCorAvg");
#endif
    

    TDirectory *avg20 = fout->mkdir("avg20");
    avg20->cd();
    int num20 = UsefulFunctions::avgSomeGraphs(graphs, 20, gavg20);
    for (int i=0; i<num20; i++){
      avg20->cd();
      gavg20[i]->Write(Form("gavg20_%d", i));
    }

    std::cout << " Done average 20 " << std::endl;

    fout->cd();
    TDirectory *avg50 = fout->mkdir("avg50");
    avg50->cd();
    int num50 = UsefulFunctions::avgSomeGraphs(graphs, 50, gavg50);
    for (int i=0; i<num50; i++){
      avg50->cd();      
      gavg50[i]->Write(Form("gavg50_%d", i));
    }
    std::cout << " Done average 50 " << std::endl;

    TDirectory *avg100 = fout->mkdir("avg100");
    avg100->cd();
    int num100 = UsefulFunctions::avgSomeGraphs(graphs, 100, gavg100);
    for (int i=0; i<num100; i++){
      avg100->cd();      
      gavg100[i]->Write(Form("gavg100_%d", i));
    }
        std::cout << " Done average 100 " << std::endl;

    TDirectory *avg200 = fout->mkdir("avg200");
    avg200->cd();
    int num200 = UsefulFunctions::avgSomeGraphs(graphs, 200, gavg200);
    for (int i=0; i<num200; i++){
      avg200->cd();      
      gavg200[i]->Write(Form("gavg200_%d", i));
    }
    std::cout << " Done average 200 " << std::endl;

    fout->Write();
    delete fout;
  }

  cout << "All averages done and saved in " << foutput << endl;
  return 0;


}



int fillGraphs(string filename){

  TFile *f = new TFile(filename.c_str(), "read");

  double deltat = 0.05E-6;
  int count = 0;
  for (int i=0; i<1000; i++){
    // cout << i << endl;
#ifdef FFTW_UTIL_EXISTS
    TGraph *gtemp = (TGraph*) f->Get(Form("graph%i", i+1));
    graphs[i] = FFTtools::getInterpolatedGraph(gtemp, deltat);
    delete gtemp;
#else
    graphs[i] = (TGraph*) f->Get(Form("graph%i", i+1));
#endif
    if(!graphs[i]) break;
    
    count++;
  }

  f->Close();

  // cout << count << endl;
  return count;
  
}
