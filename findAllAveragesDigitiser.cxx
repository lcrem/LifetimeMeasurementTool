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

  int ch;
  string filename="";
  
  if((argc!=3)&&(argc!=4)){
    std::cerr << "Usage 1: " << argv[0] << " [filename] [channel] (recreate)" << std::endl;
    return 1;
  } else {
    filename += argv[1];
    ch = atoi(argv[2]);
    if (argc==3){
      recreate=false;
    } else{
      recreate = atoi(argv[3]);
    }
    if (recreate==true) cout << "Recreating averages ... " << endl;
  }

  string finput = filename + ".root";
  string foutput = filename + "_averages" + Form("ch%d", ch) + ".root";

  TGraph *gavg25[40];
  TGraph *gavg50[20];
  TGraph *gavg100[10];
  TGraph *gavg200[5];
  
  
  TFile *fintemp = new TFile(foutput.c_str(), "read");
   
  if (fintemp->IsZombie() || recreate){
    

    TFile *filein = new TFile(finput.c_str(), "read");


    TFile *fout = new TFile(foutput.c_str(), "recreate");

    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    TDirectory *avg25 = fout->mkdir("avg25");
    avg25->cd();
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    int num25 = UsefulFunctions::avgSomeGraphs(filein, 25, gavg25, ch);
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    for (int i=0; i<num25; i++){
      avg25->cd();
      gavg25[i]->Write(Form("gavg25_%d", i));
    }
    fout->cd();
    avg25->Write();

    std::cout << " Done average 25 " << std::endl;
    
    fout->cd();
    TDirectory *avg50 = fout->mkdir("avg50");
    avg50->cd();
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    int num50 = UsefulFunctions::avgSomeGraphs(gavg25, 2, gavg50);
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    for (int i=0; i<num50; i++){
      avg50->cd();      
      gavg50[i]->Write(Form("gavg50_%d", i));

      delete gavg25[i*2];
      delete gavg25[i*2+1];
    }
    std::cout << " Done average 50 " << std::endl;

    TDirectory *avg100 = fout->mkdir("avg100");
    avg100->cd();
    int num100 = UsefulFunctions::avgSomeGraphs(gavg50, 2, gavg100);
    for (int i=0; i<num100; i++){
      avg100->cd();      
      gavg100[i]->Write(Form("gavg100_%d", i));
      delete gavg50[i*2];
      delete gavg50[i*2+1];
    }
    std::cout << " Done average 100 " << std::endl;
    
    // TDirectory *avg200 = fout->mkdir("avg200");
    // avg200->cd();
    // int num200 = UsefulFunctions::avgSomeGraphs(gavg100, 2, gavg200);
    // for (int i=0; i<num200; i++){
    //   avg200->cd();      
    //   gavg200[i]->Write(Form("gavg200_%d", i));
    //   delete gavg100[i*2];
    //   delete gavg100[i*2+1];
    // }
    // std::cout << " Done average 200 " << std::endl;

    TDirectory *avg200 = fout->mkdir("avg200");
    avg200->cd();
    int num200 = UsefulFunctions::avgSomeGraphs(filein, 200, gavg200, ch);
    for (int i=0; i<num200; i++){
      avg200->cd();      
      gavg200[i]->Write(Form("gavg200_%d", i));
      // delete gavg100[i*2];
      // delete gavg100[i*2+1];
    }
    std::cout << " Done average 200 " << std::endl;

    fout->cd();
    avg200->Write();

    cout << "At line " << __LINE__ << endl;
    //TGraph *justAvg      = UsefulFunctions::justAverage(100, gavg100);                                                      
    //std::cout << " Done just average " << std::endl;                                                                                      
    //justAvg->SetTitle(";Time [ns];Amplitude [mV]");  
    //justAvg->Write("justAvg");
    

    fout->cd();
    
    cout << "At line " << __LINE__ << endl;
    TGraph *justAvg      = UsefulFunctions::justAverage(100, gavg100);                                                      
    std::cout << " Done just average " << std::endl;                                                                                      
    justAvg->SetTitle(";Time [ns];Amplitude [mV]");  
    justAvg->Write("justAvg");


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
  for (int i=1; i<=1000; i++){
    // cout << i << endl;
// #ifdef FFTW_UTIL_EXISTS
//     TGraph *gtemp = (TGraph*) f->Get(Form("g_ch1_%i", i+1));
//     graphs[i] = FFTtools::getInterpolatedGraph(gtemp, deltat);
//     delete gtemp;
// #else
    graphs[i] = (TGraph*) f->Get(Form("g_ch4_%i", i+1));
// #endif
    if(!graphs[i]) break;
    
    count++;
  }

  f->Close();

  // cout << count << endl;
  return count;
  
}
