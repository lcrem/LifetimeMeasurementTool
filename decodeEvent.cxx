#include "RawDigitiser.h"
#include "LifetimeConventions.h"
#include "UsefulFunctions.h"

#include "TFile.h"
#include "TGraph.h"
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

double deltat = 2.;      // 2ns delta t
double ADCconv = 0.122;  // ADC conversion
int meanNum = 1000;        // Use the first 50 samples to calculate the baseline

int main(int argc, char *argv[]){

  bool recreate = true;
  string folder="";
  int nmax=-1;

  if((argc!=2)&&(argc!=3)){
    std::cerr << "Usage 1: " << argv[0] << " [folder] (nmax)" << std::endl;
    return 1;
  } else {
    folder += argv[1];
    if (argc==2){
      nmax = 100000;
      //      recreate=false;
    } else{
      //      recreate = atoi(argv[2]);
      nmax = atoi(argv[2]);
    }
    if (recreate==true) cout << "Recreating rootfile ... " << endl;
  }

  string finput = folder + "/Binary.bin";
  string foutput = folder + "/RootifiedFromBinary.root";           

  int ch;

  int MaxNChannels=8;

  uint32_t size;

  FILE *fp;

  char inname[180];

  // Unzip binary file                                                                                                                    
  char unzipfilecmd[1000];
  sprintf(unzipfilecmd, "gzip -d %s.gz", finput.c_str());
  system(unzipfilecmd);
  fp = fopen (finput.c_str(), "rb");

  ch = -1;

  int count = 0;

  TFile *fout = new TFile(foutput.c_str(), "recreate");

  int num[8] = { 0, 0, 0, 0, 0 ,0 ,0 ,0};

  uint16_t DataChannel[2000000];
  while ( !feof(fp) ) {

    fread(&ch, sizeof(ch), 1, fp);
    fread(&size, sizeof(size), 1, fp);

    // cout << __LINE__ << endl;                                                                                                          
    double mean=0;

    //    printf("Channel %i, size %i, number %i  \n", ch, size, num[ch]);
    fread(&DataChannel, sizeof(uint16_t)*size, 1, fp);
    //    cout << DataChannel[0] << " " << DataChannel[100] << " " << DataChannel[2000] << endl;
    // cout << __LINE__ << endl;                                                                                                          
    TGraph *g = new TGraph(size);
    
    for (int i=0; i<size; i++){
      g->GetX()[i] = i*deltat;
      g->GetY()[i] = DataChannel[i]*ADCconv;
      if (i<meanNum) mean += DataChannel[i]*ADCconv;
      //      cout << x[i] << " " << y[i] << endl;                                                                                        
    }

    mean = mean*1.0/meanNum;
    for (int i=0; i<size; i++){
      g->GetY()[i]-=mean;
    }
    //     for (int i=0; i<size; i++) printf("%i ", DataChannel[i] );                                                                     
    //    printf("\n");                                                                                                                  

    num[ch]++;

    g->Write(Form("g_ch%i_%i", ch, num[ch]));

    delete g;
    if (ch==0) cout << "\rDone  " <<  num[ch]*100./1000. << "% "<< flush;

    if (num[ch]>nmax) break;
  }

  fout->Close();

  cout << "Wrote everything in " << folder << endl;



  return 0;


}

