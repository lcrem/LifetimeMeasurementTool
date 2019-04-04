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

  string filename[2];
  for (int itdc=0; itdc<2; itdc++){
    filename[itdc] = findFilename(dirname, itdc+1);
  }
  

  string line;
    
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
