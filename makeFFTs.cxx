#include "UsefulFunctions.h"
#include "LifetimeConventions.h"

#include "FFTtools.h"
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
#include "TPaveText.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;

double xmin = 0.05E-3;
double xmax = 0.015E-3;
double ymin = -0.8;
double ymax = +0.8;

void getFields (string s, double fields[3]);


TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 1;
  string basename;
  int whichPrM=0;
  int runNumber;

  if((argc!=4)){
    std::cerr << "Usage : " << argv[0] << " [basename] [runNumber] [whichPrM]" << std::endl;
    return 1;
  } else {
    basename += argv[1];
    runNumber = atoi(argv[2]);
    whichPrM = atoi(argv[3]);
  }

  basename += Form("/Run%03d/", runNumber);

  string outfile     = basename + Form("PrM%i", whichPrM) + "_exampleFFTs.root";

  string chname[2];
  string chnamenice[2] = {"anode", "cathode"};

  if (whichPrM==1){
    chname[0] += "ch1";
    chname[1] += "ch2";
  } else if (whichPrM==2) {
    chname[0] += "ch4";
    chname[1] += "ch5";
  } else {
    cout << " I don't know this purity monitor number " << endl;
    return 1;
  }
  
  cout << "Output file for FFT is " << outfile << endl;
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  TGraph *gdiff[2];
  TGraph *gint[2];
  TGraph *gfil[2];
  TGraph *gfiltmp[2];

  TGraph *power[2];

  string f1 = basename + "/RootifiedFromBinary.root";
  
  TFile *file1 = new TFile(f1.c_str(), "read");
  cout << file1->IsOpen() << endl;
  
  for (int igraph=1; igraph<=nnum; igraph++){
    
    for (int ich=0; ich<2; ich++){

      
      string gname = Form("g_%s_%d", chname[ich].c_str(), igraph);
      
      TGraph *g1 = (TGraph*)file1->Get(gname.c_str());
      cout << g1 << " " <<  gname << endl;
      g1->SetName("g1");        
      
      
      
      out->cd();

      for (int i=0; i<g1->GetN(); i++){
	g1->GetX()[i]*=1.e-9;
	//	g1->GetY()[i]*=1.e-3;
      }
	
      power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(g1);
      power[ich]->SetTitle("Power Spectrum;Frequency [MHz];Amplitude");
      
      power[ich]->Write(Form("%s_fft", gname.c_str()));
    }   
    //	delete g1;
  }
  
  out->Close();
  
  // fclose (pFile);
  
  return 0;
  
}
  




  


void getFields (string s, double fields[3]){
  

  std::string ext("Vcm");
  
  if (s.find(ext) != std::string::npos) {
    
    int place = s.find(ext);
    // if so then strip them off                                                                                                              
    s = s.substr(0, place);
    //      cout << s << endl;                                                                                                                
  }
  
  cout << s << endl;

  string torm = "PrM2 field ";
  
  s = s.substr(torm.size(), s.size());

  cout << s << endl;

  std::replace( s.begin(), s.end(), '.', ' ');
  cout << s << endl;
  
  std::stringstream ss(s);
  ss >> fields[0] >> fields[1] >> fields[2];
  
  cout << " The three fields are " << fields[0] << " " << fields[1] << " " << fields[2] << " V/cm" << endl;
}

TGraph *smoothGraph(TGraph *g, int nnn){

  int n = g->GetN();
  double *x = g->GetX();
  double *y = g->GetY();
  double *newy = new double [1000000];

  int count=0;
  int insidecount=0;
  
  for (int i=nnn; i<(n-nnn); i++){
    newy[count]=0;
    
    insidecount=0;
    for (int j=i-nnn; j<i+nnn; j++){
      newy[count]+=(y[j]);
      insidecount++;
    }
    newy[count]/=(insidecount*1.);
    count++;
  }

  TGraph *gnew = new TGraph(count, x, newy);

  return gnew;

}
