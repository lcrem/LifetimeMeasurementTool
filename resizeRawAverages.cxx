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


// string howManyAvg[5] = {"justAvg", "avg25", "avg50", "avg100", "avg200"};
// int howManyAvgInt[5] = {1000, 25, 50, 100, 200};
// int howManyGraphs[5] = {1, 40, 20, 10, 5};

string howManyAvg[] = {"justAvg", "avg50", "avg100"};
int howManyAvgInt[] = {1000, 50, 100};
int howManyGraphs[] = {1, 20, 10};

// string howManyAvg[] = {"justAvg","avg100"};
// int howManyAvgInt[] = {1000, 100};
// int howManyGraphs[] = {1, 10};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 3;
  string basename;
  int whichPrM=0;
  int runNumber;
  int channel;


  if((argc!=4)){
    std::cerr << "Usage : " << argv[0] << " [basename] [runNumber] [channel]" << std::endl;
    return 1;
  } else {
    basename += argv[1];
    runNumber = atoi(argv[2]);
    channel = atoi(argv[3]);
  }

  basename += Form("/Run%03d/", runNumber);

  string outfile     = basename + Form("RawAverages_ch%i", channel) + "_lil.root";

  cout << "Output file is " << outfile << endl;
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  TGraph *gint[3];

  string f1 = basename + "/RawAverages_" + Form("ch%d", channel) + ".root";
  
  TFile *file1 = new TFile(f1.c_str(), "read");
  //cout << file1->IsOpen() << endl;
  
  for (int inum=0; inum<nnum; inum++){

    out->cd();
    TDirectory *dir;
    if (inum>0){
      dir  = out->mkdir(howManyAvg[inum].c_str());
    }
    
    for (int igraph=0; igraph<howManyGraphs[inum]; igraph++){
      
      string gname ;
      if (howManyGraphs[inum]==1){
	gname+= "justAvg";
	if (!file1->GetListOfKeys()->Contains(gname.c_str())) continue;
      } else {
	gname += Form("%s/g%s_%d", howManyAvg[inum].c_str(), howManyAvg[inum].c_str(), igraph);	
      }

      TGraph *g1 = (TGraph*)file1->Get(gname.c_str());
      //      cout << g1 << " " <<  gname << endl;
      g1->SetName("g1");        
      
     
      out->cd();

      if (inum==0){
	g1->Write("justAvg");
      } else {
	dir->cd();
	g1->Write(Form("g%s_%d", howManyAvg[inum].c_str(), igraph));
      }
      
      
      
      cout << "\r " << igraph+1 << " / " << howManyGraphs[inum] << flush;

      delete g1;
    }
    
    cout << endl;
  }
  
  file1->Close();  
  out->Close();
  
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

  nnn/=2;

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
