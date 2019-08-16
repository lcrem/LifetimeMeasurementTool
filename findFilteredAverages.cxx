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

// string howManyAvg[] = {"justAvg", "avg50", "avg100", "avg200"};
// int howManyAvgInt[] = {1000, 50, 100, 200};
// int howManyGraphs[] = {1, 20, 10, 5};

string howManyAvg[] = {"justAvg","avg100"};
int howManyAvgInt[] = {1000, 100};
int howManyGraphs[] = {1, 10};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 2;
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

  string outfile     = basename + Form("PrM%i", whichPrM) + "_filtAvg.root";

  FILE *felog;
  felog = fopen(Form("%s/Elog.txt", basename.c_str()), "r");

  char elog[10000];
  fgets(elog, 100000, felog);

  string selog = Form("%s", elog);
  fclose(felog);


  string chname[3];
  string chnamenice[3] = {"anode", "cathode", "diff"};
  chname[2] += "diff";

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
  
  cout << "Output file for filtered avg is " << outfile << endl;
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  TGraph *gint[3];
  TGraph *gfil[3];
  TGraph *gfiltmp[3];

  TGraph *power[3];

  bool saveCanvas = false;

  TCanvas *c = new TCanvas("c");

  for (int inum=0; inum<nnum; inum++){

    // if(inum==0) saveCanvas=true;
    // else saveCanvas=false;

    double finalNumbers[2][3]; // [0 anode, 1 cathode] [0 amplitude, 1 start time, 2 peak time]

    out->cd();
    TDirectory *dir;
    if (inum>0){
      dir  = out->mkdir(howManyAvg[inum].c_str());
    }
    
    for (int igraph=0; igraph<howManyGraphs[inum]; igraph++){
      
      

      for (int ich=0; ich<2; ich++){
	string f1 = basename + "/RawAverages_" + chname[ich] +".root";
        
	TFile *file1 = new TFile(f1.c_str(), "read");
	//cout << file1->IsOpen() << endl;
     
	string gname ;
	if (howManyGraphs[inum]==1) gname+= "justAvg";
	else gname += Form("%s/g%s_%d", howManyAvg[inum].c_str(), howManyAvg[inum].c_str(), igraph);

	TGraph *g1 = (TGraph*)file1->Get(gname.c_str());
	//	cout << g1 << " " <<  gname << endl;
	g1->SetName("g1");        

	file1->Close();
        
	out->cd();

	g1 = UsefulFunctions::translateGraph(g1, -(0.05*g1->GetN()*2));

	UsefulFunctions::zeroBaseline(g1);
	
	// convert from nanoseconds to seconds
	for (int i=0; i<g1->GetN(); i++) g1->GetX()[i]*=1.e-9;
	

	  // if (inum==0){
	  //   g1->Write(Form("gorig_%s", chname[ich].c_str()));
	  //   power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(g1);
	  //   power[ich]->Write(Form("gpow_orig_%s", chname[ich].c_str()));
	  // }
	
	// Interpolate
	gint[ich] = FFTtools::getInterpolatedGraph(g1, 32*2.e-9);

      }
      
      gint[2] = new TGraph(gint[0]->GetN());
      
      for (int i=0; i<gint[0]->GetN(); i++) gint[2]->SetPoint(i, gint[0]->GetX()[i], gint[0]->GetY()[i]-gint[1]->GetY()[i]);

      for (int ich=0; ich<3; ich++){
	// Filter out everything above 200kHz
	gfil[ich] = FFTtools::simplePassBandFilter(gint[ich], 0., 200000000.);
	
	// Apply kernel 100
	gfil[ich] = smoothGraph(gfil[ich], 100);

	// gfiltmp[ich] = FFTtools::simplePassBandFilter(gint[ich], 0., 50000000.);
	// gfiltmp[ich] = smoothGraph(gfiltmp[ich], 100);

	// gdiff[ich] = (TGraph*)g1->Clone();	
	

	//	delete g1;

	if (inum==0){
	  // gdiff[ich]->Write(Form("gsmoo_%d", ich));
	  // gint[ich]->Write(Form("gint_%d", ich));
	  gfil[ich]->Write(Form("gfil_%s", chname[ich].c_str()));

	  // gfiltmp[ich]->Write(Form("gfiltmp_%s", chname[ich].c_str()));
	  // power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gdiff[ich]);
	  // power[ich]->Write(Form("gpow_smooth%d", ich));
	  // power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gint[ich]);
	  // power[ich]->Write(Form("gpow_int%d", ich));
	  // power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gfil[ich]);
	  // power[ich]->Write(Form("gpow_fil_%s", chname[ich].c_str()));
	  //   power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gfiltmp[ich]);
	  //   power[ich]->Write(Form("gpow_filtmp%d", ich));

	} else {
	  dir->cd();
	  gfil[ich]->Write(Form("gfil_%s_%d", chname[ich].c_str(), igraph));
	}

      }
      
      cout << "\r " << igraph+1 << " / " << howManyGraphs[inum] << flush;
    }
    
    if (inum==0){
      c->cd();
      TPaveText *pav = new TPaveText(0.25, 0.85, 0.99, 0.95, "NB NDC");
      pav->SetFillColor(kWhite);
      pav->AddText(elog);
      double min = TMath::MinElement(gfil[0]->GetN(), gfil[0]->GetY());
      double max = TMath::MaxElement(gfil[1]->GetN(), gfil[1]->GetY());
      gfil[0]->GetYaxis()->SetRangeUser(min*1.1, max*1.1);
      gfil[0]->SetTitle(Form("PrM%d, Filtered Averages;Time [s];Amplitude [mV]", whichPrM));
      gfil[0]->Draw("Al");
      gfil[1]->SetLineColor(kRed);
      gfil[1]->Draw("l");
      pav->Draw();
      c->Print(Form("%s/PurityMonitor%d_filAvg.png", basename.c_str(), whichPrM));
    }
    
    cout << endl;
  }


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
