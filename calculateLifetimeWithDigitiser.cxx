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


string howManyAvg[5] = {"justAvg", "avg25", "avg50", "avg100", "avg200"};
int howManyAvgInt[5] = {1000, 25, 50, 100, 200};
int howManyGraphs[5] = {1, 40, 20, 10, 5};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 5;
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

  string outfile2    = basename + Form("PrM%i", whichPrM) + "_lifeInfo.root";
  string outfile     = basename + Form("PrM%i", whichPrM) + "_filtAvg.root";
  string outtxtfile  = basename + Form("PrM%i", whichPrM) + "_lifetime_withErrors.txt";
 

  FILE *felog;
  felog = fopen(Form("%s/Elog.txt", basename.c_str()), "r");
  
  char elog[10000];
  fgets(elog, 100000, felog);
  
  string selog = Form("%s", elog);
  fclose(felog);


  double fields[3], distance[3], tTheory[3];
  getFields (elog, fields);

  string chname[2];
  string chnamenice[2] = {"anode", "cathode"};

  if (whichPrM==1){
    distance[0] = PrM1distance[0];
    distance[1] = PrM1distance[1];
    distance[2] = PrM1distance[2];
    chname[0] += "ch1";
    chname[1] += "ch2";
  } else if (whichPrM==2) {
    distance[0] = PrM2distance[0];
    distance[1] = PrM2distance[1];
    distance[2] = PrM2distance[2];
    chname[0] += "ch4";
    chname[1] += "ch5";
  } else {
    cout << " I don't know this purity monitor number " << endl;
    return 1;
  }
  
  cout << "Output file for filtered avg is " << outfile << endl;
  cout << "Output file for tree of lifetime info is " << outfile2 << endl;
  
  // double timedelay=0;
  
  // string stimedelay = basename + fieldname  + ".ch1.traces_averages.root";
  // TFile *ftimedelay = new TFile (stimedelay.c_str(), "read");
  // TGraph *gtimedelay = (TGraph*)ftimedelay->Get("justAvg");
  // double *xTimeDelay = gtimedelay->GetX();
  // double *yTimeDelay = gtimedelay->GetY();
  // double baseY = yTimeDelay[0];
  // double step = TMath::MaxElement(gtimedelay->GetN(), yTimeDelay)*0.9;
  // double timeStep = xTimeDelay[1]-xTimeDelay[0];
  // for (int ip=0; ip<gtimedelay->GetN(); ip++){
  //   if (yTimeDelay[ip]>(step)){
  //     timedelay = xTimeDelay[ip];
  //     break;
  //   }
  // }
  // ftimedelay->Close();
  // cout << "The time delay is " << timedelay << endl;

  double timeStep = 2.E-9; // 2ns
  for (int i=0; i<3; i++){
    tTheory[i] = distance[i]/UsefulFunctions::ICARUSpolynomial(fields[i]);
    // cout << distance[i] << " " << fields[i] << " " << tTheory[i] << " " << timeStep << " " << UsefulFunctions::getSmoothingNumber(timeStep, tTheory[i]) <<endl;
  }
  int smoothing[2];
  smoothing[0] = 100;// UsefulFunctions::getSmoothingNumber(timeStep, tTheory[2]);
  smoothing[1] = 100;// UsefulFunctions::getSmoothingNumber(timeStep, tTheory[0]);
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  TGraph *gdiff[2];
  TGraph *gint[2];
  TGraph *gfil[2];

  TGraph *power[2];

  bool saveCanvas = false;

  
  TCanvas *c = new TCanvas("c");

  double lifeApprox, lifetime;
  double QK, QA, QKcorr, QAcorr;
  double t1, t2, t3;
  int numAveraged;
  

  TFile *outLife = new TFile(outfile2.c_str(), "recreate");
  TTree *lifeTree = new TTree("lifeTree", Form("Lifetime info for Purity Monitor %d", whichPrM));
  lifeTree->Branch("run",          &runNumber    );
  lifeTree->Branch("numAveraged",  &numAveraged  );
  lifeTree->Branch("lifeApprox",   &lifeApprox   );
  lifeTree->Branch("lifetime",     &lifetime     );
  lifeTree->Branch("QK",           &QK           );
  lifeTree->Branch("QA",           &QA           );
  lifeTree->Branch("QKcorr",       &QKcorr       );
  lifeTree->Branch("QAcorr",       &QAcorr       );
  lifeTree->Branch("t1",           &t1           );
  lifeTree->Branch("t2",           &t2           );
  lifeTree->Branch("t3",           &t3           );


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

	string gname ;
	if (howManyGraphs[inum]==1) gname+= "justAvg";
	else gname += Form("%s/g%s_%d", howManyAvg[inum].c_str(), howManyAvg[inum].c_str(), igraph);

	TGraph *g1 = (TGraph*)file1->Get(gname.c_str());
	//	cout << f1 << " " <<  gname << endl;
	g1->SetName("g1");
        
	file1->Close();
        
	out->cd();

	g1 = UsefulFunctions::translateGraph(g1, -(0.05*g1->GetN()*2));

	UsefulFunctions::zeroBaseline(g1);
	
	// convert from nanoseconds to seconds
	for (int i=0; i<g1->GetN(); i++) g1->GetX()[i]*=1.e-9;
	
	// g1->Write(Form("gorig_%d", ich));
	// power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(g1);
	// power[ich]->Write(Form("gpow_orig%d", ich));

	// Interpolate
	gint[ich] = FFTtools::getInterpolatedGraph(g1, 32*2.e-9);

	// Filter out everything above 1 GHz
	gfil[ich] = FFTtools::simplePassBandFilter(gint[ich], 0., 1000000000.);

	// gdiff[ich] = (TGraph*)g1->Clone();	
	// gdiff[ich] = smoothGraph(gdiff[ich], smoothing[ich]);

	if (inum==0){
	  // gdiff[ich]->Write(Form("gsmoo_%d", ich));
	  // gint[ich]->Write(Form("gint_%d", ich));
	  gfil[ich]->Write(Form("gfil_%d", ich));
	} else {
	  dir->cd();
	  gfil[ich]->Write(Form("%s/gfil_ch%d_%d", howManyAvg[inum].c_str(), ich, inum));
	}
	
	// power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gdiff[ich]);
	// power[ich]->Write(Form("gpow_smooth%d", ich));

	// power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gint[ich]);
	// power[ich]->Write(Form("gpow_int%d", ich));

	// power[ich] = FFTtools::makePowerSpectrumVoltsSeconds(gfil[ich]);
	// power[ich]->Write(Form("gpow_fil%d", ich));

	//	delete g1;
      }

      double tlifetime[10];
      
      int ok = UsefulFunctions::calculateLifetime(gfil[0], gfil[1],  whichPrM-1, tTheory, tlifetime, saveCanvas);
      
      numAveraged=howManyAvgInt[inum];
      
      lifeApprox = tlifetime[0];
      lifetime   = tlifetime[1];
      QK         = tlifetime[2];
      QA         = tlifetime[3];
      QKcorr     = tlifetime[4];
      QAcorr     = tlifetime[5];
      t1         = tlifetime[6];
      t2         = tlifetime[7];
      t3         = tlifetime[8];

      outLife->cd();
      lifeTree->Fill();

      cout << "\r " << igraph << " / " << howManyGraphs[inum] << flush;

    }
    if (inum==0){
      c->cd();
      TPaveText *pav = new TPaveText(0.25, 0.85, 0.99, 0.95, "NB NDC");
      pav->SetFillColor(kWhite);
      pav->AddText(elog);
      double min = TMath::MinElement(gfil[0]->GetN(), gfil[0]->GetY());
      double max = TMath::MaxElement(gfil[1]->GetN(), gfil[1]->GetY());
      gfil[0]->GetYaxis()->SetRangeUser(min*1.1, max*1.1);
      gfil[0]->SetTitle(Form("PrM%d, Filtered Averages;Time [ns];Amplitude [mV]", whichPrM));
      gfil[0]->Draw("Al");
      gfil[1]->Draw("l");
      pav->Draw();
      c->Print(Form("%s/PurityMonitor%d_filAvg.png", basename.c_str(), whichPrM));
    }
    
    cout << endl;
  }

  outLife->cd();
  lifeTree->Write();
  outLife->Close();


  
  // FILE * pFile;
  
  // cout << "Writing these info to " << outtxtfile << endl;
  
  
  // pFile = fopen (outtxtfile.c_str(),"w");
  
  // fprintf(pFile, "Is this the problem? \n");  

  // cout << "At line " << __LINE__ << endl;

  // for (int inum=1; inum<nnum; inum++){

  // cout << "At line " << __LINE__ << endl;

  //   double avgLifetime = hpurity[inum]->GetMean();
  //   double errLifetime = hpurity[inum]->GetMeanError();
  //   double rmsLifetime = hpurity[inum]->GetRMS();
  //   double rmserrLifetime = hpurity[inum]->GetRMSError();
  

  //   cout << "THIS IS OUR PURITY AND ERROR: " << avgLifetime << " " << errLifetime << endl;


  //   // TF1 *fgaus = new TF1("fgaus", "gaus");
  //   // fgaus->SetParameter(1, avgLifetime);
  //   // hpurity[inum]->Fit("fgaus");
  //   // cout << "SOMENUMBERS " << fgaus->GetParameter(0) << " " << fgaus->GetParameter(1) << " " << fgaus->GetParameter(2) << endl;

    
  //   fprintf(pFile, "%s \n", howManyAvg[inum].c_str());
  //   fprintf(pFile, "%8.3e %8.3e %8.3e %8.3e \n", avgLifetime , errLifetime, rmsLifetime, rmserrLifetime);
  //   // fprintf(pFile, "tK     : %12.4e \n",  tK  );
  //   // fprintf(pFile, "tGK    : %12.4e \n",  tGK );
  //   // fprintf(pFile, "tGA    : %12.4e \n",  tGA );
  //   // fprintf(pFile, "tA     : %12.4e \n",  tA  );
  //   // fprintf(pFile, "t1     : %12.4e \n",  t1  );
  //   // fprintf(pFile, "t2     : %12.4e \n",  t2  );
  //   // fprintf(pFile, "t3     : %12.4e \n",  t3  );
  //   // fprintf(pFile, "QA     : %12.4e \n",  QA  );
  //   // fprintf(pFile, "QK     : %12.4e \n",  QK  );
  //   // fprintf(pFile, "QA corr: %12.4e \n",  newQA  );
  //   // fprintf(pFile, "QK corr: %12.4e \n",  newQK  );
  //   // fprintf(pFile, "R      : %12.4e \n",  R   );
  //   // fprintf(pFile, "lifetime : %12.4e \n",  lifetime  );
  //   // fprintf(pFile, "lifetime2: %12.4e \n",  lifetime2 );
  
  // cout << "At line " << __LINE__ << endl;
  //   out->cd();
  //   hpurity[inum]->Write(Form("hpurity_%d", inum));
  // }
  // cout << "At line " << __LINE__ << endl;

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
