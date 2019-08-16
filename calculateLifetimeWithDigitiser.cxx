#include "UsefulFunctions.h"
#include "LifetimeConventions.h"

#include "FFTtools.h"
#include "TFile.h"
#include "TF1.h"
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
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>

using namespace std;

double xmin = 0.05E-3;
double xmax = 0.015E-3;
double ymin = -0.8;
double ymax = +0.8;

void getFields (string s, double fields[3]);

double getRawVoltageRMS(string basename, int irun, string chname);

double getLifetimeFromCathode(TGraph *g1, double parlife[10]);

Double_t greenFunction(Double_t *x, Double_t *par);

// string howManyAvg[5] = {"justAvg", "avg25", "avg50", "avg100", "avg200"};
// int howManyAvgInt[5] = {1000, 25, 50, 100, 200};
// int howManyGraphs[5] = {1, 40, 20, 10, 5};

// string howManyAvg[] = {"justAvg", "avg50", "avg100", "avg200"};
// int howManyAvgInt[] = {1000, 50, 100, 200};
// int howManyGraphs[] = {1, 20, 10, 5};

string howManyAvg[] = {"justAvg", "avg100"};
int howManyAvgInt[] = {1000, 100};
int howManyGraphs[] = {1, 10};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 2;
  string basename;
  int whichPrM=0;
  int runNumber;
  int runNoise=-1;

  if((argc!=4 && argc!=5)){
    std::cerr << "Usage : " << argv[0] << " [basename] [runNumber] [whichPrM] (noiseOnlyRunNumber)" << std::endl;
    return 1;
  } else {
    basename += argv[1];
    runNumber = atoi(argv[2]);
    whichPrM = atoi(argv[3]);
    if (argc==5) runNoise = atoi(argv[4]);
  }


  string runname = basename + Form("/Run%03d/", runNumber);

  string infile      = runname + Form("PrM%i", whichPrM) + "_filtAvg.root";
  string outfile     = runname + Form("PrM%i", whichPrM) + "_lifeInfo.root";
  string outtxtfile  = runname + Form("PrM%i", whichPrM) + "_avgLifeInfo.txt";
  string lampFile    = runname + "/RawAverages_ch0.root";



  FILE * outFile;                                                                                                                                                                 
  cout << "Writing some info to " << outtxtfile << endl;                                                                                                                       
  outFile = fopen (outtxtfile.c_str(),"w");   

  FILE *felog;
  felog = fopen(Form("%s/Elog.txt", runname.c_str()), "r");
  
  char elog[10000];
  fgets(elog, 100000, felog);
  
  string selog = Form("%s", elog);
  fclose(felog);


  struct stat t_stat;
  stat(Form("%s", runname.c_str()), &t_stat);
  UInt_t timestamp = t_stat.st_ctime;

  double fields[3], distance[3], tTheory[3];
  getFields (elog, fields);
  cout << " The three fields are " << fields[0] << " " << fields[1] << " " << fields[2] << " V/cm" << endl;
  if (fields[0]==0){
    cout << "I do not usually process fields at 0Vcm, are you sure that is what you want to do?" << endl;
    cout << "... Well sorry I am going to have to quit this program now. Bye!!" << endl;
    return -1;
  }



  if (runNoise==-1){
    cout << "Looking for the 0Vcm run number closest in time " << endl;
    Int_t minDiff=999999999;
    Int_t tmpDiff;
    UInt_t timestampt;
    double fieldst[3];
    

    for  (int irun=runNumber-5; irun<runNumber+5; irun++){

      //      cout << Form("%s/Run%03d/Elog.txt", basename.c_str(), irun) << endl;
      FILE *felogt;

      if (( felogt = fopen(Form("%s/Run%03d/Elog.txt", basename.c_str(), irun), "r"))) {
	
	char elogt[10000];
	fgets(elogt, 100000, felogt);
	
	string selogt = Form("%s", elogt);
	fclose(felogt);
	struct stat t_statt;
	stat(Form("%s/Run%03d/RawAverages_ch0.root", basename.c_str(), irun), &t_statt);
	timestampt = t_statt.st_ctime;
	
	
	string noiseFilet   = basename + Form("/Run%03d/", irun) + Form("PrM%i", whichPrM) + "_filtAvg.root";

	//	cout << noiseFilet << endl;

	struct stat buffer;
	if (stat(noiseFilet.c_str(), &buffer)!=0) continue;
	//	cout << "Here or not?" << endl;


	getFields (elogt, fieldst);
	//	cout << fieldst[0] << " " << fieldst[1] << " " << fieldst[2] << endl;

	if (fieldst[0]>0) continue;

	tmpDiff = int(timestampt) - int(timestamp);
	tmpDiff = TMath::Abs(tmpDiff);
	cout << tmpDiff << " " << minDiff << endl;
	if (tmpDiff < minDiff ) {
	  minDiff=tmpDiff;
	  runNoise=irun;
	  cout << "Best run noise so far " << runNoise << " " << fieldst[0] << " " << minDiff  << endl;
	}
	
	
      }

    }

  }

  cout << "Using run noise " << runNoise << endl;

  string noiseFile   = basename + Form("/Run%03d/", runNoise) + Form("PrM%i", whichPrM) + "_filtAvg.root";


  string chname[2];
  string chnamenice[2] = {"anode", "cathode"};

  TCanvas *c = new TCanvas("c");


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
  
  cout << "Input file for filtered avg is " << infile << endl;
  cout << "Output file for tree of lifetime info is " << outfile << endl;
  

  TFile *flamp = new TFile(lampFile.c_str(), "read");
  TGraph *glamp = (TGraph*)flamp->Get("justAvg");
  double lampIntegral = glamp->Integral();
  double lampMax = TMath::Abs(TMath::MaxElement(glamp->GetN(), glamp->GetY()));
  double lampMin = TMath::Abs(TMath::MinElement(glamp->GetN(), glamp->GetY()));
  double lampPeak;
  if (lampMax>lampMin) lampPeak = lampMax;
  else lampPeak = lampMin;

  flamp->Close();

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

  cout << "T theory in us : " ;
  for (int i=0; i<3; i++){
    tTheory[i] = distance[i]/UsefulFunctions::ICARUSpolynomial(fields[i]);
    cout << tTheory[i]*1e6 <<  " ";
    // cout << distance[i] << " " << fields[i] << " " << tTheory[i] << " " << timeStep << " " << UsefulFunctions::getSmoothingNumber(timeStep, tTheory[i]) <<endl;
  }
  cout << endl;
  int smoothing[2];
  smoothing[0] = 100;// UsefulFunctions::getSmoothingNumber(timeStep, tTheory[2]);
  smoothing[1] = 100;// UsefulFunctions::getSmoothingNumber(timeStep, tTheory[0]);
  
  TFile *fin = new TFile(infile.c_str(), "read");
  if (!fin) {
    cout << "Input file " << infile.c_str() << " does not exist. Please run findFilteredAverages and try again " << endl;
    return -1;
  }

  TFile *fnoise;
  if (runNoise>-1){
    fnoise = new TFile(noiseFile.c_str(), "read");
    if (!fnoise) {
      cout << "Noise file " << noiseFile.c_str() << " does not exist. Please run findFilteredAverages and try again " << endl;
      return -1;
    }
  }


  TGraph *gfil[3];
  TGraph *gfiltmp[2];
  TGraph *power[2];
  TGraph *gnoise[2];

  bool saveCanvas = false;

  double lifeCathodeOnly, lifeApprox, lifetime;
  double QK, QA, QKcorr, QAcorr;
  double t1, t2, t3;
  int numAveraged;
  
  double rawVoltageRMSk = getRawVoltageRMS(basename, runNumber, chname[0]);
  double rawVoltageRMSa = getRawVoltageRMS(basename, runNumber, chname[1]);

  TFile *outLife = new TFile(outfile.c_str(), "recreate");
  TTree *metaTree = new TTree("metaTree", "Metadata tree");
  metaTree->Branch("run",             &runNumber         );
  metaTree->Branch("runNoise",        &runNoise          );
  metaTree->Branch("timestamp",       &timestamp         );
  metaTree->Branch("fields[3]",       fields             );
  metaTree->Branch("lampIntegral",    &lampIntegral      );
  metaTree->Branch("lampPeak",        &lampPeak          );
  metaTree->Branch("rawVoltageRMSk",   &rawVoltageRMSk   );
  metaTree->Branch("rawVoltageRMSa",   &rawVoltageRMSa   );

  metaTree->Fill();
  metaTree->Write();


  TTree *lifeTree = new TTree("lifeTree", Form("Lifetime info for Purity Monitor %d", whichPrM));
  lifeTree->Branch("run",             &runNumber       );
  lifeTree->Branch("numAveraged",     &numAveraged     );
  lifeTree->Branch("lifeCathodeOnly", &lifeCathodeOnly );
  lifeTree->Branch("lifeApprox",      &lifeApprox      );
  lifeTree->Branch("lifetime",        &lifetime        );
  lifeTree->Branch("QK",              &QK              );
  lifeTree->Branch("QA",              &QA              );
  lifeTree->Branch("QKcorr",          &QKcorr          );
  lifeTree->Branch("QAcorr",          &QAcorr          );
  lifeTree->Branch("t1",              &t1              );
  lifeTree->Branch("t2",              &t2              );
  lifeTree->Branch("t3",              &t3              );


  for (int inum=0; inum<nnum; inum++){

    // if(inum==0) saveCanvas=true;
    // else saveCanvas=false;

    double finalNumbers[2][3]; // [0 anode, 1 cathode] [0 amplitude, 1 start time, 2 peak time]
    
    for (int igraph=0; igraph<howManyGraphs[inum]; igraph++){
      
      for (int ich=0; ich<2; ich++){
	
	if (inum==0) gfil[ich] = (TGraph*)fin->Get(Form("gfil_%s", chname[ich].c_str()));
	else gfil[ich] = (TGraph*) fin->Get(Form("%s/gfil_%s_%d", howManyAvg[inum].c_str(), chname[ich].c_str(), igraph)); 
	

	if (runNoise>-1){
	  
	  if (inum==0) gnoise[ich] = (TGraph*)fnoise->Get(Form("gfil_%s", chname[ich].c_str()));
	  else gnoise[ich] = (TGraph*) fnoise->Get(Form("%s/gfil_%s_%d", howManyAvg[inum].c_str(), chname[ich].c_str(), igraph)); 

	  for (int ip=0; ip<gfil[ich]->GetN(); ip++){
	    gfil[ich]->GetY()[ip] -= gnoise[ich]->GetY()[ip];
	  }
	  
	  
	  if (inum==0){
	    outLife->cd();
	    gfil[ich]->Write(Form("gfin_%s", chname[ich].c_str()));

	    
	    // TGraph *gpower = FFTtools::makePowerSpectrumVoltsSeconds(gfil[ich]);
	    // gpower->Write(Form("gpow_fin_%d", ich));


	    if (ich==1){
	      TGraph *gdiff = new TGraph(gfil[0]->GetN());
	      for (int ip=0; ip<gfil[0]->GetN(); ip++) gdiff->SetPoint(ip, gfil[0]->GetX()[ip], gfil[0]->GetY()[ip]-gfil[1]->GetY()[ip]);
	      gdiff->Write("gdiff");
	    }

	  }
	}
      }

      if (inum==0){
	c->cd();
	TPaveText *pav = new TPaveText(0.25, 0.11, 0.99, 0.15, "NB NDC");
	pav->SetFillColor(kWhite);
	pav->AddText(elog);
	double min = TMath::MinElement(gfil[0]->GetN(), gfil[0]->GetY());
	double max = TMath::MaxElement(gfil[1]->GetN(), gfil[1]->GetY());
	if (max<5) max = 5.;
	gfil[0]->GetYaxis()->SetRangeUser(min*1.2, max*1.2);
	gfil[0]->SetTitle(Form("PrM%d, Filtered Averages and Noise subtracted;Time [s];Amplitude [mV]", whichPrM));
	gfil[0]->Draw("Al");
	gfil[1]->SetLineColor(kRed);
	gfil[1]->Draw("l");
	pav->Draw();
	c->Print(Form("%s/PurityMonitor%d_filAvg_subNoise.png", runname.c_str(), whichPrM));
      }
     
      double tlifetime[10];
      
      int ok = UsefulFunctions::calculateLifetime(gfil[0], gfil[1],  whichPrM-1, tTheory, tlifetime, saveCanvas);
      
      numAveraged=howManyAvgInt[inum];
 
      lifeCathodeOnly = -999.;
      lifeApprox = tlifetime[0];
      lifetime   = tlifetime[1];
      QK         = tlifetime[2];
      QA         = tlifetime[3];
      QKcorr     = tlifetime[4];
      QAcorr     = tlifetime[5];
      t1         = tlifetime[6];
      t2         = tlifetime[7];
      t3         = tlifetime[8];



      if (numAveraged==1000){

	lifeCathodeOnly = getLifetimeFromCathode(gfil[0], tlifetime);

	printf("Lifetime from Cathode [us] : %8.2f \n", lifeCathodeOnly*1e6);
	printf("Lifetime [us]              : %8.2f \n", lifetime*1e6);
	printf("QK [mV]                    : %8.2f \n", QK);
	printf("QA [mV]                    : %8.2f \n", QA);
	printf("t1 [us]                    : %8.2f \n", t1*1e6);
	printf("t2 [us]                    : %8.2f \n", t2*1e6);
	printf("t3 [us]                    : %8.2f \n", t3*1e6);

	fprintf(outFile, "Lifetime from Cathode [us] : %8.2f \n", lifeCathodeOnly*1e6);
	fprintf(outFile, "Lifetime [us]              : %8.2f \n", lifetime*1e6);
	fprintf(outFile, "QK [mV]                    : %8.2f \n", QK);
	fprintf(outFile, "QA [mV]                    : %8.2f \n", QA);
	fprintf(outFile, "t1 [us]                    : %8.2f \n", t1*1e6);
	fprintf(outFile, "t2 [us]                    : %8.2f \n", t2*1e6);
	fprintf(outFile, "t3 [us]                    : %8.2f \n", t3*1e6);

	if (inum==0){
	  c->cd();
	  gStyle->SetOptFit(1);
	  double min = TMath::MinElement(gfil[0]->GetN(), gfil[0]->GetY());
	  double max = TMath::MaxElement(gfil[1]->GetN(), gfil[1]->GetY());
	  if (max<5) max = 5.;
	  gfil[0]->GetYaxis()->SetRangeUser(min*1.2, max*1.2);
	  gfil[0]->SetTitle(Form("PrM%d, Filtered Averages and Noise subtracted;Time [ns];Amplitude [mV]", whichPrM));
	  gfil[0]->Draw("Al");
	  c->Print(Form("%s/PurityMonitor%d_cathodeOnlyFit.png", runname.c_str(), whichPrM));
	}

	
      }
      

      outLife->cd();
      lifeTree->Fill();

      cout << "\r " << igraph+1 << " / " << howManyGraphs[inum] << flush;

    }
      
    cout << endl;
  }

  outLife->cd();
  lifeTree->Write();
  outLife->Close();

  fclose (outFile);
  
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
  
  // cout << s << endl;

  string torm = "field";
  
  size_t position = s.find(torm.c_str());   

  s = s.substr(position + torm.size(), s.size()); 

  // cout << s << endl;

  std::replace( s.begin(), s.end(), '.', ' ');
  // cout << s << endl;
  
  std::stringstream ss(s);
  ss >> fields[0] >> fields[1] >> fields[2];
  
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


double getRawVoltageRMS(string basename, int irun, string chname){
  string filename = Form("%s/Run%03i/RootifiedFromBinary.root", basename.c_str(), irun);
  struct stat st;
 
  double rms =-1;
 
  if ( stat(filename.c_str(), &st)==0 ){
    
    TFile *ftemp = new TFile(filename.c_str(), "read");
    TGraph *gtemp = (TGraph*)ftemp->Get(Form("g_%s_1", chname.c_str()));
    rms = TMath::RMS(gtemp->GetN(), gtemp->GetY());
    delete gtemp;
    delete ftemp;
    
  }

  return rms;
}


Double_t greenFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec
                                                                                                                                                              
  double zero = par[4];
  double t = x[0]-zero;
  double tauEl = par[0]*1e-6;
  double T1 = par[1];
  double min = par[2];
  double tauLife = par[3]*1e-6;

  double y;
  double factor = 1/(1/tauLife-1/tauEl);
  double GQ = min*T1/factor/(exp(-T1/tauLife)-exp(-T1/tauEl));

  if (t>=T1) y = GQ/T1*factor*(exp(T1/tauEl)-exp(T1/tauLife))*exp(-t/tauEl-T1/tauLife);
  if (t>0 && t<T1) y = GQ*factor/T1*(exp(t/tauEl)-exp(t/tauLife))*exp(-t/tauEl-t/tauLife);
  if (t<=0) y = 0;

  return y;
}


double getLifetimeFromCathode(TGraph *g1, double parlife[10]){

  double_t min_at_x = parlife[6];
  
  TF1 *func = new TF1("func",greenFunction, g1->GetX()[0], g1->GetX()[g1->GetN()-1], 5);

  double tauLife = 50.;

  double tauEl = 300.;

  double factor = 1/(1/tauLife-1/tauEl);
  double minimum = parlife[2];
  double GQ = minimum/(exp(-min_at_x/tauEl)-1)*min_at_x/factor;
  double GQmod = minimum/(exp(-min_at_x/tauEl)-exp(-min_at_x/tauLife))*min_at_x/factor;
  // cout << "factor " << factor << endl;
  // cout << "GQ " << GQ << " \t GQmod " << GQmod << endl << endl;

  func->SetParameters(tauEl, min_at_x, minimum, tauLife, 4.*1e-6);

  //  func->FixParameter(0, tauEl);                                                                                                                                                
  func->SetParLimits(0,tauEl*0.7, tauEl*1.3);

  // //  func->FixParameter(1, min_at_x);                                                                                                                                          
  func->SetParLimits(1, min_at_x*0.9, min_at_x*1.1);
  func->SetParLimits(4, -2.e-6, 6.e-6);
  //func->FixParameter(4, 4.e-6);                                                                                                                                                  


  func->SetParName(0, "TauEl (#mus)");
  func->SetParName(1, "T1 (s)");
  func->SetParName(2, "minimum (V)");
  func->SetParName(3, "tauLife (#mus)");
  func->SetParName(4, "t0 (s)");

  func->SetLineColor(kMagenta);
  func->SetLineWidth(2);

  gStyle->SetOptStat();
  gStyle->SetOptFit();
  g1->Fit("func","RQ");
  func->Draw("same");

  double outpar[5];
  for (int i=0; i<4; i++) outpar[i] = func->GetParameter(i);

  double outMin = outpar[2]*(exp(-outpar[1]/(outpar[0]*1e-6))-exp(-outpar[1]/(outpar[3]*1e-6)) )/(outpar[1]*(1/(outpar[3]*1e-6)-1/(outpar[0]*1e-6)));

  double funcMin = func->Eval(outpar[1]);

  cout << "Minimum found by hand " << minimum << endl;
  cout << "Value of fitted function at T1 " << funcMin << endl;
  cout << "Minimum of fitted function " << func->GetMinimum() << endl;
  cout << "Minimum found from the fitted GQ " << outMin << endl;

  return outpar[3]*1e-6;


}
