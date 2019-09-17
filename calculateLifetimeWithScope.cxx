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
#include "TF1.h"

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

void getFields (string fieldname, double fields[3]);

double getLifetimeFromCathode(TGraph *g1, double parlife[10]);

Double_t greenFunction(Double_t *x, Double_t *par);


string howManyAvg[4] = {"justAvg", "avg50", "avg100", "avg200"};
int howManyGraphs[4] = {1, 20, 10, 5};
int howManyAvgInt[4] = {1000, 50, 100, 200};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 1;
  string basename, fieldname;
  int whichPrM=0;

  if((argc!=4)){
    std::cerr << "Usage : " << argv[0] << " [basename] [fieldname] [whichPrM]" << std::endl;
    return 1;
  } else {
    basename += argv[1];
    fieldname += argv[2];
    whichPrM = atoi(argv[3]);
  }

  string outfile    = basename + fieldname + "_signals_withErrors.root";
  string outtxtfile = basename + fieldname + "_lifetime_withErrors.txt";
  
  double fields[3], distance[3], tTheory[3];
  getFields (fieldname, fields);

  if (whichPrM%2==0){
    distance[0] = PrM1distance[0];
    distance[1] = PrM1distance[1];
    distance[2] = PrM1distance[2];
  } else {
    distance[0] = PrM2distance[0];
    distance[1] = PrM2distance[1];
    distance[2] = PrM2distance[2];
  }
  
  cout << "Output file is " << outfile << endl;
  
  string chname[2] = {"ch3", "ch4"};
  string chnamenice[2] = {"anode", "Cathode"};
  
  double timedelay=0;
  
  string stimedelay = basename + fieldname  + ".ch1.traces_averages.root";
  TFile *ftimedelay = new TFile (stimedelay.c_str(), "read");
  TGraph *gtimedelay = (TGraph*)ftimedelay->Get("justAvg");
  double *xTimeDelay = gtimedelay->GetX();
  double *yTimeDelay = gtimedelay->GetY();
  double baseY = yTimeDelay[0];
  double step = TMath::MaxElement(gtimedelay->GetN(), yTimeDelay)*0.9;
  double timeStep = xTimeDelay[1]-xTimeDelay[0];
  for (int ip=0; ip<gtimedelay->GetN(); ip++){
    if (yTimeDelay[ip]>(step)){
      timedelay = xTimeDelay[ip];
      break;
    }
  }
  ftimedelay->Close();
  cout << "The time delay is " << timedelay << endl;

  for (int i=0; i<3; i++){
    tTheory[i] = distance[i]/UsefulFunctions::ICARUSpolynomial(fields[i]);
    cout << distance[i] << " " << fields[i] << " " << tTheory[i] << " " << timeStep << " " << UsefulFunctions::getSmoothingNumber(timeStep, tTheory[i]) <<endl;
  }
  int smoothing[2];
  smoothing[0] =  UsefulFunctions::getSmoothingNumber(timeStep, tTheory[2]);
  smoothing[1] =  UsefulFunctions::getSmoothingNumber(timeStep, tTheory[0]);
  double newy[20000];
  double newx[20000];  
  
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  TH1D *hpurity[5];

  bool saveCanvas;

  double lifeCathodeOnly, lifeApprox, lifetime;
  double QK, QA, QKcorr, QAcorr;
  double t1, t2, t3;
  int numAveraged;
  
  TCanvas *c = new TCanvas("c");

  for (int inum=0; inum<nnum; inum++){

    if(inum==0) saveCanvas=true;
    else saveCanvas=false;

    hpurity[inum]= new TH1D (Form("hpurity_%d", inum), "", 1000, 0, 0.005);
    double finalNumbers[2][3]; // [0 anode, 1 cathode] [0 amplitude, 1 start time, 2 peak time]

    for (int igraph=0; igraph<howManyGraphs[inum]; igraph++){

      TGraph *gdiff[2];
      
      for (int ich=0; ich<2; ich++){
	string f1 = basename + fieldname  + "." + chname[ich] +".traces_averages.root";
        
	TFile *file1 = new TFile(f1.c_str(), "read");

	string gname ;
	if (howManyGraphs[inum]==1) gname+= "justAvg";
	else gname += Form("%s/g%s_%d", howManyAvg[inum].c_str(), howManyAvg[inum].c_str(), igraph);

	TGraph *g1 = (TGraph*)file1->Get(gname.c_str());
	cout << f1 << " " <<  Form("%s/g%s_%d", howManyAvg[inum].c_str(), howManyAvg[inum].c_str(), igraph) << endl;
	g1->SetName("g1");
        
	file1->Close();
        
	//	g1 = UsefulFunctions::translateGraph(g1, -timedelay);

	UsefulFunctions::zeroBaseline(g1);
	
	// Filter out everything above 200kHz                                                                                     
        gdiff[ich] = FFTtools::simplePassBandFilter(g1, 0., 100000000.);
	
	//	gdiff[ich] = smoothGraph(gdiff[ich], 5);
	//	gdiff[ich] = (TGraph*)g1->Clone();

      }

      double tlifetime[10];
      
      int ok = UsefulFunctions::calculateLifetime(gdiff[1], gdiff[0],  whichPrM, tTheory, tlifetime, saveCanvas);
      
      if (ok==1) hpurity[inum]->Fill(tlifetime[0]);

      numAveraged=howManyAvgInt[inum];

      lifeCathodeOnly = -999.;

      QK         = tlifetime[2];
      QA         = tlifetime[3];
      t1         = tlifetime[6];
      t2         = tlifetime[7];
      t3         = tlifetime[8];
      lifeApprox = tlifetime[0];
      lifetime   = tlifetime[1];
      QKcorr     = tlifetime[4];
      QAcorr     = tlifetime[5];

      if (numAveraged==1000){

        lifeCathodeOnly = getLifetimeFromCathode(gdiff[1], tlifetime);

        printf("Lifetime from Cathode [us] : %8.2f \n", lifeCathodeOnly*1e6);
        printf("Lifetime [us]              : %8.2f \n", lifetime*1e6);
        printf("QK [mV]                    : %8.2f \n", QK);
        printf("QA [mV]                    : %8.2f \n", QA);
        printf("t1 [us]                    : %8.2f \n", t1*1e6);
        printf("t2 [us]                    : %8.2f \n", t2*1e6);
        printf("t3 [us]                    : %8.2f \n", t3*1e6);

        if (inum==0){
          c->cd();
          gStyle->SetOptFit(1);
          double min = TMath::MinElement(gdiff[0]->GetN(), gdiff[0]->GetY());
          double max = TMath::MaxElement(gdiff[1]->GetN(), gdiff[1]->GetY());
          if (max<5) max = 5.;
	  //          gdiff[1]->GetYaxis()->SetRangeUser(min*1.2, max*1.2);
          gdiff[1]->SetTitle(Form("PrM%d, Filtered Averages;Time [s];Amplitude [V]", whichPrM));
          gdiff[1]->Draw("Al");
	  gdiff[0]->Draw("l");
	  c->Print(Form("%s_PurityMonitor%d_cathodeOnlyFit.png", (basename+fieldname).c_str(), whichPrM));
          c->Print(Form("%s_PurityMonitor%d_cathodeOnlyFit.root", (basename+fieldname).c_str(), whichPrM));
        }

      }
  
    }
    
  }
  
  FILE * pFile;
  
  cout << "Writing these info to " << outtxtfile << endl;
  
  
  pFile = fopen (outtxtfile.c_str(),"w");
  
  
  for (int inum=1; inum<nnum; inum++){

    double avgLifetime = hpurity[inum]->GetMean();
    double errLifetime = hpurity[inum]->GetMeanError();
    double rmsLifetime = hpurity[inum]->GetRMS();
    double rmserrLifetime = hpurity[inum]->GetRMSError();
  

    cout << "THIS IS OUR PURITY AND ERROR: " << avgLifetime << " " << errLifetime << endl;


    // TF1 *fgaus = new TF1("fgaus", "gaus");
    // fgaus->SetParameter(1, avgLifetime);
    // hpurity[inum]->Fit("fgaus");
    // cout << "SOMENUMBERS " << fgaus->GetParameter(0) << " " << fgaus->GetParameter(1) << " " << fgaus->GetParameter(2) << endl;

    
    fprintf(pFile, "%s \n", howManyAvg[inum].c_str());
    fprintf(pFile, "%8.3e %8.3e %8.3e %8.3e \n", avgLifetime , errLifetime, rmsLifetime, rmserrLifetime);
    // fprintf(pFile, "tK     : %12.4e \n",  tK  );
    // fprintf(pFile, "tGK    : %12.4e \n",  tGK );
    // fprintf(pFile, "tGA    : %12.4e \n",  tGA );
    // fprintf(pFile, "tA     : %12.4e \n",  tA  );
    // fprintf(pFile, "t1     : %12.4e \n",  t1  );
    // fprintf(pFile, "t2     : %12.4e \n",  t2  );
    // fprintf(pFile, "t3     : %12.4e \n",  t3  );
    // fprintf(pFile, "QA     : %12.4e \n",  QA  );
    // fprintf(pFile, "QK     : %12.4e \n",  QK  );
    // fprintf(pFile, "QA corr: %12.4e \n",  newQA  );
    // fprintf(pFile, "QK corr: %12.4e \n",  newQK  );
    // fprintf(pFile, "R      : %12.4e \n",  R   );
    // fprintf(pFile, "lifetime : %12.4e \n",  lifetime  );
    // fprintf(pFile, "lifetime2: %12.4e \n",  lifetime2 );
  
    out->cd();
    hpurity[inum]->Write(Form("hpurity_%d", inum));
  }

  
  fclose (pFile);

  return 0;

}
  


  


void getFields (string fieldname, double fields[3]){
  
  string s = fieldname.substr(fieldname.find("_") + 1) ;
  
  std::string ext("Vcm");
  if ( s != ext &&
       s.size() > ext.size() &&
       s.substr(s.size() - ext.size()) == "Vcm" )
    {
      // if so then strip them off
      s = s.substr(0, s.size() - ext.size());
    }
  std::replace( s.begin(), s.end(), '.', ' ');
  
  std::stringstream ss(s);
  ss >> fields[0] >> fields[1] >> fields[2];
  
  cout << " The three fields are " << fields[0] << " " << fields[1] << " " << fields[2] << " V/cm" << endl;
}

TGraph *smoothGraph(TGraph *g, int nnn){

  int n = g->GetN();
  double *x = g->GetX();
  double *y = g->GetY();
  double newy[100010];
  double newx[100010];

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
    newx[count] = x[i];
    count++;
  }

  TGraph *gnew = new TGraph(count, newx, newy);

  return gnew;

}


Double_t greenFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec                                                           
                                                                                                                                 \

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

  double tauLife = 500.;

  double tauEl = 100.;

  double factor = 1/(1/tauLife-1/tauEl);
  double minimum = parlife[2];
  double GQ = minimum/(exp(-min_at_x/tauEl)-1)*min_at_x/factor;
  double GQmod = minimum/(exp(-min_at_x/tauEl)-exp(-min_at_x/tauLife))*min_at_x/factor;
  // cout << "factor " << factor << endl;                                                                                         
  // cout << "GQ " << GQ << " \t GQmod " << GQmod << endl << endl;                                                                

  func->SetParameters(tauEl, min_at_x, minimum, tauLife, 0.);
                                                                                                                                 
  //  func->SetParLimits(0,tauEl*0.5, tauEl*1.5);

  func->SetParLimits(0, 10., 350.);

  func->SetParLimits(3, 0, 10000);

  // //  func->FixParameter(1, min_at_x);                                                                                        \
                                                                                                                                  
  func->SetParLimits(1, min_at_x*0.7, min_at_x*1.3);
  func->SetParLimits(4, -10.e-6, 50.e-6);
  //func->FixParameter(4, 4.e-6);                                                                                                \
                                           


  func->SetParName(0, "TauEl (#mus)");
  func->SetParName(1, "T1 (s)");
  func->SetParName(2, "minimum (V)");
  func->SetParName(3, "tauLife (#mus)");
  func->SetParName(4, "t0 (s)");

  func->SetLineColor(kMagenta);
  func->SetLineWidth(2);

  gStyle->SetOptStat();
  gStyle->SetOptFit();
  int status = g1->Fit("func","R");
  func->Draw("same");

  if (status!=0) return 0.;


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

