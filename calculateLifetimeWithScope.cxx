#include "UsefulFunctions.h"
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


string howManyAvg[4] = {"justAvg", "avg50", "avg100", "avg200"};
int howManyGraphs[4] = {1, 20, 10, 5};

TGraph *smoothGraph(TGraph *g, int nnn);

int main(int argc, char *argv[]){

  int nnum = 3;
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
  TH1D *hratio[5];
  TH1D *hpurityPlus[5];
  TH1D *hpurityMinus[5];

  bool saveCanvas;
  
  for (int inum=0; inum<nnum; inum++){

    if(inum==0) saveCanvas=true;
    else saveCanvas=false;

    hpurity[inum]      = new TH1D (Form("hpurity_%d", inum), "", 1000, 0, 0.1);
    hratio[inum]       = new TH1D (Form("hratio_%d", inum), "", 1000, -2, 2);
    hpurityPlus[inum]  = new TH1D (Form("hpurityPlus_%d", inum), "", 1000, 0, 0.1);
    hpurityMinus[inum] = new TH1D (Form("hpurityMinus_%d", inum), "", 1000, 0, 0.1);
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
	
	//        g1 = smoothGraph(g1, smoothing[ich]);

	gdiff[ich] = (TGraph*)g1->Clone();

      }

      double lifetime[3]; // 0: approx formula for lifetime, 1: full formula for lifetime, 2: ratio of used for lifetime calculation
      
      int ok = UsefulFunctions::calculateLifetime(gdiff[1], gdiff[0],  whichPrM, tTheory, lifetime, saveCanvas);
      
      if (ok==1){
	hpurity[inum]->Fill(lifetime[1]);
	hratio[inum]->Fill(lifetime[2]);
	hpurityPlus[inum]->Fill(lifetime[1]*TMath::Log(lifetime[2])/TMath::Log(lifetime[2]*gainErrorPlus));
	hpurityMinus[inum]->Fill(lifetime[1]*TMath::Log(lifetime[2])/TMath::Log(lifetime[2]*gainErrorMinus));
      }
  
    }
    
  }
  
  FILE * pFile;
  
  cout << "Writing these info to " << outtxtfile << endl;


  double avgLifetime, avgRatio;
  double statErr, statErrRatio;
  double rmsLifetime, rmserrLifetime;
  
  double avgLifetimePlus;
  double statErrPlus;
  double avgLifetimeMinus;
  double statErrMinus;
  
  double absErrPlus;
  double absErrMinus;
  
  double totSyst, totError;
  double totSystRatio, totRatioError;
  
  pFile = fopen (outtxtfile.c_str(),"w");
  
  
  for (int inum=0; inum<nnum; inum++){

    avgLifetime = hpurity[inum]->GetMean();
    statErr = hpurity[inum]->GetMeanError();
    rmsLifetime = hpurity[inum]->GetRMS();
    rmserrLifetime = hpurity[inum]->GetRMSError();
  
    avgLifetimePlus  = hpurityPlus[inum]->GetMean();
    statErrPlus  = hpurityPlus[inum]->GetMeanError();
    avgLifetimeMinus = hpurityMinus[inum]->GetMean();
    statErrMinus = hpurityMinus[inum]->GetMeanError();

    absErrPlus  = TMath::Abs(avgLifetimePlus - avgLifetime);
    absErrMinus = TMath::Abs(avgLifetimeMinus - avgLifetime);

    totSyst  = absErrPlus ? absErrMinus : absErrPlus > absErrMinus;

    avgRatio     = hratio[inum]->GetMean();
    statErrRatio = hratio[inum]->GetMeanError();
    totSystRatio = avgRatio*(gainErrorPlus-gainErrorMinus)*0.5;
    
    if (inum==0){
      totError      = totSyst;
      totRatioError = totSystRatio;
    } else {
      totError      = TMath::Sqrt(totSyst*totSyst + statErr*statErr);
      totRatioError = TMath::Sqrt(totSystRatio*totSystRatio + statErrRatio*statErrRatio);
    }
    cout << "THIS IS OUR PURITY, STAT and SYST ERRORS: " << avgLifetime << " " << statErr << " " << totSyst << endl;


    // TF1 *fgaus = new TF1("fgaus", "gaus");
    // fgaus->SetParameter(1, avgLifetime);
    // hpurity[inum]->Fit("fgaus");
    // cout << "SOMENUMBERS " << fgaus->GetParameter(0) << " " << fgaus->GetParameter(1) << " " << fgaus->GetParameter(2) << endl;

    
    fprintf(pFile, "%s \n", howManyAvg[inum].c_str());
    fprintf(pFile, "%8.3e %8.3e %8.3e %8.3e \n", avgLifetime , totError,      statErr,      totSyst      );
    fprintf(pFile, "%8.3f %8.3f %8.3f %8.3f \n", avgRatio    , totRatioError, statErrRatio, totSystRatio );

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
    hratio[inum]->Write(Form("hratio_%d", inum));
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
