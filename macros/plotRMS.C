#include <sys/stat.h>
void plotRMS(string basename="/data/PurityMonitor/Filling/"){

  int runNumber;

  int maxRun=0;
  FILE *frun;
  if (( frun = fopen(Form("%s/LastRun", basename.c_str()), "r")) ){
    fscanf(frun,"%d", &maxRun);
    printf("The last run number was %d\n", maxRun);
    fclose(frun);
  }

  double rms;
  double x[500], y[500];
  int n=0;

  for (int irun=50; irun<=maxRun;irun++){
    string filename = Form("%s/Run%03i/RootifiedFromBinary.root", basename.c_str(), irun);
    struct stat st;
    if (stat(filename.c_str(), &st)==0 ){
      
      TFile *ftemp = new TFile(filename.c_str(), "read");
      TGraph *gtemp = (TGraph*)ftemp->Get("g_ch4_1");
      if (!gtemp) continue;
      rms = TMath::RMS(gtemp->GetN(), gtemp->GetY());
      cout << filename << " " << irun << " " << rms << endl;
      x[n] = irun;
      y[n] = rms;
      n++;
      delete gtemp;
      delete ftemp;

    }
  }


  TGraph *g = new TGraph(n, x, y);

  g->Draw("Al");


}
