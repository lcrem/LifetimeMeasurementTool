#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

// Filling ended roughly with run 341
// From run 342 Cryostat is full

void makePlotsVsTime(int whichPrM=2, double whichField=50, int firstRun=342,  int maxRun=-1, string addString="" ){

  string basename="/data/PurityMonitor/Filling/";
  int whichAvg=1000;		     

  string outdir = Form("plots/%s/Field%dVcm", addString.c_str(), int(whichField));

  system(Form("mkdir -p %s", outdir.c_str()));

  string title;

  double fields[3];  
  int runMeta, runTree;
  UInt_t timestamp;
  double lifeCathodeOnly, lifeApprox, lifetime;
  double QK, QA, QKcorr, QAcorr;
  double t1, t2, t3, lampIntegral, lampPeak;
  double rawVoltageRMSk, rawVoltageRMSa;
  int numAveraged;

  int ndiv;
  
  FILE *frun;
  if (maxRun==-1){
    if (( frun = fopen(Form("%s/LastRun", basename.c_str()), "r")) ){
      fscanf(frun,"%d", &maxRun);
      printf("The last run number was %d\n", maxRun);
      fclose(frun);
    }
  }

  TChain *lifeTree = new TChain("lifeTree");
  TChain *metaTree = new TChain("metaTree");
  
  for (int irun=firstRun; irun<=maxRun;irun++){
    string filename = Form("%s/Run%03i/PrM%i_lifeInfo.root", basename.c_str(), irun, whichPrM);
    struct stat st;
    
    // check if the bad flag file exists
    if (stat(Form("%s/Run%03i/BadFlag", basename.c_str(), irun), &st)==0) continue;
    if (stat(Form("%s/Run%03i/BadFlagPrM%d", basename.c_str(), irun, whichPrM), &st)==0) continue;

    if (stat(filename.c_str(), &st)==0 ){
      lifeTree->Add(filename.c_str());
      metaTree->Add(filename.c_str());
      cout << "Adding " << filename << endl;
    }
  }


  double varDouble[20];

  string var[]     = {"rawVoltageRMSk", "rawVoltageRMSa",  "lampPeak", "lampIntegral", "lifeCathodeOnly", "lifeApprox", "lifetime", "QK", "QA", "QKcorr", "QAcorr", "t1", "t2", "t3"};

  string niceVar[] = {"Raw Waveform RMS cathode", "Raw Waveform RMS anode", "Lamp Peak", "Lamp Integral", "Lifetime from Cathode fit", "Lifetime approximation", "Lifetime", "QK [mV]", "QA [mV]", "QK corrected [mV]", "QA corrected [mV]", "t1 [s]", "t2 [s]", "t3 [s]"};

  int nvar = sizeof(var)/sizeof(var[0]);

  cout << lifeTree->GetEntries() << endl;

  int runNoise;

  metaTree->SetBranchAddress("run",             &runMeta         );
  metaTree->SetBranchAddress("runNoise",        &runNoise        );
  metaTree->SetBranchAddress("timestamp",       &timestamp       );
  metaTree->SetBranchAddress("fields[3]",       fields           );
  metaTree->SetBranchAddress(var[0].c_str(),    &varDouble[0]    );
  metaTree->SetBranchAddress(var[1].c_str(),    &varDouble[1]    );
  metaTree->SetBranchAddress(var[2].c_str(),    &varDouble[2]    );
  metaTree->SetBranchAddress(var[3].c_str(),    &varDouble[3]    );

  metaTree->BuildIndex("run");

  lifeTree->SetBranchAddress("run",             &runTree         );
  lifeTree->SetBranchAddress("numAveraged",     &numAveraged     );
  for (int ivar=4; ivar<nvar; ivar++){
    lifeTree->SetBranchAddress(var[ivar].c_str(),    &varDouble[ivar]    );
  }


  int count=0;

  std::vector<double> xValues;
  std::vector<double> yValues[20];

  FILE * outFile;
  outFile = fopen (Form("%s/PrM%i_summary.csv", outdir.c_str(), whichPrM),"w");
  
  fprintf(outFile, "Run, Timestamp, Field K-GK [Vcm], Field GK-GA [Vcm], Field [Vcm], ");

  for (int ivar=0; ivar<nvar; ivar++) fprintf(outFile, "%s, ", niceVar[ivar].c_str());
  fprintf(outFile, "\n");
  
  bool doneTitle=false;

  
  for (int ie=0; ie<lifeTree->GetEntries(); ie++){
    lifeTree->GetEntry(ie);

    if (numAveraged != whichAvg ) continue;

    metaTree->GetEntryWithIndex(runTree);
    
    if (fields[0] != whichField) continue;
    if (fields[1] != whichField*2) continue;
    if (fields[2] != whichField*4) continue;
    if (!doneTitle) {
      title+=Form("PrM %d, field %i-%i-%i V/cm", whichPrM, int(fields[0]), int(fields[1]), int(fields[2]));
      doneTitle=true;
    }
    
    xValues.push_back(timestamp);

    for (int ivar=0; ivar<nvar; ivar++){
      yValues[ivar].push_back(TMath::Abs(varDouble[ivar]));
      //      cout << niceVar[ivar] << " " << varDouble[ivar] << " " << yValues[ivar].at(count) << endl;
    }

    cout << runMeta << " " << runTree << " " << runNoise << " " << timestamp << " " << yValues[7].at(count) << " " << yValues[4].at(count) << " " << yValues[6].at(count) << endl;

    fprintf(outFile, "%d, %d, %d, %d, %d, ", runMeta, timestamp, int(fields[0]), int(fields[1]), int(fields[2]));
    for (int ivar=0; ivar<nvar; ivar++) fprintf(outFile, "%8.3e, ", yValues[ivar].at(count) );
    fprintf(outFile, "\n");
    count++;

  }

  
  fclose (outFile);

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "", 1000, 500); 

  TFile *fout = new TFile(Form("%s/PrM%i_LatestSummaryGraphs.root", outdir.c_str(), whichPrM), "recreate");


  double firstDay = floor(xValues[0]/(3600*24))*3600*24 + 600;
  
  double lastDay = (floor(xValues[count-1]/(3600*24)) + 1 )*3600*24 + 600;

  ndiv =  (lastDay-firstDay)/(3600*24)/2;

  double ymin, ymax;

  for (int ivar=0; ivar<nvar; ivar++){
    
    TGraph *gK = new TGraph(count, &xValues[0], &yValues[ivar][0]);
  
    ymin = TMath::MinElement(gK->GetN(), gK->GetY());
    ymax = TMath::MaxElement(gK->GetN(), gK->GetY());

    TH2D *htemp = new TH2D("htemp", "", ndiv, firstDay, lastDay, 10, ymin*0.99, ymax*1.01);

    //    ndiv = floor (0.5 + ( gK->GetXaxis()->GetBinUpEdge(gK->GetXaxis()->GetNbins()) - gK->GetXaxis()->GetBinLowEdge(1))/(3600*24) );

    htemp->SetTitle(Form("%s: %s;; %s", title.c_str(), niceVar[ivar].c_str(), niceVar[ivar].c_str()));
    
    htemp->GetXaxis()->SetTimeDisplay(1);
    htemp->GetXaxis()->SetTimeOffset(0, "GMT");
    //    htemp->GetXaxis()->SetTimeFormat("#splitline{%d}{%b %H}");
    htemp->GetXaxis()->SetTimeFormat("%d-%b");
    gK->SetMarkerStyle(8);

    htemp->GetXaxis()->SetNdivisions(ndiv, 2, 1, kFALSE);

    htemp->Draw();
    gK->Draw("lp");
    
    c->SetGridx();
    c->SetGridy();

    c->Print(Form("%s/PrM%i_%s.png", outdir.c_str(), whichPrM, var[ivar].c_str()));
    
    fout->cd();
    gK->Write(var[ivar].c_str());
    delete gK;
    delete htemp;
  }


  
}
