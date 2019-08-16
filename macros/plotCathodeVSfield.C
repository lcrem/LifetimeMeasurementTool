void plotCathodeVSfield(){

  string basename="/data/PurityMonitor/Filling/";

  // int runs[] = {48, 49, 50, 51, 52, 53};
  // int lowField[] = {70, 50, 40, 30, 20, 10};

  int runs[] = {164, 167, 166, 165, 163, 168};
  int lowField[] = { 0, 20, 30, 40, 50, 60};

  int nruns = sizeof(runs)/sizeof(runs[0]);

  TGraph *g[10];

  int color = 52;

  TCanvas *c = new TCanvas("c");

  TLegend *leg = new TLegend(0.5, 0.11, 0.89, 0.5, "Field");

  for (int irun=0; irun<nruns;irun++){
    
    TFile *fin = new TFile(Form("%s/Run%03i/PrM2_filtAvg.root", basename.c_str(), runs[irun]), "read");
    
    g[irun] = (TGraph*)fin->Get("gfil_ch4");
    g[irun]->SetLineColor(color);
    g[irun]->SetLineWidth(2);

    color +=8;
    

    if (irun==0){
      g[irun]->GetYaxis()->SetRangeUser(-40, 10);

      g[irun]->SetTitle("Cathode scan;Time [s]; Cathode signal [mV]");
      g[irun]->Draw("Al");
    } else g[irun]->Draw("l");

    leg->AddEntry(g[irun], Form("%i-%i-%i Vcm", lowField[irun], lowField[irun]*2, lowField[irun]*4), "l");

  }

  leg->Draw();

  c->Print("CathodeScan.png");
}
