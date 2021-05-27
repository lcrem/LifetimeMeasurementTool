Double_t ICARUSpolynomial(Double_t *x, Double_t *par);
Double_t Schottky(Double_t *x, Double_t *par);

void plotCathodeVSfield(){

  string basename="/data/PurityMonitor/Filling/";

  // int runs[] = {48, 49, 50, 51, 52, 53};
  // int lowField[] = {70, 50, 40, 30, 20, 10};

  int runs[] =     {630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640};
  int lowField[] = { 20,  40,  50,  60,  80, 100, 120, 150, 200, 400, 600};  

  // First Purity Monitor machined at UCL (PrM1)
  const Double_t PrM1distance[3] = {0.01823, 0.16424, 0.00985};
  // ICARUS Purity Monitor refurbished (PrM2)
  const Double_t PrM2distance[3] = {0.018, 0.15, 0.01};

  double kDistance[2] = {0.01823, 0.018};

  double x[2][100], t1[2][100], v[2][100], time, qk[2][100], q0k[2][100], QKvalue, QKcorr, errx[2][100], erry[2][100];

  int nruns = sizeof(runs)/sizeof(runs[0]);

  TGraph *g[2][30];

  int color = 52;

  TCanvas *c = new TCanvas("c");

  TLegend *leg = new TLegend(0.7, 0.11, 0.89, 0.75, "Field");

  string chname[2] = {"ch1", "ch4"};

  for (int iprm=0; iprm<2; iprm++){

    for (int irun=0; irun<nruns;irun++){
      
      TFile *fin = new TFile(Form("%s/Run%03i/PrM%d_lifeInfo.root", basename.c_str(), runs[irun], iprm+1), "read");
    
      TTree *tree = (TTree*)fin->Get("lifeTree");
      tree->SetBranchAddress("t1", &time);
      tree->SetBranchAddress("QK", &QKvalue);
      tree->SetBranchAddress("QKcorr", &QKcorr);
      tree->GetEntry(0);
      
      
      g[iprm][irun] = (TGraph*)fin->Get(Form("gfin_%s", chname[iprm].c_str()));
      g[iprm][irun]->SetLineColor(color);
      g[iprm][irun]->SetLineWidth(2);

      if (irun%2==0) color +=5;
      else color +=4;
      
      if (irun==0){
	g[iprm][irun]->GetXaxis()->SetRangeUser(-0.05e-3, 0.4e-3);
	g[iprm][irun]->GetYaxis()->SetRangeUser(-130, 10);

	g[iprm][irun]->SetTitle("Cathode scan;Time [s]; Cathode signal [mV]");
	g[iprm][irun]->Draw("Al");
	color = 52;
      } else g[iprm][irun]->Draw("l");
      
      if (iprm==0) leg->AddEntry(g[iprm][irun], Form("%i-%i-%i Vcm", lowField[irun], 0, 0), "l");
      
      x[iprm][irun] = lowField[irun];
      t1[iprm][irun] = time;
      v[iprm][irun]  = kDistance[iprm]/time;
      qk[iprm][irun] = -QKvalue;
      q0k[iprm][irun]= -QKcorr;
      errx[iprm][irun]=0;
      erry[iprm][irun]=1;
      cout << x[iprm][irun] << " " << t1[iprm][irun] << " " << v[iprm][irun] << endl;
    }
    
    leg->Draw();

    c->Print(Form("CathodeScan_PrM%d_1.png", iprm+1));
  
  }

  gStyle->SetOptFit(1);

  TGraphErrors *gqk1 = new TGraphErrors(nruns, x[0], qk[0], errx[0], erry[0]);
  TGraphErrors *gq0k1 = new TGraphErrors(nruns, x[0], q0k[0], errx[0], erry[0]);
  TGraph *gtime1 = new TGraph(nruns, x[0], t1[0]);
  TGraph *gv1 = new TGraph(nruns, x[0], v[0]);

  TGraphErrors *gqk2 = new TGraphErrors(nruns, x[1], qk[1], errx[1], erry[1]);
  TGraphErrors *gq0k2 = new TGraphErrors(nruns, x[1], q0k[1], errx[1], erry[1]);
  TGraph *gtime2 = new TGraph(nruns, x[1], t1[1]);
  TGraph *gv2 = new TGraph(nruns, x[1], v[1]);

  gtime1->SetMarkerStyle(22);
  gv1->SetMarkerStyle(22);
  gqk1->SetMarkerStyle(22);
  gtime2->SetMarkerStyle(22);
  gv2->SetMarkerStyle(22);
  gqk2->SetMarkerStyle(22);
  
  gqk2->SetMarkerColor(kRed);
  gq0k2->SetMarkerColor(kRed);
  gtime2->SetMarkerColor(kRed);
  gv2->SetMarkerColor(kRed);

  gv1->SetTitle("Velocity;Field [V/cm];Velocity [m/s]");
  gv1->GetYaxis()->SetRangeUser(0, 2000);
  gv1->Draw("Ap");
  gv2->Draw("p");


  cout << "Velocity fit to PrM1" << endl;
  TF1 *funcUS1 = new TF1("funcUS1", "pol5", 0, 1000);
  funcUS1->SetLineColor(kBlack);
  gv1->Fit(funcUS1);

  cout << "Velocity fit to PrM2" << endl;
  TF1 *funcUS2 = new TF1("funcUS2", "pol5", 0, 1000);
  gv2->Fit(funcUS2);
  funcUS2->SetLineColor(kRed);

  TF1 *funcICARUS = new TF1("funcICARUS", ICARUSpolynomial, 0, 1000, 1); 
  funcICARUS->Draw("l same");
  funcICARUS->SetLineColor(kBlue);

  double t1_50avg = 2.97142e-05;
  double t1_50avgerr = 1.39949e-07; 
  double t1_60avg = 2.55282e-05;
  double t1_60avgerr = 1.67714e-07; 
  double v1_50avg = kDistance[0]/t1_50avg;
  double v1_60avg = kDistance[0]/t1_60avg;
  

  // TLine *l1 = new TLine(0, v1_50avg, 650, v1_50avg);
  // l1->SetLineWidth(2);
  // l1->SetLineColor(kCyan);

  // TLine *l2 = new TLine(0, v1_60avg, 650, v1_60avg);
  // l2->SetLineWidth(2);
  // l2->SetLineColor(kViolet);

  // l1->Draw("same");
  // l2->Draw("same");

  TLegend *leg1 = new TLegend(0.6, 0.11, 0.89, 0.25);
  leg1->AddEntry(funcUS1, "pol5 fit to PrM1", "l");
  leg1->AddEntry(funcUS2, "pol5 fit to PrM2", "l");
  leg1->AddEntry(funcICARUS, "ICARUS polynomial", "l");
  // leg1->AddEntry(l1,   "K-GK velocity at 50.100.200Vcm", "l");
  // leg1->AddEntry(l2,   "K-GK velocity at 60.120.240Vcm", "l");
  leg1->Draw();


  cout << " 50 Vcm : " << v1_50avg << " , " << funcUS1->GetX(v1_50avg) << " , " << funcICARUS->GetX(v1_50avg) << endl;
  cout << " 60 Vcm : " << v1_60avg << " , " << funcUS1->GetX(v1_60avg) << " , " << funcICARUS->GetX(v1_60avg) << endl;


  c->SetGridx();
  c->SetGridy();

  c->Update();
  TPaveStats *stat = (TPaveStats*)(gv1->FindObject("stats"));
  TPaveStats *stat1 = (TPaveStats*)(gv2->FindObject("stats"));
  if(stat && stat1) {
    stat->SetTextColor(kBlack);
    stat1->SetTextColor(kRed);
    float height = stat1->GetY2NDC() - stat1->GetY1NDC();
    stat1->SetY1NDC(stat->GetY1NDC() - height);
    stat1->SetY2NDC(stat->GetY1NDC() );
    stat1->Draw();
  }  

  c->Print("CathodeScan_2.png");


  gqk1->GetYaxis()->SetRangeUser(0, 130);

  gqk1->SetMarkerStyle(22);

  gqk1->SetTitle(";Field [V/cm];Cathode Amplitude");
  gqk1->Draw("Ape");
  gqk2->Draw("pe");

  cout << "QK fit to PrM1" << endl;
  TF1 *funcUSk1 = new TF1("funcUSk1", Schottky, 0, 300, 2);
  funcUSk1->SetLineColor(kBlack);
  gqk1->Fit(funcUSk1, "r");

  cout << "QK fit to PrM2" << endl;
  TF1 *funcUSk2 = new TF1("funcUSk2", Schottky, 0, 300, 2);
  funcUSk2->SetLineColor(kRed);
  gqk2->Fit(funcUSk2, "r");

  c->Update();
  TPaveStats *statt = (TPaveStats*)(gqk1->FindObject("stats"));
  TPaveStats *statt1 = (TPaveStats*)(gqk2->FindObject("stats"));

  if(statt && statt1) {
    statt->SetTextColor(kBlack);
    statt1->SetTextColor(kRed);
    float height = statt1->GetY2NDC() - statt1->GetY1NDC();
    statt1->SetY1NDC(statt->GetY1NDC() - height);
    statt1->SetY2NDC(statt->GetY1NDC() );
    statt1->Draw();
  }

  c->Print("CathodeScan_3.png");

  gtime1->SetMarkerStyle(22);
  gtime2->SetMarkerStyle(22);
  gtime1->GetYaxis()->SetRangeUser(0, 0.13e-3);
  gtime1->SetTitle(";Field [V/cm];t1 [s]");
  gtime1->Draw("Ape");
  gtime2->Draw("pe");
  c->Print("CathodeScan_4.png");
  

  TFile *fout = new TFile ("plots/CathodeScanFunctions.root", "recreate");
  funcUSk1->Write("funcUSqk1");
  funcUSk2->Write("funcUSqk2");
  funcUS1->Write("funcUSv1");
  funcUS2->Write("funcUSv2");
  funcICARUS->Write("funcICARUS");
  gqk1->Write("gqk1");
  gqk2->Write("gqk2");

  fout->Close();
}

Double_t Schottky(Double_t *x, Double_t *par){

  return par[0] + par[1]*TMath::Sqrt(x[0]);
}


Double_t ICARUSpolynomial(Double_t *x, Double_t *par){

  // transform E in kV/cm
  double E = x[0]*1e-3;
  
  Double_t Et = 0.5;
  Double_t T  = 89;      // Kelvin
  Double_t T0 = 90.371;  // Kelvin
  Double_t p0 = -0.03229;
  Double_t p1 = 6.231;
  Double_t p2 = -10.62;
  Double_t p3 = 12.74;
  Double_t p4 = -9.112;
  Double_t p5 = 2.83;
  Double_t w1 = -0.01481;
  Double_t w2 = -0.0075;
  Double_t w3 = 0.141;
  Double_t w4 = 12.4;
  Double_t w5 = 1.627;
  Double_t w6 = 0.317;

  Double_t K1 = p0 + p1*Et + p2*Et*Et + p3*Et*Et*Et + p4*Et*Et*Et*Et + p5*Et*Et*Et*Et*Et;
  Double_t K2 = ( w1*(T-T0) + 1 )*( w3*Et*TMath::Log(1+w4/Et) + w5*TMath::Power(Et,w6) ) + w2*(T-T0);

  Double_t vE = (p0 + p1*E + p2*E*E + p3*E*E*E + p4*E*E*E*E + p5*E*E*E*E*E)*(K2/K1);

  // convert vE from mm/us to m/s

  vE *= (1e-3/1e-6);
  
  return vE;

}
