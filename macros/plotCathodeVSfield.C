Double_t ICARUSpolynomial(Double_t *x, Double_t *par);

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

  double x[100], t1[100], v[100], time;

  int nruns = sizeof(runs)/sizeof(runs[0]);

  TGraph *g[30];

  int color = 52;

  TCanvas *c = new TCanvas("c");

  TLegend *leg = new TLegend(0.7, 0.11, 0.89, 0.75, "Field");

  for (int irun=0; irun<nruns;irun++){
    
    TFile *fin = new TFile(Form("%s/Run%03i/PrM1_lifeInfo.root", basename.c_str(), runs[irun]), "read");
    
    TTree *tree = (TTree*)fin->Get("lifeTree");
    tree->SetBranchAddress("t1", &time);
    tree->GetEntry(0);


    g[irun] = (TGraph*)fin->Get("gfin_ch1");
    g[irun]->SetLineColor(color);
    g[irun]->SetLineWidth(2);

    if (irun%2==0) color +=5;
    else color +=4;

    if (irun==0){
      g[irun]->GetXaxis()->SetRangeUser(-0.05e-3, 0.4e-3);
      g[irun]->GetYaxis()->SetRangeUser(-130, 10);

      g[irun]->SetTitle("Cathode scan;Time [s]; Cathode signal [mV]");
      g[irun]->Draw("Al");
    } else g[irun]->Draw("l");

    leg->AddEntry(g[irun], Form("%i-%i-%i Vcm", lowField[irun], 0, 0), "l");

    x[irun] = lowField[irun];
    t1[irun] = time;
    v[irun]  = PrM1distance[0]/time;
    cout << x[irun] << " " << t1[irun] << " " << v[irun] << endl;
  }

  leg->Draw();

  c->Print("CathodeScan.png");
  
  gStyle->SetOptFit(1);

  TGraph *gtime = new TGraph(nruns, x, t1);
  gtime->SetMarkerStyle(22);
  TGraph *gv = new TGraph(nruns, x, v);
  gv->SetMarkerStyle(22);

  gv->SetTitle("Velocity;Field [V/cm];Velocity [m/s]");
  gv->GetYaxis()->SetRangeUser(0, 2000);
  gv->Draw("Ap");
  TF1 *funcUS = new TF1("funcUS", "pol5", 0, 700);
  gv->Fit(funcUS);

  TF1 *funcICARUS = new TF1("funcICARUS", ICARUSpolynomial, 0, 700, 1); 
  funcICARUS->Draw("l same");
  funcICARUS->SetLineColor(kBlue);

  double t1_50avg = 2.97142e-05;
  double t1_50avgerr = 1.39949e-07; 
  double t1_60avg = 2.55282e-05;
  double t1_60avgerr = 1.67714e-07; 
  double v1_50avg = PrM1distance[0]/t1_50avg;
  double v1_60avg = PrM1distance[0]/t1_60avg;
  

  TLine *l1 = new TLine(0, v1_50avg, 650, v1_50avg);
  l1->SetLineWidth(2);
  l1->SetLineColor(kCyan);

  TLine *l2 = new TLine(0, v1_60avg, 650, v1_60avg);
  l2->SetLineWidth(2);
  l2->SetLineColor(kViolet);

  l1->Draw("same");
  l2->Draw("same");

  TLegend *leg1 = new TLegend(0.6, 0.11, 0.89, 0.4);
  leg1->AddEntry(funcUS, "Our polynomial fit", "l");
  leg1->AddEntry(funcICARUS, "ICARUS polynomial", "l");
  leg1->AddEntry(l1,   "K-GK velocity at 50.100.200Vcm", "l");
  leg1->AddEntry(l2,   "K-GK velocity at 60.120.240Vcm", "l");
  leg1->Draw();


  cout << " 50 Vcm : " << v1_50avg << " , " << funcUS->GetX(v1_50avg) << " , " << funcICARUS->GetX(v1_50avg) << endl;
  cout << " 60 Vcm : " << v1_60avg << " , " << funcUS->GetX(v1_60avg) << " , " << funcICARUS->GetX(v1_60avg) << endl;


  c->SetGridx();
  c->SetGridy();

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
