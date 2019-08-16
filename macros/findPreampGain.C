double getMinimum(TGraph *g1);
double getMaximum(TGraph *g1);

Double_t greenFunction(Double_t *x, Double_t *par);
Double_t greenFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec
  double t = x[0];
  double GQ = par[0];
  double tau = par[1]*1e-6;
  double rise = par[2]*1e-6;
  double shift = par[3]*1e-6;

  double heaviside;
  if (t<0) heaviside=0;
  else heaviside = 1;

  //double y = (-GQ)*exp(-t/(tau))*heaviside;

  double y;
  double tmp1 = GQ*(-1-erf((0+shift)/rise));
  double tmp2 = -GQ*exp(-0/tau);
  if (t>0) y = -GQ*exp(-t/tau);
  else y = GQ*(-1-erf((t+shift)/rise))*tmp2/tmp1;

  return y;
}


void findPreampGain(){

  bool isCathode=true;

  string ch;
  string chNice;
  if (isCathode){
    ch+="ch1";
    chNice+="Cathode";
  } else {
    ch+="ch3";
    chNice+="Anode";
  }

  // string name1 = "/unix/dune/purity/CERN/2019/Liquid/PrM1/Day4/AllFibres_wFilters_aCathode_bAnode_newResistors/Field_25.50.100Vcm."+ch+".traces_averages.root";
  // string name2 = "/unix/dune/purity/CERN/2019/Liquid/PrM1/Day4/AllFibres_wFilters_bCathode_aAnode_newResistors/Field_25.50.100Vcm."+ch+".traces_averages.root";


  // string name1 = "/data/PurityMonitor/GasTests/Run074/RawAverages_ch1.root"; <--- Linda's
  // string name2 = "/data/PurityMonitor/GasTests/Run075/RawAverages_ch2.root"; <--- Linda's

  // string preamp1 = "B"; <--- Linda's
  // string preamp2 = "A"; <--- Linda's

  // string whichprm = "1"; <--- Linda's

  //PrM1: ch1 = cathode, ch2 = anode; PrM2: ch4 = cathode, ch 5 = anode
  string name1 = "/data/PurityMonitor/GasTests/Run119/RawAverages_ch4.root";
  string name2 = "/data/PurityMonitor/GasTests/Run122/RawAverages_ch5.root"; 

  string preamp1 = "D";
  string preamp2 = "C";

  string whichprm = "2";

  TFile *f1 = new TFile(name1.c_str(), "read");
  TFile *f2 = new TFile(name2.c_str(), "read");

  TGraph *g1 = (TGraph*)f1->Get("justAvg");
  TGraph *g2 = (TGraph*)f2->Get("justAvg");

  for (int ip=0; ip<g1->GetN(); ip++){
    g1->GetX()[ip]*=1e-9; // convert from ns to s
    g2->GetX()[ip]*=1e-9;
    g1->GetY()[ip]*=1e-3; // convert from mV to V
    g2->GetY()[ip]*=1e-3;
  }

  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.6, 0.31, 0.89, 0.44);
  leg->AddEntry(g1, Form("PreAmp %s", preamp1.c_str() ), "l");
  leg->AddEntry(g2, Form("PreAmp %s", preamp2.c_str() ), "l");

  TCanvas *c = new TCanvas("c");
  g1->SetTitle(Form("PreAmps for PrM %s, %s signal;Time [s];Amplitude [V]", whichprm.c_str(), chNice.c_str()));
  double min1 = TMath::MinElement(g1->GetN(), g1->GetY());
  double max1 = TMath::MaxElement(g1->GetN(), g1->GetY());
  double min2 = TMath::MinElement(g2->GetN(), g2->GetY());
  double max2 = TMath::MaxElement(g2->GetN(), g2->GetY());
  
  double min, max;
  if (min1<min2) min = min1;
  else min = min2;
  if (max1>max2) max = max1;
  else max = max2;
  
  g1->GetYaxis()->SetRangeUser(min*1.1, max*1.1);

  g1->Draw("Al");
  g2->Draw("l");
  leg->Draw();

  double k1, k2;
  if (isCathode){
    k1 = getMinimum(g1);
    k2 = getMinimum(g2);
  } else {
    k1 = getMaximum(g1);
    k2 = getMaximum(g2);
  }

  double kmean1=0;
  double kmean2=0;
  double krms1=0;
  double krms2=0;

  double somek1[10];
  double somek2[10];

  int ngraphs=10;
  for (int i=0; i<ngraphs; i++){
    TGraph *g1t = (TGraph*)f1->Get(Form("avg100/gavg100_%i", i)); 
    if (isCathode) somek1[i] = getMinimum(g1t);
    else somek1[i] = getMaximum(g1t);
    delete g1t;
    TGraph *g2t = (TGraph*)f2->Get(Form("avg100/gavg100_%i", i));
    if (isCathode) somek2[i] = getMinimum(g2t);
    else somek2[i] = getMaximum(g2t);
    delete g2t;
    kmean1 += somek1[i]/ngraphs;
    kmean2 += somek2[i]/ngraphs;
  }
 
  for (int i=0; i<ngraphs; i++){
    krms1 += (somek1[i]-kmean1)*(somek1[i]-kmean1)/ngraphs;
    krms2 += (somek2[i]-kmean2)*(somek2[i]-kmean2)/ngraphs;
  }
  krms1 = TMath::Sqrt(krms1);
  krms2 = TMath::Sqrt(krms2);
  cout << k1 << " " << k2 << " " << k1/k2 << endl;


  double ratio, ratioerr;
  if (isCathode) ratio = kmean1/kmean2;
  else ratio = kmean2/kmean1;
  ratioerr = ratio*TMath::Sqrt(krms1*krms1/(kmean1*kmean1)+krms2*krms2/(kmean2*kmean2));

  TPaveText *paves = new TPaveText(0.6, 0.11, 0.89, 0.3, "ndc");
  paves->AddText(Form("Preamp %s : %4.3f +/- %4.3f", preamp1.c_str(), kmean1, krms1));
  paves->AddText(Form("Preamp %s : %4.3f +/- %4.3f", preamp2.c_str(), kmean2, krms2));
  if (isCathode) paves->AddText(Form("Ratio %s/%s : %4.3f +/- %4.3f", preamp1.c_str(), preamp2.c_str(), ratio, ratioerr));
  else  paves->AddText(Form("Ratio %s/%s : %4.3f +/- %4.3f", preamp2.c_str(), preamp1.c_str(), ratio, ratioerr));
  paves->SetFillColor(kWhite);
  paves->Draw();

  cout << "Preamp " << preamp1 << " : " << kmean1 << " +/- " << krms1 << endl;
  cout << "Preamp " << preamp2 << " : " << kmean2 << " +/- " << krms2 << endl;
  if (isCathode) cout << "Ratio " <<  preamp1 << "/" << preamp2 << " : " << ratio << " +/- " << ratioerr << endl;
  else cout << "Ratio " <<  preamp2 << "/" << preamp1 << " : " << ratio << " +/- " << ratioerr << endl;

  c->Print(Form("plots/filtered/ProtoDUNE_gasArgon_GainCalibrationForPrM%s_from%s.png", whichprm.c_str(), chNice.c_str()));


  TF1 *func = new TF1("func",greenFunction,0, 1e-3,4);
  
  //par[0] = Q*G, par[1]=tau_el                                                                                                               
  func->SetParameters(1., 300, 0.6, 0.6);
  func->SetParName(0, "Gain x Q ");
  func->SetParName(1, "Tau (#mus)");
  func->SetParLimits(1, 250, 300);
  func->SetParName(2, "Rise (#mus)");
  func->SetParName(3, "Shift (#mus)");

  //g1->Fit("func","R");

}


double getMinimum(TGraph *g1){

  int n1     = g1->GetN();
  double *x1 = g1->GetX();
  double *y1 = g1->GetY();
  double k1=0.;
  for (int i=0; i<n1; i++){
    if (y1[i]<k1 && x1[i]>0.02E-3){
      k1=y1[i];
    }
  }

  return k1;
}



double getMaximum(TGraph *g1){

  int n1     = g1->GetN();
  double *x1 = g1->GetX();
  double *y1 = g1->GetY();
  double k1=0.;
  for (int i=0; i<n1; i++){
    if (y1[i]>k1 && x1[i]>0.02E-3){
      k1=y1[i];
    }
  }

  return k1;
}
