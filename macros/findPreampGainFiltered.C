double getMinimum(TGraph *g1);
double getMaximum(TGraph *g1);

//Laura's functions
Double_t findMaxPosition(TGraph *G) {
  Double_t x, y;
  G->GetPoint(0, x, y);
  Double_t max_x = x, max_y = y;
  for(int i = 1; i < G->GetN(); i++)
    {
      G->GetPoint(i, x, y);
      //If you want to find the peak of a positive function (anode)
      //if(y > max_y) {
      //If you want to find the peak of a negative function (cathode)
      if(y < max_y) {
        max_x = x;
        max_y = y;
      }
    }
  return max_x;
}

Double_t greenFunction(Double_t *x, Double_t *par);
Double_t greenFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec
  double t = x[0];
  double tauEl = par[0]*1e-6;
  double T1 = par[1];
  double min = par[2];

  double heaviside_const;
  if (t >= (T1)) heaviside_const = 1.0;
  else heaviside_const = 0;
  double heaviside;
  if (t >= 0) heaviside = 1.0;
  else heaviside = 0;

  double y;
  double factor = -tauEl;
  double GQ = min*T1/factor/(1-exp(-T1/tauEl));

  if(t<=0) y = 0;
  if(t>=T1) y = GQ/T1*factor*(exp(T1/tauEl)-1)*exp(-t/tauEl);
  if (t>0 && t<T1) y = GQ*factor/T1*(1-exp(-t/tauEl));

  return y;
}

void findPreampGainFiltered(){

  string name1 = "/data/PurityMonitor/GasTests/Run129/PrM2_filtAvg.root"; // filtered cathode
  string name2 = "/data/PurityMonitor/GasTests/Run131/PrM2_filtAvg.root"; // filtered cathode
 
  string preamp1 = "D";
  string preamp2 = "C";

  string whichprm = "2"; //State PrM here
  int whichPrM = 2;      //And here!!! Until I remember how to convert a string to an int...

  bool isCathode=true;

  string ch[2];
  string chNice;
  if (isCathode){
    if(whichPrM==1){
      ch[0]+="ch1";
      ch[1]+="ch2";
    } else if (whichPrM==2){
      ch[0]+="ch4";
      ch[1]+="ch5";
    }
    chNice+="Cathode";
  } else {
    if(whichPrM==1){
      ch[0]+="ch2";
      ch[1]+="ch1";
    } else if (whichPrM==2){
      ch[0]+="ch5";
    ch[1]+="ch4";
    }
    chNice+="Anode";
  }
  
  cout << "First channel should be cathode: " << ch[0] << ", and second channel should be inverted cathode: " << ch[1] << endl;

  TFile *f1 = new TFile(name1.c_str(), "read");
  TFile *f2 = new TFile(name2.c_str(), "read");

  TGraph *g1 = (TGraph*)f1->Get(Form("gfil_%s", ch[0].c_str()));
  TGraph *g2 = (TGraph*)f2->Get(Form("gfil_%s", ch[1].c_str()));

  for (int ip=0; ip<g1->GetN(); ip++){
    // g1->GetX()[ip]*=1e-9; // convert from ns to s
    // g2->GetX()[ip]*=1e-9;
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

  double thousand_avg_min = getMinimum(g1);
  cout << "Cathode minimum from 1000 averaged waveforms: " << thousand_avg_min*1000 << " mV"  << endl;

  double x_value = 0;
  double y_value = 0;
  g1->GetPoint(0, x_value, y_value);
  cout << "X_value: " << x_value << ", and y_value: " << y_value << endl;

  double x2_value = 0;
  double y2_value = 0;
  g2->GetPoint(0, x2_value, y2_value);
  cout << "X2_value: " << x2_value << ", and y2_value: " << y2_value << endl;

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
    TGraph *g1t = (TGraph*)f1->Get(Form("avg100/gfil_%s_%i", ch[0].c_str(), i)); //if using filtered averages file
    if (isCathode) {
      somek1[i] = (getMinimum(g1t));
      //cout << somek1[i] << " ";
    }
    else somek1[i] = getMaximum(g1t);
    delete g1t;
    TGraph *g2t = (TGraph*)f2->Get(Form("avg100/gfil_%s_%i", ch[1].c_str(), i));
    if (isCathode){
      somek2[i] = (getMinimum(g2t));
      //cout << somek2[i] << endl;
    } else somek2[i] = getMaximum(g2t);
    delete g2t;
    kmean1 += somek1[i]/ngraphs;
    kmean2 += somek2[i]/ngraphs;
  }
 
  for (int i=0; i<ngraphs; i++){
    krms1 += (somek1[i]-kmean1)*(somek1[i]-kmean1)/((ngraphs-1)*ngraphs*1.0);
    krms2 += (somek2[i]-kmean2)*(somek2[i]-kmean2)/((ngraphs-1)*ngraphs*1.0);
    //cout << krms1 << " " << krms2 << endl;
  }
  krms1 = TMath::Sqrt(krms1);
  krms2 = TMath::Sqrt(krms2);
  cout << k1 << " " << k2 << " " << k1/k2 << " " << ngraphs << endl;


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

  //Fitting cathode 
  //3 param, 1st is tauEl in us, 2nd is x pos of the peak, 3rd is the height of the peak

  //Fit "standard" cathode:
  TF1 *func = new TF1("func",greenFunction,x_value, 0.94e-3, 3);
  double_t min_at_x = findMaxPosition(g1);
  cout << "Minimum at x: " << min_at_x << endl;
  cout << "Minimum is: " << getMinimum(g1) << endl;
  
  func->SetParameter(0, 300);
  func->SetParameter(1, min_at_x);
  func->SetParameter(2, getMinimum(g1));
  //func->SetParameter(3, x_value);
  //func->FixParameter(3, x_value);

  func->SetParName(0, "TauEl (#mus)");
  func->SetParName(1, "T1 (#mus)");
  func->SetParName(2, "min (V)");

  func->SetLineColor(kMagenta);
  func->SetLineWidth(2);
  g1->Fit("func","R");

  //Fit "inverted" cathode: 
  TF1 *func2 = new TF1("func2",greenFunction,x2_value, 0.94e-3,3);
  double_t min_at_x2 = findMaxPosition(g2);
  cout << "Minimum at x2: " << min_at_x2 << endl;
  cout << "Minimum is: " << getMinimum(g2) << endl;

  func2->SetParameter(0, 300);
  func2->SetParameter(1, min_at_x2);
  func2->SetParameter(2, getMinimum(g2));

  func2->SetParName(0, "TauEl (#mus)");
  func2->SetParName(1, "T1 (#mus)");
  func2->SetParName(2, "min (V)");

  func2->SetLineColor(kCyan);
  func2->SetLineWidth(2);
  g2->Fit("func2","R");

}

double getMinimum(TGraph *g1){

  int n1     = g1->GetN();
  double *x1 = g1->GetX();
  double *y1 = g1->GetY();
  double k1=0.;
  for (int i=0; i<n1; i++){
    if (y1[i]<k1 && x1[i]>0.002E-3){ 
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
    if (y1[i]>k1 && x1[i]>0.002E-3){
      k1=y1[i];
    }
  }

  return k1;
}
