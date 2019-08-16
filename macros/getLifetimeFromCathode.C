//This macro is designed to fit a cathode signal and extract the lifetime from it, when the cathode and grid cathode only are submerged. It works on filtered waveforms (although is easily adaptable for unfiltered) and fits one purity monitor at a time.
#include <TROOT.h>
#include <TStyle.h>

double getMinimum(TGraph *g1);

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

Double_t greenFunction2(Double_t *x, Double_t *par);
Double_t greenFunction2(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec                                                                                          
  double zero = par[4];                                      
  //double t = x[0];
  double t = x[0]-zero;                 
  double tauEl = par[0]*1e-6;
  double T1 = par[1];
  double GQ = par[2];
  double tauLife = par[3]*1e-6;
  double factor = 1/(1/tauLife-1/tauEl);

  double y;

  if (t>=T1) y = GQ/T1*factor*(exp(T1/tauEl)-exp(T1/tauLife))*exp(-t/tauEl-T1/tauLife);
  if (t>0 && t<T1) y = GQ*factor/T1*(exp(t/tauEl)-exp(t/tauLife))*exp(-t/tauEl-t/tauLife);
  //LINDA if (t>=T1) y = GQ*(exp(T1/tauEl)-exp(T1/tauLife))*exp(-t/tauEl-T1/tauLife);
  //LINDA if (t>0 && t<T1) y = GQ*(exp(t/tauEl)-exp(t/tauLife))*exp(-t/tauEl-t/tauLife);
  if (t<=0) y = 0;

  return y;
}

Double_t ilsooFunction(Double_t *x, Double_t *par);
Double_t ilsooFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec                                                                                                             
  double t = x[0];
  double tauEl = par[0]*1e-6;
  double T1 = par[1];
  double GQ = par[2];
  double tauLife = par[3]*1e-6;
  // double boff = par[4];

  double y;

  if (t>=T1) y = (GQ*(exp(T1/tauEl)-exp(T1/tauLife))) *exp(-(t-T1)/tauEl);
  if (t>0 && t<T1) y = GQ*(exp(-t/tauEl)-exp(-t/tauLife));
  if (t<=0) y = 0;

  return y;
}

void getLifetimeFromCathode(){

  //string name1 = "/data/PurityMonitor/ColdGas/Run025/PrM2_filtAvg.root"; // filtered cathode
  // string name1 = "/data/PurityMonitor/Filling/Run018/run018_prm2_subtracted.root"; // background subtracted cathode          
  string name1 = "/data/PurityMonitor/Filling/Run295/PrM2_lifeInfo.root";

  string preamp1 = "D"; //preamp B for PrM1 and preamp D for PrM2

  string whichprm = "2"; //State PrM here                                                                                                                                                           
  int whichPrM = atoi(whichprm.c_str());      //And here!!! Until I remember how to convert a string to an int...    

  string ch;
  if (whichPrM==1) ch = "ch1";
  if (whichPrM==2) ch = "ch4";

  TFile *f1 = new TFile(name1.c_str(), "read");
  //TGraph *g1 = (TGraph*)f1->Get(Form("gfil_%s", ch.c_str()));
  // TGraph *g1 = (TGraph*)f1->Get("Graph");  
  TGraph *g1 = (TGraph*)f1->Get("gfin_ch4");  

  // // subtracting anode
  // TGraph *g5 = (TGraph*)f1->Get("gfin_ch5");  
  // for (int ip=0; ip<g1->GetN(); ip++) g1->GetY()[ip]-=g5->GetY()[ip];


  for (int ip=0; ip<g1->GetN(); ip++){
    // g1->GetX()[ip]*=1e-9; // convert from ns to s
    g1->GetY()[ip]*=1e-3; // convert from mV to V
  }

  TCanvas *c = new TCanvas("c");
  g1->SetTitle(Form("PrM %s cathode;Time [s];Amplitude [V]", whichprm.c_str()));
  double min1 = TMath::MinElement(g1->GetN(), g1->GetY());
  double max1 = TMath::MaxElement(g1->GetN(), g1->GetY());

  g1->GetYaxis()->SetRangeUser(min1*1.1, max1*5.);

  g1->Draw("Al");
  
  //For the fitting function:
  double thousand_avg_min = getMinimum(g1);
  cout << "Cathode minimum from 1000 averaged waveforms: " << thousand_avg_min*1000 << " mV"  << endl;

  //3 param, 1st is tauEl in us, 2nd is x pos of the peak, 3rd is the height of the peak, 4th is tauLife in us                                               

  double_t min_at_x = findMaxPosition(g1);
  TF1 *func = new TF1("func",greenFunction, 0, g1->GetX()[g1->GetN()-1], 5);

  cout << "Minimum at x: " << min_at_x << endl;
  cout << "Minimum is: " << thousand_avg_min << endl;

  double tauLife = 50.; 
  
  double tauEl;
  if (whichPrM == 2) {
    tauEl = 303.;
  } 
  if (whichPrM == 1) {
    tauEl = 262.;
  }     
  
  double factor = 1/(1/tauLife-1/tauEl);
  double minimum = thousand_avg_min;
  double GQ = minimum/(exp(-min_at_x/tauEl)-1)*min_at_x/factor;
  double GQmod = minimum/(exp(-min_at_x/tauEl)-exp(-min_at_x/tauLife))*min_at_x/factor;
  cout << "factor " << factor << endl;
  cout << "GQ " << GQ << " \t GQmod " << GQmod << endl << endl;

  func->SetParameters(tauEl, min_at_x, minimum, tauLife, 4.*1e-6);
  
  //  func->FixParameter(0, tauEl);
  func->SetParLimits(0,tauEl*0.7, tauEl*1.3);

  // //  func->FixParameter(1, min_at_x);
  func->SetParLimits(1, min_at_x*0.9, min_at_x*1.1);

  // func->SetParLimits(2, minimum*0.8, minimum*1.2);

  // //  func->SetParameter(3, tauLife);
  // func->SetParLimits(3, 0, 1000);
  func->SetParLimits(4, 2.e-6, 6.e-6);
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
  g1->Fit("func","R");
  func->Draw("same");

  //LAURA
  //This is to get tauEl, min_at_x, and minimum
  //double *x = graph->GetX(); double *y = graph->getY();
  //dx = x[0];
  //for ( int i = 0; i<g1->GetN(); i++) x[i] = x[i]-dx;
  

  // double zero = 4.*1e-6;
  // TF1 *func1 = new TF1("func1",greenFunction, min_at_x, g1->GetX()[g1->GetN()-1], 5);
  // func1->SetParameters(tauEl, min_at_x, minimum, tauLife, zero);
  // func1->SetParameter(0, 270);
  // func1->FixParameter(3, 1*1e10);

  // func1->SetParName(0, "TauEl (#mus)");
  // func1->SetParName(1, "T1 (s)");
  // func1->SetParName(2, "min (V)");
  // func1->SetParName(3, "tauLife (#mus)");

  // func1->SetLineColor(kBlack);
  // func1->SetLineWidth(5);
  // g1->Fit("func1","R");
  // func1->Draw("same");

  // tauEl = func1->GetParameter(0);
  // min_at_x = func1->GetParameter(1);
  // minimum = func1->GetParameter(2);

  // TF1 *func2 = new TF1("func2",greenFunction, min_at_x, g1->GetX()[g1->GetN()-1], 4);
  // func2->SetParameters(tauEl, min_at_x, minimum, tauLife, zero);
  
  // func2->SetParName(0, "TauEl (#mus)");
  // func2->SetParName(1, "T1 (s)");
  // func2->SetParName(2, "min (V)");
  // func2->SetParName(3, "tauLife (#mus)");

  // func2->SetLineColor(kRed);
  // func2->SetLineWidth(3);
  // g1->Fit("func2","R");
  // func2->Draw("same");

  // TF1 *func3 = new TF1("func3",greenFunction, 0, min_at_x, 4);
  // func3->SetParameters(tauEl, min_at_x, minimum, tauLife, zero);
  // func3->FixParameter(0, func2->GetParameter(0));
  // func3->FixParameter(1, func2->GetParameter(1));
  // func3->FixParameter(2, func2->GetParameter(2));
  
  // func3->SetParName(0, "TauEl (#mus)");
  // func3->SetParName(1, "T1 (s)");
  // func3->SetParName(2, "min (V)");
  // func3->SetParName(3, "tauLife (#mus)");

  // func3->SetLineColor(kGreen);
  // func3->SetLineWidth(2);
  // g1->Fit("func3","R");
  // func3->Draw("same");

  double outpar[5];
  for (int i=0; i<4; i++) outpar[i] = func->GetParameter(i);
  
  //  double GQ = minimum/(exp(-min_at_x/tauEl)-1)*min_at_x/factor;
  //  minimum = GQ*(exp(-min_at_x/tauEl)-1)*factor/min_at_x

  //LINDA double outMin = outpar[2]*(1/(outpar[3]*1e-6)-1/(outpar[0]*1e-6))*(exp(-outpar[1]/(outpar[0]*1e-6))-exp(-outpar[1]/(outpar[3]*1e-6)) );
  double outMin = outpar[2]*(exp(-outpar[1]/(outpar[0]*1e-6))-exp(-outpar[1]/(outpar[3]*1e-6)) )/(outpar[1]*(1/(outpar[3]*1e-6)-1/(outpar[0]*1e-6)));

  double funcMin = func->Eval(outpar[1]);

  cout << "Minimum found by hand " << minimum << endl;
  cout << "Value of fitted function at T1 " << funcMin << endl;
  cout << "Minimum of fitted function " << func->GetMinimum() << endl;
  cout << "Minimum found from the fitted GQ " << outMin << endl;

  TFile *fout = new TFile("ForLaura.root", "recreate");
  c->Write();


  //  c->Print(Form("plots/filtered/cathode_signal_PrM%s.png", whichprm.c_str()));
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
