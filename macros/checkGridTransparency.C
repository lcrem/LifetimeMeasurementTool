#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>

Double_t ICARUSpolynomial(Double_t E);

void checkGridTransparency(string infolder="/data/PurityMonitor/Filling/", int firstRun=600, int lastRun=608){

  string outname="GridTransparency_betterpurity";
    
  double errQK = 2;
  double errQA = 2;

  double tauLife = 540.e-6;

  double distance[2][3];

  distance[0][0] =  0.01823;
  distance[1][0] =  0.018;
  distance[0][1] =  0.16424;
  distance[1][1] =  0.15;

  double t1, t2, t3, QKcorr, QAcorr, lifetime;

  double QK[2][100], QA[2][100], R[2][100], x[2][100];
  double t1arr[2][100], life[2][100], errlife[2][100];
  double errR[2][100], errx[2][100], errV[2][100], errt[2][100];
  double fields[3];

  double tTheory[3][3];
  double nomRatio[2];

  int count=0;
  
  string kch[2]={"ch1", "ch4"};
  string ach[2]={"ch2", "ch5"};

  for (int iprm=0; iprm<2; iprm++){

    count=0;
    for (int irun=firstRun; irun<lastRun; irun++){

      TFile *fin = new TFile(Form("/data/PurityMonitor/Filling/Run%03d/PrM%d_lifeInfo.root", irun, iprm+1));
      TTree *metaTree = (TTree*)fin->Get("metaTree");
      metaTree->SetBranchAddress("fields[3]", fields);
      metaTree->GetEntry(0);
      cout << fields[0] << " " << fields[1] << " " << fields[2] << endl;

      x[iprm][count]  = fields[1]*1.0/fields[0] + 0.01*iprm;
      errx[iprm][count]=0;
      for (int ifield=0; ifield<3; ifield++) tTheory[iprm][ifield] = distance[iprm][ifield]/ICARUSpolynomial(fields[ifield]);  

      TTree *lifeTree = (TTree*)fin->Get("lifeTree");
      lifeTree->SetBranchAddress("t1", &t1);
      lifeTree->SetBranchAddress("t2", &t2);
      lifeTree->SetBranchAddress("t3", &t3);
      lifeTree->SetBranchAddress("QK", &QKcorr);
      lifeTree->SetBranchAddress("QA", &QAcorr);
      lifeTree->SetBranchAddress("lifetime", &lifetime);
      lifeTree->GetEntry(0);

      TGraph *gK = (TGraph*)fin->Get(Form("gfin_%s", kch[iprm].c_str()));
      QK[iprm][count] = TMath::MinElement(gK->GetN(), gK->GetY());
      TGraph *gA = (TGraph*)fin->Get(Form("gfin_%s", ach[iprm].c_str()));
      QA[iprm][count] = TMath::MaxElement(gA->GetN(), gA->GetY());

      QK[iprm][count] = -QKcorr;
      QA[iprm][count] = QAcorr;
      t1arr[iprm][count] = t1*1e6;
      errt[iprm][count]=1;
      life[iprm][count]=lifetime*1e6;
      errlife[iprm][count]=lifetime*1e6*0.1;

      cout <<  QK[iprm][count] << " " <<  QA[iprm][count] << " "  ;

      QA[iprm][count] *= TMath::Exp((t2+(t1+t3)*0.5)/tauLife);

      R[iprm][count]  = -QA[iprm][count]/QK[iprm][count];
    
      errV[iprm][count] = errQA;

      cout <<  QA[iprm][count] << " " <<  R[iprm][count]  << endl;

      errR[iprm][count] = R[iprm][count]*TMath::Sqrt( (errQK/QK[iprm][count])*(errQK/QK[iprm][count]) + (errQA/QA[iprm][count])*(errQA/QA[iprm][count])  );

      if (fields[1]/fields[0]==2) nomRatio[iprm] =  R[iprm][count];
      count++;

      delete fin;
    }

  }

  

  // renormalise to case of 2
  for (int i=0; i<count; i++){
    cout << x[1][i] << " " << R[1][i] << " " ;
    R[0][i]/=nomRatio[0];
    R[1][i]/=nomRatio[1];

    errR[0][i]/=nomRatio[0];
    errR[1][i]/=nomRatio[1];

    cout << R[1][i] << endl;
  }

  TCanvas *c = new TCanvas("c");
  c->SetGridy();

  TGraphErrors *g1 = new TGraphErrors (count, x[0], R[0], errx[0], errR[0]);
  g1->GetYaxis()->SetRangeUser(0, 1.5);
  g1->SetMarkerStyle(20);
  g1->SetTitle("Grid transparency check;E2/E1;-QA/QK (normalised to 1)");
  g1->Draw("Ap");

  TGraphErrors *g2 = new TGraphErrors (count, x[1], R[1], errx[1], errR[1]);
  g2->SetMarkerStyle(20);
  g2->SetMarkerColor(kRed);
  g2->SetLineColor(kRed);
  g2->Draw("p");

  TLegend *leg = new TLegend(0.6, 0.11, 0.89, 0.2);
  leg->AddEntry(g1, "Purity Monitor 1", "lp");
  leg->AddEntry(g2, "Purity Monitor 2", "lp");
  leg->Draw();

  c->Print(Form("%s.png", outname.c_str()));


  TGraphErrors *g1QK = new TGraphErrors (count, x[0], QK[0], errx[0], errV[0]);
  g1QK->SetMarkerStyle(20);
  TGraphErrors *g2QK = new TGraphErrors (count, x[1], QK[1], errx[1], errV[1]);
  g2QK->SetMarkerStyle(20);
  g2QK->SetMarkerColor(kRed);
  g2QK->SetLineColor(kRed);
  TGraphErrors *g1t1 = new TGraphErrors (count, x[0], t1arr[0], errx[0], errt[0]);
  g1t1->SetMarkerStyle(20);
  TGraphErrors *g2t1 = new TGraphErrors (count, x[1], t1arr[1], errx[1], errt[1]);
  g2t1->SetMarkerStyle(20);
  g2t1->SetMarkerColor(kRed);
  g2t1->SetLineColor(kRed);

  TGraphErrors *g1life = new TGraphErrors (count, x[0], life[0], errx[0], errlife[0]);
  g1life->SetMarkerStyle(20);
  TGraphErrors *g2life = new TGraphErrors (count, x[1], life[1], errx[1], errlife[1]);
  g2life->SetMarkerStyle(20);
  g2life->SetMarkerColor(kRed);
  g2life->SetLineColor(kRed);

  g1QK->GetYaxis()->SetRangeUser(10, 55);
  g1QK->SetTitle("Grid transparency check;E2/E1;QK [mV]");
  g1QK->Draw("Ap");
  g2QK->Draw("p");
  leg->Draw();
  c->Print(Form("%s_QK.png", outname.c_str()));

  
  g1t1->GetYaxis()->SetRangeUser(25,50);
  g1t1->SetTitle("Grid transparency check;E2/E1;t1 [#mus]");
  g1t1->Draw("Ap");
  g2t1->Draw("p");
  leg->Draw();
  tTheory[0][0]*=1e6;
  tTheory[1][0]*=1e6;
  double y1 =  (tTheory[0][0]-25.)/25.;
  double y2 =  (tTheory[1][0]-25.)/25.;

  cout << tTheory[0][0] << " " << tTheory[1][0] << " " << y1 << " " << y2 << endl;

  TLine *l1 = new TLine(0, y1, 1, y1);
  l1->SetLineColor(kBlue);
  l1->SetLineWidth(2);
  l1->Draw("ndc");
  TLine *l2 = new TLine(0, y2, 1, y2);
  l2->SetLineColor(kRed);
  l2->SetLineWidth(2);
  l2->Draw();
  c->Print(Form("%s_t1.png", outname.c_str()));

  g1life->GetYaxis()->SetRangeUser(0, 800);
  g1life->SetTitle("Grid transparency check;E2/E1;Lifetime [#mus]");
  g1life->Draw("Ap");
  g2life->Draw("p");
  leg->Draw();
  c->Print(Form("%s_lifetime.png", outname.c_str()));


}




Double_t ICARUSpolynomial(Double_t E){

  // transform E in kV/cm

  E*=1e-3;
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
