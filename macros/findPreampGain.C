double getMinimum(TGraph *g1);
double getMaximum(TGraph *g1);

void findPreampGain_each(bool isCathode, string filter, string PrM);

void findPreampGain(){

  findPreampGain_each(true,  "wFilters",  "1");
  findPreampGain_each(false, "wFilters",  "1");
  findPreampGain_each(true,  "noFilters", "1");
  findPreampGain_each(false, "noFilters", "1");
  findPreampGain_each(true,  "wFilters",  "2");
  findPreampGain_each(false, "wFilters",  "2");
  findPreampGain_each(true,  "noFilters", "2");
  findPreampGain_each(false, "noFilters", "2");

}


void findPreampGain_each(bool isCathode, string filter, string PrM){

  string ch;
  string chNice;
  string preamp1;
  string preamp2;
  string postname1, postname2;

  if (isCathode){
    ch+="ch4";
    chNice+="Cathode";
    if (PrM=="1"){
      preamp1 += "A";
      preamp2 += "B";
      postname1 += "aCathode_bAnode";
      postname2 += "bCathode_aAnode";
    } else if (PrM=="2"){
      preamp1 += "C";
      preamp2 += "D";
      postname1 += "cCathode_dAnode";
      postname2 += "dCathode_cAnode";
    }
  } else {
    ch+="ch3";
    chNice+="Anode";
    if (PrM=="1"){
      preamp1 += "B";
      preamp2 += "A";
      postname1 += "aCathode_bAnode";
      postname2 += "bCathode_aAnode";
    } else if (PrM=="2"){
      preamp1 += "D";
      preamp2 += "C";
      postname1 += "cCathode_dAnode";
      postname2 += "dCathode_cAnode";
    }
  }


  string name1 = "/unix/dune/purity/CERN/2019/Liquid/PrM"+PrM+"/Day4/AllFibres_"+filter+"_"+postname1+"_newResistors/Field_25.50.100Vcm."+ch+".traces_averages.root";
  string name2 = "/unix/dune/purity/CERN/2019/Liquid/PrM"+PrM+"/Day4/AllFibres_"+filter+"_"+postname2+"_newResistors/Field_25.50.100Vcm."+ch+".traces_averages.root";

  string whichprm = PrM+filter;

  TFile *f1 = new TFile(name1.c_str(), "read");
  TFile *f2 = new TFile(name2.c_str(), "read");

  TGraph *g1 = (TGraph*)f1->Get("justAvg");
  TGraph *g2 = (TGraph*)f2->Get("justAvg");

  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.6, 0.31, 0.89, 0.44);
  leg->AddEntry(g1, Form("PreAmp %s", preamp1.c_str() ), "l");
  leg->AddEntry(g2, Form("PreAmp %s", preamp2.c_str() ), "l");

  TCanvas *c = new TCanvas("c");
  g1->SetTitle(Form("PreAmps for PrM %s, %s signal;Time [s];Amplitude [V]", whichprm.c_str(), chNice.c_str()));
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
    krms1 += (somek1[i]-kmean1)*(somek1[i]-kmean1);
    krms2 += (somek2[i]-kmean2)*(somek2[i]-kmean2);
  }
  krms1 = TMath::Sqrt(krms1)/ngraphs;
  krms2 = TMath::Sqrt(krms2)/ngraphs;
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

  c->Print(Form("plots/GainCalibrationForPrM%s_from%s.png", whichprm.c_str(), chNice.c_str()));
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
