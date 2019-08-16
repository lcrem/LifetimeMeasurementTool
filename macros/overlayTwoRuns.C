void overlayTwoRuns(string baseline="/data/PurityMonitor/Filling/", int whichPrM=2, int run1=241, int run2=243, int run3=243, string title="Field 50.100.200Vcm", string leg1="July 29th", string leg2="July 30th", string leg3="Ignore"){

  string chname[2];

  if (whichPrM==1){
    chname[0] += "ch1";
    chname[1] += "ch2";
  } else if (whichPrM==2) {
    chname[0] += "ch4";
    chname[1] += "ch5";
  }



  TGraph *gK[3];
  TGraph *gA[3];

  TFile *f1 = new TFile(Form("%s/Run%03d/PrM%i_lifeInfo.root", baseline.c_str(), run1, whichPrM), "read");
  gK[0] = (TGraph*)f1->Get(Form("gfin_%s", chname[0].c_str()));
  gA[0] = (TGraph*)f1->Get(Form("gfin_%s", chname[1].c_str()));

  TFile *f2 = new TFile(Form("%s/Run%03d/PrM%i_lifeInfo.root", baseline.c_str(), run2, whichPrM), "read");
  gK[1] = (TGraph*)f2->Get(Form("gfin_%s", chname[0].c_str()));
  gA[1] = (TGraph*)f2->Get(Form("gfin_%s", chname[1].c_str()));

  TFile *f3 = new TFile(Form("%s/Run%03d/PrM%i_filtAvg.root", baseline.c_str(), run3, whichPrM), "read");
  gK[2] = (TGraph*)f3->Get(Form("gfil_%s", chname[0].c_str()));
  gA[2] = (TGraph*)f3->Get(Form("gfil_%s", chname[1].c_str()));

  gK[0]->SetTitle(Form("PrM%d: %s, Cathode;Time [s];Amplitude [mV]", whichPrM, title.c_str()));
  gK[1]->SetLineColor(kRed);
  gK[2]->SetLineColor(kBlue);
  gK[0]->SetLineWidth(2);
  gK[1]->SetLineWidth(2);
  gK[2]->SetLineWidth(2);

  double minK[3], maxK[3];
  double min=0;
  double max=0;
  for (int i=0; i<3; i++){
    minK[i] = TMath::MinElement(gK[i]->GetN(), gK[i]->GetY());
    maxK[i] = TMath::MaxElement(gK[i]->GetN(), gK[i]->GetY());
    if (minK[i]<min) min = minK[i];
    if (maxK[i]>max) max = maxK[i];
  }

  gA[0]->SetTitle(Form("PrM%d: %s, Anode;Time [s];Amplitude [mV]", whichPrM, title.c_str()));
  gA[1]->SetLineColor(kRed);
  gA[2]->SetLineColor(kBlue);
  gA[0]->SetLineWidth(2);
  gA[1]->SetLineWidth(2);
  gA[2]->SetLineWidth(2);

  TCanvas *c = new TCanvas("c");

  gK[0]->GetYaxis()->SetRangeUser(min*1.1, max*1.1);
  gK[0]->Draw("Al");
  gK[1]->Draw("l");
  // gK[2]->Draw("l");

  TLegend *leg = new TLegend(0.6, 0.75, 0.89, 0.89);
  leg->AddEntry(gK[0], leg1.c_str(), "l");
  leg->AddEntry(gK[1], leg2.c_str(), "l");
  // leg->AddEntry(gK[2], leg3.c_str(), "l");
  leg->Draw();


  c->Print(Form("%s/Run%03d/PrM%i_CompareCathodeToRun%03d.png", baseline.c_str(), run1, whichPrM, run2));

  TCanvas *c2 = new TCanvas("c2");
  gA[0]->Draw("Al");
  gA[1]->Draw("l");
  // gA[2]->Draw("l");
  leg->Draw();
  c2->Print(Form("%s/Run%03d/PrM%i_CompareAnodeToRun%03d.png", baseline.c_str(), run1, whichPrM, run2));

}
