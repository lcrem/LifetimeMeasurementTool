void subtractBaseline(){

  TFile *outfile = new TFile("/data/PurityMonitor/Filling/Run016/run016_prm2_subtracted.root", "recreate");

  TFile *f1 = new TFile("/data/PurityMonitor/Filling/Run016/PrM2_filtAvg.root", "read");
  TGraph *g1 = (TGraph*)f1->Get("gfil_ch4");
  g1->SetTitle("PrM 2 at 50.100.200 V/cm");
  g1->GetYaxis()->SetTitle("[mV]");
  g1->GetXaxis()->SetTitle("[s]");  

  leg = new TLegend(0.522923, 0.177215, 0.879656, 0.377637);
  leg->AddEntry(g1, "HV on", "l");

  TGraph *g2 = (TGraph*)f1->Get("gfil_ch5");

  TFile *f2 = new TFile("/data/PurityMonitor/Filling/Run017/PrM2_filtAvg.root", "read");
  TGraph *g3 = (TGraph*)f2->Get("gfil_ch4");
  g3->SetLineColor(2);  
  leg->AddEntry(g3, "HV off", "l");

  TGraph *g4 = (TGraph*)f2->Get("gfil_ch5");
  g4->SetLineColor(2);

  TGraph *gSub_cathode = new TGraph(g1->GetN());
  TGraph *gSub_anode = new TGraph(g3->GetN());

  for (int i=0;i<g1->GetN(); i++){

    gSub_cathode->SetPoint(i, g1->GetX()[i], g1->GetY()[i]-g3->GetY()[i]);
    gSub_anode->SetPoint(i, g2->GetX()[i], g2->GetY()[i]-g4->GetY()[i]); 
 }

  gSub_cathode->SetLineColor(4);
  gSub_anode->SetLineColor(4);

  leg->AddEntry(gSub_cathode, "Subtracted", "l");

  c1 = new TCanvas();

  g1->Draw();
  g2->Draw("same");
  g3->Draw("same");
  g4->Draw("same");
  gSub_cathode->Draw("same");
  gSub_anode->Draw("same");
  leg->SetLineColor(10);
  leg->Draw("same");

  c1->Update();
  c1->Print("plots/PrM2_background_subtracted.png");

  outfile->cd();
  gSub_cathode->Write();
  outfile->Close();
}
