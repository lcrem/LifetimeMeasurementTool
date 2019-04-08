void findPreampGain(string preamp1="C", string preamp2="D"){

  string name1 = "/unix/dune/purity/CERN/2019/Pulser/PreAmp"+preamp1+"/Pulser.ch3.traces_averages.root";
  string name2 = "/unix/dune/purity/CERN/2019/Pulser/PreAmp"+preamp2+"/Pulser.ch3.traces_averages.root";

  TFile *f1 = new TFile(name1.c_str(), "read");
  TFile *f2 = new TFile(name2.c_str(), "read");

  TGraph *g1 = (TGraph*)f1->Get("justAvg");
  TGraph *g2 = (TGraph*)f2->Get("justAvg");

  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.6, 0.11, 0.89, 0.44);
  leg->AddEntry(g1, Form("PreAmp %s", preamp1.c_str()), "l");
  leg->AddEntry(g2, Form("PreAmp %s", preamp2.c_str()), "l");

  TCanvas *c = new TCanvas("c", "", 650, 800);
  c->Divide(1,2);
  c->cd(1);
  g1->SetTitle("PreAmps;Time [s];Amplitude [V]");
  g1->Draw("Al");
  g2->Draw("l");
  leg->Draw();

  c->cd(2);
  
  double x[10000];
  double y[10000];
  int n = g1->GetN();
  for (int i=0; i<n; i++){
    x[i] = g1->GetX()[i];
    y[i] = g1->GetY()[i]/g2->GetY()[i];
  }
  TGraph *gratio = new TGraph(n,x,y);
  gratio->SetTitle(Form("Preamp %s/%s ratio;Time [s];Ratio", preamp1.c_str(), preamp2.c_str()));
  gratio->Draw("Al");

  c->Print(Form("PreampComparison_%s_%s.png", preamp1.c_str(), preamp2.c_str()));
  
  double integral1 = g1->Integral();
  double integral2 = g2->Integral();
  double peak1 = TMath::Abs(TMath::MinElement(g1->GetN(), g1->GetY()));
  double peak2 = TMath::Abs(TMath::MinElement(g2->GetN(), g2->GetY()));
  cout << "Integrals, peaks and ratios " << endl;
  cout << "Preamp " << preamp1 << " " << integral1 << " " << peak1 << endl;
  cout << "Preamp " << preamp2 << " " << integral2 << " " << peak2 << endl;
  cout << "Ratio " << integral1/integral2 << " " << peak1/peak2 << endl;
}
