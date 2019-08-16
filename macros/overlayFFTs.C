void overlayFFTs(string folder="/data/PurityMonitor/Filling/", int run1=217, int run2=203, int run3=213, string name1="Standard config", string name2="SHV filter 160kHz", string name3="BNC filter 160kHz"){

  string outsuffix="Avg";

  // TFile *fin1 = new TFile(Form("%s/Run%03d/PrM2_exampleFFTs.root", folder.c_str(), run1));
  // TGraph *g1 = (TGraph*)fin1->Get("g_ch4_1_fft");
  // fin1->Close();

  // TFile *fin2 = new TFile(Form("%s/Run%03d/PrM2_exampleFFTs.root", folder.c_str(), run2));
  // TGraph *g2 = (TGraph*)fin2->Get("g_ch4_1_fft");
  // fin2->Close();

  TFile *fin1 = new TFile(Form("%s/Run%03d/PrM2_filtAvg.root", folder.c_str(), run1));
  TGraph *g1 = (TGraph*)fin1->Get("gpow_orig_ch4");
  fin1->Close();

  TFile *fin2 = new TFile(Form("%s/Run%03d/PrM2_filtAvg.root", folder.c_str(), run2));
  TGraph *g2 = (TGraph*)fin2->Get("gpow_orig_ch4");
  fin2->Close();

  TFile *fin3 = new TFile(Form("%s/Run%03d/PrM2_filtAvg.root", folder.c_str(), run3));
  TGraph *g3 = (TGraph*)fin3->Get("gpow_orig_ch4");
  fin3->Close();

  TLegend *leg = new TLegend(0.6, 0.11, 0.89, 0.25);

  cout << g1->GetN() << endl;

  TGraph *gratio2 = new TGraph(g1->GetN());
  TGraph *gratio3 = new TGraph(g1->GetN());

  for (int ip=0;ip<g1->GetN(); ip++){
    gratio2->GetX()[ip] = g1->GetX()[ip];
    gratio3->GetX()[ip] = g1->GetX()[ip];

    gratio2->GetY()[ip] = g2->GetY()[ip]/g1->GetY()[ip];
    gratio3->GetY()[ip] = g3->GetY()[ip]/g1->GetY()[ip];

  }
  



  TCanvas *c = new TCanvas("c", "", 800, 800);
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd(); 
  pad1->SetLogy();

  g1->SetTitle("Power Spectrum;Frequency [MHz];Amplitude");
  gratio2->SetTitle(";Frequency [MHz];Amplitude Ratio");

  gratio2->GetXaxis()->SetTitleSize(20);
  gratio2->GetXaxis()->SetTitleFont(43);
  gratio2->GetXaxis()->SetTitleOffset(4.);
  gratio2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  gratio2->GetXaxis()->SetLabelSize(15);

  g1->SetLineWidth(2);
  g2->SetLineWidth(2);
  g3->SetLineWidth(2);

  g2->SetLineColor(kRed);
  g3->SetLineColor(kBlue);


  gratio2->SetLineWidth(2);
  gratio3->SetLineWidth(2);
  gratio2->SetLineColor(kRed);
  gratio3->SetLineColor(kBlue);
  
  g1->Draw("Al");
  g2->Draw("l");
  g3->Draw("l");

  leg->AddEntry(g1, name1.c_str(), "l");
  leg->AddEntry(g2, name2.c_str(), "l");
  leg->AddEntry(g3, name3.c_str(), "l");
  leg->Draw();

  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();  
  gratio2->GetYaxis()->SetRangeUser(0, 1);
  gratio2->Draw("Al");
  gratio3->Draw("l");

  c->Print(Form("%s/Run%03d/FFTsComparison%s.png", folder.c_str(), run1, outsuffix.c_str()));

  pad1->cd();
  g1->GetXaxis()->SetRangeUser(0, 50);
  gratio2->GetXaxis()->SetRangeUser(0, 50);
  gratio2->GetYaxis()->SetRangeUser(0, 1);
  g1->Draw("Al");
  g2->Draw("l");
  g3->Draw("l");
  leg->Draw();
  pad2->cd();  
  gratio2->Draw("Al");
  gratio3->Draw("l");
  c->Print(Form("%s/Run%03d/FFTsComparison_upto50MHz%s.png", folder.c_str(), run1, outsuffix.c_str()));

  pad1->cd();
  g1->GetXaxis()->SetRangeUser(0, 10);
  gratio2->GetXaxis()->SetRangeUser(0, 10);
  gratio2->GetYaxis()->SetRangeUser(0, 1);
  g1->Draw("Al");
  g2->Draw("l");
  g3->Draw("l");
  leg->Draw();
  pad2->cd();  
  gratio2->Draw("Al");
  gratio3->Draw("l");
  c->Print(Form("%s/Run%03d/FFTsComparison_upto10MHz%s.png", folder.c_str(), run1, outsuffix.c_str()));

  pad1->cd();
  g1->GetXaxis()->SetRangeUser(0, 1);
  gratio2->GetXaxis()->SetRangeUser(0, 1);
  gratio2->GetYaxis()->SetRangeUser(0, 1);
  g1->Draw("Al");
  g2->Draw("l");
  g3->Draw("l");
  leg->Draw();
  pad2->cd();  
  gratio2->Draw("Al");
  gratio3->Draw("l");
  c->Print(Form("%s/Run%03d/FFTsComparison_upto1MHz%s.png", folder.c_str(), run1, outsuffix.c_str()));

  pad1->cd();
  g1->GetXaxis()->SetRangeUser(0, 0.2);
  gratio2->GetXaxis()->SetRangeUser(0, 0.2);
  gratio2->GetYaxis()->SetRangeUser(0, 1);
  g1->Draw("Al");
  g2->Draw("l");
  g3->Draw("l");
  leg->Draw();
  pad2->cd();  
  gratio2->Draw("Al");
  gratio3->Draw("l");
  c->Print(Form("%s/Run%03d/FFTsComparison_upto0.2MHz%s.png", folder.c_str(), run1, outsuffix.c_str()));

}
