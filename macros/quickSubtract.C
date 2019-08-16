void quickSubtract(string fname="/data/PurityMonitor/ColdGas//Run035/PrM2_filtAvg.root"){


  TFile *fin = new TFile(fname.c_str(),"read");

  TGraph *g1 = (TGraph*)fin->Get("gfil_ch4");
  TGraph *g2 = (TGraph*)fin->Get("gfil_ch5");

  TGraph *gSub = new TGraph(g1->GetN());

  for (int i=0;i<g1->GetN(); i++){

    gSub->SetPoint(i, g1->GetX()[i], g1->GetY()[i]-g2->GetY()[i]);
  }

  gSub->Draw();

}
