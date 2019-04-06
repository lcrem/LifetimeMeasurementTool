#include "UsefulFunctions.h"
#include "LifetimeConventions.h"

#include "TFile.h"
#include "TF1.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace std;

double xmin = 0.05E-3;
double xmax = 0.015E-3;
double ymin = -0.8;
double ymax = +0.8;


int main(int argc, char *argv[]){

  
  string basename, fieldname, preamp;

  if((argc!=4)){
    std::cerr << "Usage : " << argv[0] << " [basename] [fieldname] [preamp]" << std::endl;
    return 1;
  } else {
    basename += argv[1];
    fieldname += argv[2];
    preamp += argv[3];
  }
  
  string outfile    = basename + fieldname + "_signals.root";
  string outcanvas    = basename + fieldname + "_tauelec";
  
  cout << "Output file is " << outfile << endl;
  
  string chname[2] = {"ch3", "ch4"};
  string chnamenice[2] = {"anode", "Cathode"};
  

  double timedelay=0;
  
  string stimedelay = basename + fieldname + ".ch1.traces_averages.root";
  TFile *ftimedelay = new TFile (stimedelay.c_str(), "read");
  TGraph *gtimedelay = (TGraph*)ftimedelay->Get("justAvg");
  double *xTimeDelay = gtimedelay->GetX();
  double *yTimeDelay = gtimedelay->GetY();
  double baseY = yTimeDelay[0];
  double step = TMath::MaxElement(gtimedelay->GetN(), yTimeDelay)*0.2;
  double timeStep = xTimeDelay[1]-xTimeDelay[0];
  for (int ip=0; ip<gtimedelay->GetN(); ip++){
    if (yTimeDelay[ip]>(step)){
      timedelay = xTimeDelay[ip];
      break;
    }
  }
  ftimedelay->Close();
  cout << "The time delay is " << timedelay << endl;
  
  double newy[20000];
  double newx[20000];
  
  TFile *out = new TFile(outfile.c_str(), "recreate");
  
  
  double finalNumbers[2][3]; // [0 anode, 1 cathode] [0 amplitude, 1 start time, 2 peak time]
  
  for (int ich=1; ich<2; ich++){
    
    string f1 = basename + fieldname + "." + chname[ich] +".traces_averages.root";
    
    TFile *file1 = new TFile(f1.c_str(), "read");
    
    TGraph *g1 = (TGraph*)file1->Get("justAvg");
    g1->SetName("g1");
    
    file1->Close();
    
    g1 = UsefulFunctions::translateGraph(g1, -timedelay);
        
    int N = g1->GetN();
    double *x = g1->GetX();
    double *y1 = g1->GetY();
    for (int i=0; i<N; i++){
      newy[i] = y1[i];
      newx[i] = x[i];
    }
    
    TGraph *gdiff = new TGraph(N, newx, newy);
      
    TCanvas *c = new TCanvas("c");
    c->Divide(1,2);
    c->cd(1);
    g1->SetLineColor(kBlue);
    g1->GetXaxis()->SetRangeUser(xmin, xmax);
    g1->GetYaxis()->SetRangeUser(ymin, ymax);
    g1->Draw("Al");
    c->cd(2);
    gdiff->GetXaxis()->SetRangeUser(xmin, xmax);
    gdiff->GetYaxis()->SetRangeUser(ymin, ymax);
    gdiff->Draw("Al");
    
    
    out->cd();
    c->Write((chnamenice[ich]+"_canvas").c_str());
    gdiff->Write((chnamenice[ich]).c_str());
    
    UsefulFunctions::zeroBaseline(gdiff);
    // TGraph *gdiff2 = UsefulFunctions::smoothGraph(gdiff, 10);
    
    // gdiff2->Write((chnamenice[ich]+"_smoothed").c_str());
    
    gStyle->SetOptFit(1);
    gStyle->SetStatY(0.4);
    gStyle->SetStatX(0.9);
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.2);
    
    TCanvas *ctemp = new TCanvas ("ctemp");
    gdiff->SetTitle(Form("%s;Time (s);Amplitude (V)", preamp.c_str()));
    gdiff->Draw("Al");
    // TF1 *func = new TF1("func",greenFunction2,-0.1E-3,0.7E-3,4);
        
    // func->SetParameters(0.01, 0.1, 500, 6);
    // func->SetParName(0, "Q (pC)");
    // //	func->SetParLimits(0, 1e-6, 2e-6);
    // func->SetParName(1, "C (pF)");
    // // func->FixParameter(1, 0.1);
    // func->SetParLimits(1, 0.1, 0.2);
    // func->SetParName(2, "R (M#Omega)");
    // //func->SetParLimits(2, 490, 510);
    // func->SetParName(3, "t_0 (#mus)");
        
    TF1 *func = new TF1("func",UsefulFunctions::greenFunction,-0.30E-3,1.2E-3,4);
        
    //par[0] = Q*G, par[1]=tau_el
    func->SetParameters(2.6, 43, 20,3);
    func->SetParName(0, "Gain x Q ");
    //func->SetParLimits(0, 2, 3);
        
    func->SetParName(1, "Tau (#mus)");
    //func->SetParLimits(1, 40, 45);

    func->SetParName(2, "Rise (#mus)");
    func->SetParName(3, "Shift (#mus)");
        
    //func->SetParName(2, "t_0 (#mus)");
    //func->SetParLimits(2, -10, 10);
        
    gdiff->Fit("func","R");
    func->Draw("same");
        
    double tauelec = func->GetParameter(1)*func->GetParameter(2);
        
    ctemp->Update();
    TLatex *text = new TLatex();
    //text->DrawTextNDC(0.5, 0.5, Form("#tau_{elec} = %3.2f #mus", tauelec));
        
    ctemp->Print((outcanvas+".png").c_str());
    ctemp->Print((outcanvas+".pdf").c_str());
    ctemp->Print((outcanvas+".C").c_str());
        
  }
      
       

  return 0;

}
