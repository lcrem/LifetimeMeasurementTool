#include "UsefulFunctions.h"
#include "LifetimeConventions.h"

#include "TFile.h"
#include "TMath.h"
#include "TAxis.h"
#include "TF1.h"
#include "TAxis.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TCanvas.h"
#include <iostream>

TGraph *UsefulFunctions::justAverage(Int_t numGraphs, TGraph **grPtrPtr)
{
  //Assume they are all at same sampling rate

  // Can't correlate and average if there's only one graph.
  // So return 0
  if(numGraphs<2) return NULL;

  // TGraph *grA = grPtrPtr[0];
  TGraph *grA = (TGraph*) grPtrPtr[0]->Clone(); // Make copy of graph rather than using first graph.

  // Copy times from grA into new array.
  Int_t numPoints=grA->GetN();  
  Double_t *timeVals= grA->GetX();
  Double_t *safeTimeVals = new Double_t[numPoints];
  Double_t *sumVolts = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) 
    safeTimeVals[i]=timeVals[i];  
  

  // Loop over graph array.
  int countWaves=1;
  for(int graphNum=1;graphNum<numGraphs;graphNum++) {
    
    TGraph *grB = grPtrPtr[graphNum];
    if(grB->GetN()<numPoints)
      numPoints=grB->GetN();
 
    Double_t *aVolts = grA->GetY();
    Double_t *bVolts = grB->GetY();

    for(int ind=0;ind<numPoints;ind++) {
      sumVolts[ind]=(aVolts[ind]+bVolts[ind]);
    }
  

    TGraph *grComAB = new TGraph(numPoints,safeTimeVals,sumVolts);

    //    delete grB;
    //    if(graphNum>1)
    delete grA;
    grA=grComAB;
    countWaves++;

  }
  for(int i=0;i<numPoints;i++) {
    sumVolts[i]/=countWaves;
  }
  Double_t meanVal=TMath::Mean(50,sumVolts);
  for(int i=0;i<numPoints;i++) {
    sumVolts[i]-=meanVal;
  }
  delete grA;
  TGraph *grRet = new TGraph(numPoints,safeTimeVals,sumVolts);
  delete [] safeTimeVals;
  delete [] sumVolts;
  return grRet;
}


TGraph *UsefulFunctions::getZeroedAverage(Int_t numGraphs, TGraph **graphs){
  
  Double_t newY[1000000];
  Int_t numPoints;
  
  TGraph *graphsZeroed[1000];

  int count = 0;
  
  for (int i=0; i<numGraphs; i++){
    
    numPoints = graphs[i]->GetN();

    /* TGraph *gtemp = smoothGraph(graphs[i], 20); */
    /* Double_t minAll = TMath::Abs(TMath::MinElement(gtemp->GetN(), gtemp->GetY()));   */
    /* if (minAll < 0.1 ) continue;   */
    /* delete gtemp; */
    
    Double_t meanVal=TMath::Mean(50, graphs[i]->GetY());

    for(int ip=0;ip<numPoints;ip++) {
      newY[ip] = graphs[i]->GetY()[ip] - meanVal;
    }
    graphsZeroed[count] = new TGraph(numPoints,graphs[i]->GetX(), newY);
    count++;
  }

  
  TGraph *zeroedAverage = justAverage( count,
				       graphsZeroed );

  Double_t mean=TMath::Mean(50, zeroedAverage->GetY());

  for(int ip=0;ip<numPoints;ip++) {
    newY[ip] = zeroedAverage->GetY()[ip] - mean;
  }

  zeroedAverage = new TGraph(numPoints,zeroedAverage->GetX(), newY);


  return zeroedAverage;
}


Int_t UsefulFunctions::avgSomeGraphs(TGraph **graphs, int nmax, TGraph **g){

  // TFile *f = new TFile(filename.c_str(), "read");
  // TGraph *graphs[1000];
  int count = 0;
  int ntot = 0;

  double newy[10005];
  double newx[10005];

  TGraph *gbatch[1000];
  for (int i=0; i<1000; i++){
    // cout << i << endl;
    if(!graphs[i]) break;

    gbatch[count] = (TGraph*) graphs[i]->Clone();

    count++;
    if (count==nmax){
      
      g[ntot] = getZeroedAverage( count, gbatch );
      zeroBaseline(g[ntot]);
      /* g[ntot]  = smoothGraph(g[ntot],  5);  */
      
      count = 0;
      ntot++;
      // delete gtemp;
    }
  }

  // f->Close();

  // cout << ntot << endl;
  return ntot;
  
}


void UsefulFunctions::zeroBaseline(TGraph *g){

  double *y = g->GetY();
  Double_t meanVal=TMath::Mean(50,g->GetY());
  for(int ip=0;ip<g->GetN();ip++) {
    y[ip] -=  meanVal;
  }

}


Double_t UsefulFunctions::ICARUSpolynomial(Double_t E){

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

Double_t UsefulFunctions::getCorrectionFactor(double fittedDriftTime, double tauelec, double taulife){
  
  double laura = (1.-exp(-(fittedDriftTime)*(1./tauelec + 1./taulife)))/(fittedDriftTime*(1./tauelec + 1./taulife));
  double alan  = (exp(-fittedDriftTime/tauelec)-exp(-fittedDriftTime/taulife))/(fittedDriftTime*(1/taulife - 1/tauelec));
  
  std::cout << "Fitted drift time " << fittedDriftTime << std::endl;
  std::cout << "Tauelec " << tauelec << std::endl;
  std::cout << "Taulife " << taulife << std::endl;
  std::cout << "Correction laura " << laura << "; alan " << alan << std::endl;
  return laura;
  
}


Double_t UsefulFunctions::fittingFunction(Double_t *x, Double_t *par){
  
  double t = x[0]*1E6 - par[3];
  double sigma = par[0];
  double tau   = par[1];
  double a     = par[2];

  double y = a*exp((sigma*sigma - 2*t*tau)/(2*tau*tau))*(1 + erf((-sigma*sigma+t*tau)/(sqrt(2)*tau*sigma)));
  
  return y;
}



Int_t UsefulFunctions::calculateLifetime(TGraph *gK, TGraph *gA, int whichPrM, double tTheory[3], double lifetime[2], bool saveCanvas){
      
  int nK      = gK->GetN();
  double *xK  = gK->GetX();
  double *yK  = gK->GetY();
      
  Double_t tauelecK, tauelecA, gainAoverK;

  switch(whichPrM){
  case 0:
    tauelecK   = PrM1preamp[0];
    tauelecA   = PrM1preamp[1];
    gainAoverK = PrM1preamp[2];
    break;
  case 1:
    tauelecK   = PrM2preamp[0];
    tauelecA   = PrM2preamp[1];
    gainAoverK = PrM2preamp[2];
    break;
  case 2:
    tauelecK   = PrM1preamp[1];
    tauelecA   = PrM1preamp[0];
    gainAoverK = 1./PrM1preamp[2];
    break;
  case 3:
    tauelecK   = PrM2preamp[1];
    tauelecA   = PrM2preamp[0];
    gainAoverK = 1./PrM2preamp[2];
    break;
  }
  
  double fittedK;
  double fittedKtime;
  double fittedKstartTime = 0;
      
  double fittedA;
  double fittedAtime;
  double fittedAstartTime = 0;
      
  int loc;
  double peak=0;
        
  std::cout << " Cathode field " << tTheory[0] << std::endl;
  
  for (int ip=nK-2; ip>0; ip--){
    if (xK[ip]<0) break;
    if (yK[ip]<peak && (xK[ip]>(tTheory[0]-100.E-6)) && (xK[ip]<(tTheory[0]+100.E-6)) ){
      peak = yK[ip];
      loc = ip;
    }
  }
      
  fittedKtime = xK[loc];
  fittedK = yK[loc];
  std::cout << "K time and K " << fittedKtime << " " << fittedK << std::endl;
  fittedKstartTime = 0.;
      
  // TF1 *funcK = new TF1("funcK",fittingFunction,0.,0.9E-3,4);
  // funcK->SetParameters(5,  tauelecK*1E6, -0.002, (tTheory[0])*1e6);
  // // funcK->FixParameter(0, 5);
  // // funcK->FixParameter(1, 43);
  // // funcK->FixParameter(2, 0.002);
  // // funcK->FixParameter(3, (tTheory[0]+tTheory[1])*1e6);
  // // funcK->SetParLimits(0, 1, 50);
  // funcK->SetParLimits(1, tauelecK*1E6*0.95, tauelecK*1E6*1.05);
  // // funcK->SetParLimits(2, 0, 0.1);
  // // funcK->SetParLimits(3, 10, 500);
  // funcK->SetParName(0, "#sigma");
  // funcK->SetParName(1, "#tau_{D}");
  // funcK->SetParName(2, "a");
  // funcK->SetParName(3, "t_0");
  // gK->Fit("funcK", "R", "", 0, tTheory[0]+tTheory[1]-0.1E-3);

  
  
    
  double tempx, tempy;
  // for (int ip=0; ip<nK; ip++){
  //   tempx = xK[ip];
  //   tempy = funcK->Eval(xK[ip]);
  //   if (tempy<0.01*fittedK){
  //     std::cout << ip << " " << tempx << " " << tempy << " " << std::endl;
  //     fittedKstartTime = tempx;
  //     break;
  //   }
  // }
  // if (fittedKstartTime<0) fittedKstartTime=0;
  // fittedKtime = funcK->GetMinimumX();
  // fittedK = funcK->GetMinimum();
  

  
  // cout << " Anode field " << fields[2] << endl;
  int nA      = gA->GetN();
  double *xA  = gA->GetX();
  double *yA  = gA->GetY();


  loc = TMath::LocMax(nA,yA);
  peak = -99999.;
  for (int ip=nA-2; ip>0; ip--){
    if (xA[ip]<tTheory[0]) break;
    if (yA[ip]>peak){
      peak = yA[ip];
      loc = ip;
    }
  }
	
      
  fittedAtime = xA[loc];
  fittedA = yA[loc];

  for (int ip=loc; ip>0; ip--){
    tempx = xA[ip];
    tempy = yA[ip];
    if (xA[ip]<tTheory[0]) break;
    if (tempy<0.01*fittedA){
      // cout << ip << " " << tempx << " " << tempy << " " << endl;
      fittedAstartTime = tempx;
      break;
    }
  }


  bool fitAnode = false;

  if (fitAnode){
      
    TF1 *funcA = new TF1("funcA",fittingFunction,0.,0.9E-3,4);
    funcA->SetParameters(5, tauelecA, 0.002, (tTheory[0]+tTheory[1])*1e6);
    // funcA->FixParameter(0, 5);
    // funcA->FixParameter(1, 43);
    // funcA->FixParameter(2, 0.002);
    // funcA->FixParameter(3, (tTheory[0]+tTheory[1])*1e6);
    // funcA->SetParLimits(0, 1, 50);
    funcA->SetParLimits(1, tauelecA*1E6*0.95, tauelecA*1E6*1.05);
    // funcA->SetParLimits(2, 0, 0.1);
    // funcA->SetParLimits(3, 10, 500);
    funcA->SetParName(0, "#sigma");
    funcA->SetParName(1, "#tau_{D}");
    funcA->SetParName(2, "a");
    funcA->SetParName(3, "t_0");
    gA->Fit("funcA", "R", "", tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1]);
    
    
    
    // TF1 *funcA = new TF1("funcA",signalFunction,tTheory[0]+tTheory[1]-0.1E-3, xK[nK-1],6);
    // funcA->SetParName(0, "Q [pC]");
    // funcA->SetParName(1, "R [MOhm]");
    // funcA->SetParName(2, "C [pF]");
    // funcA->SetParName(3, "taulife [us]");
    // funcA->SetParName(4, "tdrift [us]");
    // funcA->SetParName(5, "t0 [us]");
    // funcA->SetParameters(5, R, C, 1000, icarusTime*1E6, (tTheory[0]+tTheory[1])*1E6);
    // // funcA->FixParameter(0, 5);
    // // funcA->FixParameter(1, R);
    // // funcA->FixParameter(2, C);
    // //	  funcA->FixParameter(3, 1000);
    // // funcA->FixParameter(4, icarusTime*1E6);
    // // funcA->FixParameter(5,  (tTheory[0]+tTheory[1])*1E6);
    // gdiff2->Fit("funcA", "R", "", tTheory[0]+tTheory[1]-0.1E-3, xK[nK-1]);
    
    
    for (int ip=0; ip<nA; ip++){
      tempx = xA[ip];
      tempy = funcA->Eval(xA[ip]);
      if (tempy>0.01*fittedA){
	// cout << ip << " " << tempx << " " << tempy << " " << endl;
	fittedAstartTime = tempx;
	break;
      }
    }
    fittedA     = funcA->GetMaximum();
    fittedAtime = funcA->GetMaximumX();
  }
      
          
      
  std::cout << " Expected vs measured time at cathode " << tTheory[0] << " " << fittedKtime - fittedKstartTime << std::endl;
  std::cout << " Expected vs measured time at anode " << tTheory[2] << " " << fittedAtime - fittedAstartTime << std::endl;
      
    
  // //  TF1 *funcAll = new TF1("funcAll", [&](double*x, double *p){ return (funcA->Eval(x[0])-funcK->Eval(x[0])); }, 0,xA[nA-1],8);
  // TF1 *funcAll = new TF1("funcAll", fittingCathodeAndAnode, 0,xA[nA-1],8);
  // funcAll->SetParameters(funcA->GetParameter(0), funcA->GetParameter(1), funcA->GetParameter(2), funcA->GetParameter(3), funcK->GetParameter(0), funcK->GetParameter(1), funcK->GetParameter(2), funcK->GetParameter(3));
  // funcAll->SetParLimits(1, tauelecA*1E6*0.95, tauelecA*1E6*1.05);
  // funcAll->SetParLimits(5, tauelecK*1E6*0.95, tauelecK*1E6*1.05);
  // TGraph *gAll =  subtractGraphs(gdiff[0], gdiff[1]);
  // gAll->Draw("Al");
  // funcAll->Draw("same");
  // gAll->Fit("funcAll", "R", "", 0, xA[nA-1]);


  // cout << "BEFORE " << endl;
  // cout << fittedK << " " << fittedKtime << " " << fittedA << " " << fittedAtime << endl;

  // fittedK     = funcAll->GetMaximum(0, tTheory[0]+0.1E-6)*(-1);
  // fittedKtime = funcAll->GetMaximumX(0, tTheory[0]+0.1E-6);

  // fittedA     = funcAll->GetMaximum(tTheory[0]+tTheory[1]/2, xA[nA-1]);
  // fittedAtime = funcAll->GetMaximumX(tTheory[0]+tTheory[1]/2, xA[nA-1]);

  // cout << "AFTER " << endl;
  // cout << fittedK << " " << fittedKtime << " " << fittedA << " " << fittedAtime << endl;
  

  
  
  //  return;
  double tK, tGK, tGA, tA, t1, t2, t3;
  double QA, QK, R;
  tK  = fittedKstartTime;
  tGK = fittedKtime;
  tGA  = fittedAstartTime;
  tA = fittedAtime;
            
  QA = fittedA;
  QK = fittedK;
      
  t1 = tGK - tK;
  t2 = tGA - tGK;
  t3 = tA - tGA;
  R = 0;
  
  
  if (tK<tGK && tGK<tGA && tGA<tA){
    // All good
  } else {
    std::cout << "THESE NUMBERS DON'T MAKE SENSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << tK << " " << tGK << " " << tGA << " " << tA << std::endl;
    lifetime[0]=lifetime[1];
    return -1;
  }
  
  
  
  //  QK *= gainAoverK;
  
  double taulife = 0.001;
  double Kcorrection = 0.;
  double Acorrection = 0.;
  double newQK = QK;
  double newQA = QA;
  
  if (R<1){
    
    int count = 0;
    while(count < 20){
      
      //cout << "Catodo " << endl;
      Kcorrection = getCorrectionFactor(t1, tauelecK, taulife);
      //cout << "Anodo " << endl;
      Acorrection = getCorrectionFactor(t3, tauelecA, taulife);
      
      newQK = QK/Kcorrection;
      newQA = QA/Acorrection;
      
      R =  TMath::Abs(newQA/newQK);
      
      lifetime[0] = (1/TMath::Abs(TMath::Log(R)))*(t2 + 0.5*(t1+t3));
      
      std::cout <<" This is my lifetime " <<  lifetime[0] << std::endl;
      
      if (TMath::Abs(lifetime[0]-taulife)<0.001e-6) break;
      taulife=lifetime[0];
      count++;
    }
    std::cout << "Number of iterations : " << count << std::endl;
    
    TF1 f ("lifetime function", "([1]/[3])*(TMath::SinH([3]/(2*x))/TMath::SinH([1]/(2*x)))*TMath::Exp(-([2]+0.5*([1]+[3]))/x) - [0]", -0.1, 0.1);
    f.SetParameters(R, t1, t2, t3);
	
    ROOT::Math::WrappedTF1 wf1(f);
	
    // Create the Integrator
    ROOT::Math::BrentRootFinder brf;
	
    // Set parameters of the method
    brf.SetFunction( wf1, -0.1, 0.1 );
    brf.Solve();
	
    lifetime[1] =  TMath::Abs(brf.Root() ) ;
	

  } else {
    lifetime[0] = lifetime[1] = 0;
  }

  
  printf("tK     : %12.4e \n",  tK  );
  printf("tGK    : %12.4e \n",  tGK );
  printf("tGA    : %12.4e \n",  tGA );
  printf("tA     : %12.4e \n",  tA  );
  printf("t1     : %12.4e \n",  t1  );
  printf("t2     : %12.4e \n",  t2  );
  printf("t3     : %12.4e \n",  t3  );
  printf("QA     : %12.4e \n",  QA  );
  printf("QK     : %12.4e \n",  QK  );
  printf("QAcorr : %12.4e \n",  newQA  );
  printf("QKcorr : %12.4e \n",  newQK  );
  printf("R      : %12.4e \n",  R   );
  printf("lifetime : %12.4e \n",  lifetime[0]  );
  printf("lifetime2: %12.4e \n",  lifetime[1] );


  if (saveCanvas){
    TCanvas *c = new TCanvas("c");
    gK->GetYaxis()->SetRangeUser(QK*1.2, QA*1.2);
    gK->Draw("Al");
    gA->Draw("l");

    c->Print("Lastlifetime.png");
    c->Print("Lastlifetime.root");
  }
  
  return 1;
}

TGraph *UsefulFunctions::translateGraph(const TGraph *grWave, Double_t deltaT){
  Int_t N = grWave->GetN();
  Double_t *X=grWave->GetX();
  Double_t *Y=grWave->GetY();
  Double_t *newX = new Double_t[N];
  for(int i=0;i<N;i++)
    newX[i]=X[i]+deltaT;
  TGraph *grOut = new TGraph(N,newX,Y);
  delete [] newX;
  return grOut;
}

Int_t UsefulFunctions::getSmoothingNumber(double deltat, double tdrift){

  int nsigma = 10;
  int nsmooth = tdrift/(nsigma*2*deltat);
  if (nsmooth>20) nsmooth=20;
  return  nsmooth;

}


Double_t UsefulFunctions::greenFunction(Double_t *x, Double_t *par){

  // x[0] is in seconds, so we need to convert par[3] from musec to sec
  double t = x[0];
  double GQ = par[0];
  double tau = par[1]*1e-6;
  double rise = par[2]*1e-6;
  double shift = par[3]*1e-6;

  double heaviside;
  if (t<0) heaviside=0;
  else heaviside = 1;

  //double y = (-GQ)*exp(-t/(tau))*heaviside;

  double y;
  double tmp1 = GQ*(-1-erf((0+shift)/rise));
  double tmp2 = -GQ*exp(-0/tau);
  if (t>0) y = -GQ*exp(-t/tau);
  else y = GQ*(-1-erf((t+shift)/rise))*tmp2/tmp1;

  return y;
}
