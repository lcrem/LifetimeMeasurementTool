#include "UsefulFunctions.h"
#include "LifetimeConventions.h"

#include "TFile.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TAxis.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TFitter.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TCanvas.h"
#include <iostream>

TGraph *UsefulFunctions::justAverageFromFile(TFile *fin, int numGraphs, int ch )
{
  //Assume they are all at same sampling rate

  // Can't correlate and average if there's only one graph.
  // So return 0
  if(numGraphs<2) return NULL;

  // TGraph *grA = grPtrPtr[0];
  TGraph *grA = (TGraph*) fin->Get(Form("g_ch%d_1", ch)); // Make copy of graph rather than using first graph.

  // Copy times from grA into new array.
  Int_t numPoints=grA->GetN();  
  Double_t *timeVals= grA->GetX();
  Double_t *safeTimeVals = new Double_t[numPoints];
  Double_t *sumVolts = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) 
    safeTimeVals[i]=timeVals[i];  
  
  // Loop over graph array.
  int countWaves=1;
  for(int graphNum=2;graphNum<=numGraphs;graphNum++) {
    
    TGraph *grB = (TGraph*)fin->Get(Form("g_ch%d_%d", ch, graphNum));
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
    delete grB;
    grA=grComAB;
    countWaves++;
    if (countWaves%10==0) std::cout << "\r" << countWaves*100./numGraphs << " %" << std::endl;
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
  
  // Copy times from grA into new array.                                                                                                     
  Int_t numPoints=graphs[0]->GetN();
  Double_t *newY = new Double_t[numPoints];
  
  TGraph *graphsZeroed[1000];

  int count = 0;
  
  for (int i=0; i<numGraphs; i++){

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
    zeroedAverage->GetY()[ip] -= mean;
  }

  for (int i=0; i<count; i++){
    delete graphsZeroed[i];
  }
 
  delete [] newY;

  return zeroedAverage;
}


Int_t UsefulFunctions::avgSomeGraphs(TGraph **graphs, int nmax, TGraph **g){

  int count = 0;
  int ntot = 0;
  int ngraphs = sizeof(graphs)/sizeof(graphs[0]);
  std::cout << "NUM GRAPHS ARE " << ngraphs << std::endl;
  TGraph *gbatch[1000];
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
  for (int i=0; i<ngraphs; i++){
    // cout << i << endl;
    if(!graphs[i]) break;
    std::cout << __FUNCTION__ << " " << __LINE__ << " on avg " << i << std::endl;
    gbatch[count] = (TGraph*) graphs[i]->Clone();
    std::cout << __FUNCTION__ << " " << __LINE__ << std::endl;
    count++;
    if (count==nmax){
      
      g[ntot] = getZeroedAverage( count, gbatch );
      zeroBaseline(g[ntot]);
      /* g[ntot]  = smoothGraph(g[ntot],  5);  */

      for (int ii=0; ii<count; ii++) delete gbatch[ii];
      
      count = 0;
      ntot++;
      // delete gtemp;
    }
  }

  // f->Close();

  // cout << ntot << endl;
  return ntot;
  
}


Int_t UsefulFunctions::avgSomeGraphs(TFile *fin, int nmax, TGraph **g, int ch){

  int count = 0;
  int ntot = 0;
  //  int ngraphs = sizeof(g)/sizeof(g[0]);

  TGraph *gbatch[100];
  for (int i=0; i<1000; i++){

    //    if (!fin->GetListOfKeys()->Contains(fin->Get(Form("g_ch%d_%d", ch, i+1)))) continue;

    gbatch[count] = (TGraph*) fin->Get(Form("g_ch%d_%d", ch, i+1));
    if(!gbatch[count]) break;

    count++;
    if (count==nmax){
      
     g[ntot] = getZeroedAverage( count, gbatch );
      zeroBaseline(g[ntot]);
      /* g[ntot]  = smoothGraph(g[ntot],  5);  */
      
      for (int ii=0; ii<count; ii++) delete gbatch[ii];

      count = 0;
      ntot++;
      std::cout << "Done " << ntot << std::endl;
      // delete gtemp;
    }
  }

  // f->Close();

  // cout << ntot << endl;
  return ntot;
  
}


void UsefulFunctions::zeroBaseline(TGraph *g){

  Double_t meanVal=TMath::Mean(g->GetN()*0.02,g->GetY());
  //  std::cout << "The mean is " << meanVal << std::endl;
  for(int ip=0;ip<g->GetN();ip++) {
    g->GetY()[ip] -=  meanVal;
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
  
  // std::cout << "Fitted drift time " << fittedDriftTime << std::endl;
  // std::cout << "Tauelec " << tauelec << std::endl;
  // std::cout << "Taulife " << taulife << std::endl;
  // std::cout << "Correction laura " << laura << "; alan " << alan << std::endl;
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



Int_t UsefulFunctions::calculateLifetime(TGraph *gK, TGraph *gA, int whichPrM, double tTheory[3], double lifetime[20], double lifeErrors[20], bool saveCanvas){
      
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
  
  double fittedK, fittedA;
  double tK, tGK, tGA, tA, t1, t2, t3;
  double errtK, errtGK, errtGA, errtA, errR;
  tK=tGK=tGA=tA=0.;
  double QA, QK, R, newQA, newQK;
   
  double errt1, errt2, errt3;
  double errQK = 2.;
  double errQA = 2.;
  double errNewQK;
  double errNewQA;

  double *errx = new double [nK];
  double *erryK = new double [nK];
  double *erryA = new double [nK];
  for (int ip=0; ip<nK; ip++){
    errx[ip]=(xK[1]-xK[0])*0.5;
    erryK[ip]=errQK;
    erryA[ip]=errQA;
  }
      
 
  int loc;
  double peak=0;       

  // std::cout << " Cathode field " << tTheory[0] << std::endl;
  
  for (int ip=nK-2; ip>0; ip--){
    if (xK[ip]<0) break;
    // if (yK[ip]<peak && (xK[ip]>(tTheory[0]-100.E-6)) && (xK[ip]<(tTheory[0]+100.E-6)) ){
    if (yK[ip]<peak){  
      peak = yK[ip];
      loc = ip;
    }
  }
      
  tGK = xK[loc];
  fittedK = yK[loc];
  // std::cout << "K time and K " << fittedKtime << " " << fittedK << std::endl;
  tK = 0.;
      
  double tauLife = 50.e-6;
  double factorK = 1/(1/tauLife-1/tauelecK);
  double Q0K = fittedK*tGK/(exp(-(tGK)/tauelecK)-exp(-(tGK)/tauLife))/factorK;

  TF1 *funcK = new TF1("funcK", greenFunction, xK[0], xK[nK-1], 6);
  funcK->SetParameters(tauelecK, tauLife, tK, fittedK, tGK, 0.);
  funcK->SetParLimits(0,tauelecK*0.5, tauelecK*1.5);
  //funcK->SetParLimits(1, 0., 0.001);
  funcK->FixParameter(1, 0.001);
  funcK->SetParLimits(2, -1.e-5, +1.e-5);
  funcK->SetParLimits(3, Q0K*3, 0);
  funcK->SetParLimits(4, xK[0], xK[nK-1]);
  funcK->SetParLimits(5, -5., 5.);

  funcK->SetParName(0, "TauEl (s)");
  funcK->SetParName(1, "TauLife (s)");
  funcK->SetParName(2, "TK (s)");
  funcK->SetParName(3, "Q0_K (mV)");
  funcK->SetParName(4, "TGK (s)");
  funcK->SetParName(5, "Baseline (mV)");

  funcK->SetLineColor(kMagenta);
  funcK->SetLineWidth(2);

  //  std::cout << "Initial pars are " << tauelecK << " " << tauLife << " " << tK << " " <<  Q0K << " " <<  tGK  << std::endl;

  TGraphErrors *gKerr = new TGraphErrors(nK, xK, yK, errx, erryK);
  TFitResultPtr resultK = gKerr->Fit("funcK","QRES", "", xK[0], xK[nK-1]);
  resultK->Print();

  double lifeCathodeOnly = 0.;
  double errLifeCathode = 0.;

  if (resultK->Status()==0){
    
    //    funcK->Draw("same");
    gK->GetListOfFunctions()->Add(funcK);    

    tauLife = funcK->GetParameter(1);
    if (tauLife<0.000099){
      lifeCathodeOnly = tauLife;
      errLifeCathode = funcK->GetParError(1);
    } else {
      funcK->FixParameter(1, tauLife);
    }
    tK = funcK->GetParameter(2);
    errtK = funcK->GetParError(2);
    tGK = funcK->GetParameter(4);
    errtGK = funcK->GetParError(4);
    lifetime[10] = funcK->GetParameter(0);

    QK = funcK->GetMinimum();
    errQK = QK*(funcK->GetParError(3)/funcK->GetParameter(3));

    newQK = funcK->GetParameter(3);
    errNewQK = funcK->GetParError(3);

    // double xxx[1] = { tGK };
    // double errmin[1];  // error on the function at point x0

    // resultK->GetConfidenceIntervals(1, 1, 1, xxx, errmin, 0.683, true);//false);
    // errQK = errmin[0];
    
    // std::cout << "Error on QK is " << xxx[0] << " " << funcK->Eval(xxx[0]) << " " << errQK << std::endl;
    
    // TMatrixDSym cov = resultK->GetCovarianceMatrix(); 
    // cov.Print();
  }

    
  double tempx, tempy;
  
  // cout << " Anode field " << fields[2] << endl;
  int nA      = gA->GetN();
  double *xA  = gA->GetX();
  double *yA  = gA->GetY();


  loc = TMath::LocMax(nA,yA);
  peak = -99999.;
  for (int ip=nA-2; ip>0; ip--){
    if (xA[ip]<20.e-6) break;
    if (yA[ip]>peak){
      peak = yA[ip];
      loc = ip;
    }
  }
	
      
  tA = xA[loc];
  fittedA = yA[loc];

  for (int ip=loc; ip>0; ip--){
    tempx = xA[ip];
    tempy = yA[ip];
    if (xA[ip]<20.e-6) break;
    if (tempy<0.01*fittedA){
      // cout << ip << " " << tempx << " " << tempy << " " << endl;
      tGA = tempx;
      break;
    }
  }



  bool fitAnode = true;

  if (fitAnode){

    // TF1 *funcA = new TF1("funcA",fittingFunction,tTheory[0], xA[nA-1],4);
    // funcA->SetParameters(20, tauelecA*1E6, 4, (tTheory[0]+tTheory[1])*1e6);
    // // funcA->FixParameter(0, 5);
    // // funcA->FixParameter(1, 43);
    // // funcA->FixParameter(2, 0.002);
    // // funcA->FixParameter(3, (tTheory[0]+tTheory[1])*1e6);
    // // funcA->SetParLimits(0, 1, 50);
    // funcA->SetParLimits(1, tauelecA*1E6*0.7, tauelecA*1E6*1.3);
    // // funcA->SetParLimits(2, 0, 0.1);
    // funcA->SetParLimits(3, tTheory[0]*1e6, 1000);
    // funcA->SetParName(0, "#sigma");
    // funcA->SetParName(1, "#tau_{D}");
    // funcA->SetParName(2, "a");
    // funcA->SetParName(3, "t_0");
    // gA->Fit("funcA", "RQ", "", tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1]);    
    
    // for (int ip=0; ip<nA; ip++){
    //   tempx = xA[ip];
    //   tempy = funcA->Eval(xA[ip]);
    //   if (tempy>0.01*fittedA){
    // 	// cout << ip << " " << tempx << " " << tempy << " " << endl;
    // 	fittedAstartTime = tempx;
    // 	break;
    //   }
    // }
    // fittedA     = funcA->GetMaximum();
    // fittedAtime = funcA->GetMaximumX();
    
    tauLife = 50e-6;
    double factorA = 1/(1/tauLife-1/tauelecA);
    double Q0A = fittedA*(tA-tGA)/(exp(-(tA)/tauelecA)-exp(-(tA)/tauLife))/factorA;

    if ( (tA<=tGA) || (tGA<=tTheory[0]+tTheory[1]-0.1E-3 ) || (tA<=tTheory[0]+tTheory[1]-0.1E-3)  ){
      tGA = tTheory[0]+tTheory[1];
      tA = tTheory[0]+tTheory[1]+tTheory[2];
    }

    TF1 *funcA = new TF1("funcA", greenFunction, tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1], 6);
    funcA->SetParameters(tauelecA, tauLife, tGA, fittedA, tA, 0.);
    funcA->SetParLimits(0,tauelecA*0.5, tauelecA*1.5);
    //funcA->SetParLimits(1, 0., 0.001);
    funcA->FixParameter(1, 0.001);
    funcA->SetParLimits(2, tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1]);
    funcA->SetParLimits(3, 0., Q0A*3);
    funcA->SetParLimits(4, tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1]);
    funcA->SetParLimits(5, -5., 5.);

    funcA->SetParName(0, "TauEl (s)");
    funcA->SetParName(1, "TauLife (s)");
    funcA->SetParName(2, "TGA (s)");
    funcA->SetParName(3, "Q0_A (mV)");
    funcA->SetParName(4, "TA (s)");
    funcA->SetParName(5, "Baseline (mV)");

    funcA->SetLineColor(kMagenta);

    //    std::cout << "Initial pars are " << tauelecA << " " << tauLife << " " << tGA << " " <<  Q0A << " " <<  tA  << std::endl;

    TGraphErrors *gAerr = new TGraphErrors(nA, xA, yA, errx, erryA);
    TFitResultPtr resultA = gAerr->Fit("funcA","QRES", "", tTheory[0]+tTheory[1]-0.1E-3, xA[nA-1]);
    resultA->Print();

    if (resultA->Status()==0 && (funcA->GetParameter(4)>funcA->GetParameter(2)) ){
      
      //      funcA->Draw("same");
      gA->GetListOfFunctions()->Add(funcA);    

      tGA = funcA->GetParameter(2);
      errtGA = funcA->GetParError(2);
      tA = funcA->GetParameter(4);
      errtA = funcA->GetParError(4);

      QA = funcA->GetMaximum();
      errQA = QA*(funcA->GetParError(3)/funcA->GetParameter(3));
      
      newQA = funcA->GetParameter(3);
      errNewQA = funcA->GetParError(3);

      lifetime[11] = funcA->GetParameter(0);

      // double xxx[1] = { tA };
      // double errmax[1];

      // resultA->GetConfidenceIntervals(1, 1, 1, xxx, errmax, 0.683, false);
      // errQA = errmax[0];
      // std::cout << "Error on QA is " << xxx[0] << " " << funcA->Eval(xxx[0]) << " " << errQA << std::endl;
      // TMatrixDSym cov = resultK->GetCovarianceMatrix();
      // cov.Print();     
      // std::cout << "Fitted pars " << fittedA << " " << tA << " " << tGA << " " << funcA->GetParameter(2) << std::endl;
      
    } else {
      

      std::cout << tA << " " << tGA << std::endl;
      errtA=errtGA=3.e-6;
      QA = newQA = fittedA;
      errNewQA = errQA = 1.5;
    }
    
    
   
  }
      
          
      
  // std::cout << " Expected vs measured time at cathode " << tTheory[0] << " " << fittedKtime - fittedKstartTime << std::endl;
  // std::cout << " Expected vs measured time at anode " << tTheory[2] << " " << fittedAtime - fittedAstartTime << std::endl;
      
    
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
  t1 = tGK - tK;
  t2 = tGA - tGK;
  t3 = tA - tGA;
  R = 0;
  
  errt1 = t1*TMath::Sqrt( (errtGK/tGK)*(errtGK/tGK) + (errtK/tK)*(errtK/tK));
  errt2 = t3*TMath::Sqrt( (errtGK/tGK)*(errtGK/tGK) + (errtGA/tGA)*(errtGA/tGA));
  errt3 = t3*TMath::Sqrt( (errtGA/tGA)*(errtGA/tGA) + (errtA/tA)*(errtA/tA));

  if (t2>(tTheory[1]+100.e-6)) {
    errt2 = t2 = tTheory[1];
  }

  if (t3>(3*tTheory[2])){
    errt3 = t3 = tTheory[2];
  }



  if (tK<tGK && tGK<tGA && tGA<tA){
    // All good
  } else {
    std::cout << "THESE NUMBERS DON'T MAKE SENSE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << tK << " " << tGK << " " << tGA << " " << tA << std::endl;
    lifetime[2] = lifetime[4] = QK;
    lifetime[3] = lifetime[5] = QA;

    lifetime[6] = t1;
    lifetime[7] = t2;
    lifetime[8] = t3;
    lifetime[9] = 0.;

    lifetime[0]=lifetime[1] = 0;
    return -1;
  }
  
  
  
  //  QK *= gainAoverK;
  
  double taulife = 0.001;
  double Kcorrection = 0.;
  double Acorrection = 0.;
  
  lifetime[2] = QK;
  lifetime[3] = QA;
  lifetime[6] = t1;
  lifetime[7] = t2;
  lifetime[8] = t3;
  lifetime[9] = lifeCathodeOnly;

  lifeErrors[2] = errQK;
  lifeErrors[3] = errQA;
  lifeErrors[4] = errNewQK;
  lifeErrors[5] = errNewQA;
  lifeErrors[6] = errt1;
  lifeErrors[7] = errt2;
  lifeErrors[8] = errt3;
  lifeErrors[9] = errLifeCathode;

  if (R<1){
    
    int count = 0;

    lifetime[4] = newQK;
    lifetime[5] = newQA;
    R = TMath::Abs(newQA/newQK);
    lifetime[0] = (1/TMath::Abs(TMath::Log(R)))*(t2 + 0.5*(t1+t3)); 

    // while(count < 20){
      
    //   //cout << "Catodo " << endl;
    //   Kcorrection = getCorrectionFactor(t1, tauelecK, taulife);
    //   //cout << "Anodo " << endl;
    //   Acorrection = getCorrectionFactor(t3, tauelecA, taulife);
      
    //   newQK = QK/Kcorrection;
    //   newQA = QA/Acorrection;
      
    //   lifetime[4] = newQK;
    //   lifetime[5] = newQA;

    //   R =  TMath::Abs(newQA/newQK);
      
    //   lifetime[0] = (1/TMath::Abs(TMath::Log(R)))*(t2 + 0.5*(t1+t3));
      
    //   // std::cout <<" This is my lifetime " <<  lifetime[0] << std::endl;
      
    //   if (TMath::Abs(lifetime[0]-taulife)<0.001e-6) break;
    //   taulife=lifetime[0];
    //   count++;
    // }
    // std::cout << "Number of iterations : " << count << std::endl;
    
    TF1 f ("lifetime function", "([1]/[3])*(TMath::SinH([3]/(2*x))/TMath::SinH([1]/(2*x)))*TMath::Exp(-([2]+0.5*([1]+[3]))/x) - [0]", -0.1, 0.1);
    f.SetParameters(R, t1, t2, t3);
    errR = R*TMath::Sqrt((errNewQA/QA)*(errNewQA/QA)+ (errNewQK/QK)*(errNewQK/QK));
    
    ROOT::Math::WrappedTF1 wf1(f);
	
    // Create the Integrator
    ROOT::Math::BrentRootFinder brf;
	
    // Set parameters of the method
    brf.SetFunction( wf1, -0.1, 0.1 );
    brf.Solve();
	
    lifetime[1] =  TMath::Abs(brf.Root() ) ;
	

    // calculate syst errors, let's give it a try

    double newLife, lifeErr, tmpErrPlus, tmpErrMinus;
    lifeErr=0;

    f.SetParameters(R-errR, t1, t2, t3);
    brf.Solve();
    newLife      = TMath::Abs(brf.Root() ) ;
    tmpErrMinus  = TMath::Abs(lifetime[1]-newLife);
    
    f.SetParameters(R+errR, t1, t2, t3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrPlus  = TMath::Abs(lifetime[1]-newLife);
    lifeErr += TMath::Power( TMath::Max(tmpErrPlus, tmpErrMinus), 2);
    // std::cout << tmpErrMinus << " " << tmpErrPlus << " " << TMath::Sqrt(lifeErr) << std::endl;

    f.SetParameters(R, t1+errt1, t2, t3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrPlus  = TMath::Abs(lifetime[1]-newLife);

    f.SetParameters(R, t1-errt1, t2, t3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrMinus  = TMath::Abs(lifetime[1]-newLife);
    lifeErr += TMath::Power( TMath::Max(tmpErrPlus, tmpErrMinus), 2);
    // std::cout << tmpErrMinus << " " << tmpErrPlus << " " << TMath::Sqrt(lifeErr) << std::endl;

    f.SetParameters(R, t1, t2+errt2, t3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrPlus  = TMath::Abs(lifetime[1]-newLife);

    f.SetParameters(R, t1, t2-errt2, t3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrMinus  = TMath::Abs(lifetime[1]-newLife);
    lifeErr += TMath::Power( TMath::Max(tmpErrPlus, tmpErrMinus), 2);
    // std::cout << tmpErrMinus << " " << tmpErrPlus << " " << TMath::Sqrt(lifeErr) << std::endl;

    f.SetParameters(R, t1, t2, t3+errt3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrPlus  = TMath::Abs(lifetime[1]-newLife);

    f.SetParameters(R, t1, t2, t3-errt3);
    brf.Solve();
    newLife = TMath::Abs(brf.Root() ) ;
    tmpErrMinus  = TMath::Abs(lifetime[1]-newLife);
    lifeErr += TMath::Power( TMath::Max(tmpErrPlus, tmpErrMinus), 2);
    
    // std::cout << tmpErrMinus << " " << tmpErrPlus << " " << TMath::Sqrt(lifeErr) << std::endl;

    lifeErr = TMath::Sqrt(lifeErr);

    lifeErrors[1] = lifeErr;
    lifeErrors[0] = lifetime[0]*lifeErrors[1]/lifetime[1];

    std::cout << "R and err " << R << " +/- " << errR << std::endl;
    std::cout << "Life and err " << lifetime[1]*1e6 << " +/- " << lifeErrors[1]*1e6 << std::endl;



  } else {
    lifetime[0] = lifetime[1] = 0;
  }

  
  // printf("tK     : %12.4e \n",  tK  );
  // printf("tGK    : %12.4e \n",  tGK );
  // printf("tGA    : %12.4e \n",  tGA );
  // printf("tA     : %12.4e \n",  tA  );
  // printf("t1     : %12.4e \n",  t1  );
  // printf("t2     : %12.4e \n",  t2  );
  // printf("t3     : %12.4e \n",  t3  );
  // printf("QA     : %12.4e \n",  QA  );
  // printf("QK     : %12.4e \n",  QK  );
  // printf("QAcorr : %12.4e \n",  newQA  );
  // printf("QKcorr : %12.4e \n",  newQK  );
  // printf("R      : %12.4e \n",  R   );
  // printf("lifetime : %12.4e \n",  lifetime[0]  );
  // printf("lifetime2: %12.4e \n",  lifetime[1] );


  if (saveCanvas){
    TCanvas *c2 = new TCanvas("c2");
    gK->GetYaxis()->SetRangeUser(QK*1.2, QA*1.2);
    gK->Draw("Al");
    gA->Draw("l");

    c2->Print("Lastlifetime.png");
    c2->Print("Lastlifetime.root");
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
  if (nsmooth<0) nsmooth=10;
  return  nsmooth;

}

Double_t UsefulFunctions::greenFunction(Double_t *x, Double_t *par) {
  Double_t F_elec=par[0];
  Double_t T_life=par[1];
  Double_t t0=par[2];
  Double_t Q0=par[3];
  Double_t t1=par[4];
  Double_t baseline=par[5];
  if(x[0]<t0) return baseline;
  Double_t start=1./((t1-t0)*((1./T_life)-(1./F_elec)));
  Double_t fallExp=TMath::Exp(-1*(x[0]-t0)/F_elec);
  Double_t riseExp=TMath::Exp(-1*(x[0]-t0)/T_life);
  //Double_t Q0=par[3]/(start*(TMath::Exp(-1*(t1-t0)/F_elec) - TMath::Exp(-1*(t1-t0)/T_life)));
  if(x[0]<t1)
    return baseline+Q0*start*(fallExp-riseExp);
  Double_t bigExp=fallExp*(1 - TMath::Exp(-1*(t1-t0)*(1./T_life - 1./F_elec)));
  return baseline+Q0*start*bigExp;
}

Double_t UsefulFunctions::funcXenon(Double_t *x, Double_t *par) {
  Double_t baseline=par[5];
  Double_t Q=par[3];
  Double_t t0=par[2];
  Double_t riseTime=par[4]-t0;
  Double_t fallElec=par[0];
  //  Double_t isCathode=par[5];
  return baseline+(Q/(1. + TMath::Exp(-1*(x[0]-t0)/riseTime)))*TMath::Exp(-1*(x[0]-t0)/fallElec);
}
