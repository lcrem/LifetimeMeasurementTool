#include "UsefulFunctions.h"

#include "TFile.h"
#include "TMath.h"
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


int UsefulFunctions::avgSomeGraphs(std::string filename, int nmax, TGraph **g){

  TFile *f = new TFile(filename.c_str(), "read");
  TGraph *graphs[1000];
  int count = 0;
  int ntot = 0;

  double newy[10005];
  double newx[10005];

  for (int i=0; i<1000; i++){
    // cout << i << endl;
    graphs[count] = (TGraph*) f->Get(Form("graph%i", i+1));
    if(!graphs[count]) break;
    count++;
    if (count==nmax){
      
      g[ntot] = getZeroedAverage( count, graphs );
      zeroBaseline(g[ntot]);
      /* g[ntot]  = smoothGraph(g[ntot],  5);  */
      
      count = 0;
      ntot++;
      // delete gtemp;
    }
  }

  f->Close();

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
