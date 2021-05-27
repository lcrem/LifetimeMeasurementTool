//////////////////////////////////////////////////////////////////////////////
/////  UsefulFunctions.h        Lifetime Calculation useful functions    /////
/////                                                                    /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
//////////////////////////////////////////////////////////////////////////////


#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H
#include <string>
#include <sstream>
#include "TGraph.h"
#include "TFile.h"


namespace UsefulFunctions {
  
  //! Just average - calculates the average of n graphs
  /**
   *
   * @param  numGraphs : number of graphs to average
   * @param  grPtrPtr : array of pointers to graphs
   * @return pointer to TGraph
   */
  TGraph *justAverage(Int_t numGraphs, TGraph **grPtrPtr);

  TGraph *justAverageFromFile(TFile *fin, Int_t numGraphs, int channel);

  //! Put baseline to 0 before averaging n graphs
  /**
   *
   * @param  numGraphs : number of graphs to average
   * @param  graphs : array of pointers to graphs
   * @return pointer to TGraph
   */
  TGraph *getZeroedAverage(Int_t numGraphs, TGraph **graphs);

  //! Average batches of nmax
  /**
   *
   * @param  graphs : graphs to average
   * @param  nmax : number of graphs to average
   * @param  g : output array of averaged graphs 
   */
  Int_t avgSomeGraphs(TGraph **graphs, int nmax, TGraph **g);

  Int_t avgSomeGraphs(TFile *fin, int nmax, TGraph **g, int ch);

  //! zero baseline
  /**
   *
   * @param  g: input/output TGraph
   */
  void zeroBaseline(TGraph *g);

  //! Calculates ICARUS polinomial to find drift velocity
  /**
   *
   * @param  E: electric field in V/cm
   */
  Double_t ICARUSpolynomial(Double_t E);

  //! Get correction factor coming from preamps (see PrM paper appendix for calculation)
  /**
   *
   * @param  fittedDriftTime : fitted drift time
   * @param  tauelec : preamps decay time in s
   * @param  taulife : expected electron lifetime in s
   */
  Double_t getCorrectionFactor(double fittedDriftTime, double tauelec, double taulife);

  //! PMT fitting function
  /**
   *
   * @param  x : array of time in s
   * @param par : [0] sigma, [1] : tau electronics in us, [2] : integrated charge, [3] : t0 in s
   */
  Double_t fittingFunction(Double_t *x, Double_t *par);


  //! Calculate lifetime
  /**
   *
   * @param gK : TGraph of cathode signal
   * @param gA : TGraph of anode signal
   * @param gSum : TGraph of anode+cathode signal
   * @param whichPrM : 0 is PrM1 and 1 is PrM2
   * @param tTheory : expected drift times between K-GK, GK-GA, GA-A
   * @param lifetime : calculated lifetime [0] approximation, [1] full formula
   * @param saveCanvas : boolean variable, save canvas or not
   */
  Int_t calculateLifetime(TGraph *gK, TGraph *gA, TGraph *gSum, int whichPrM, double tTheory[3], double lifetime[20], double lifeErrors[20], bool saveCanvas);

  //! Translate graph
  /**
   *
   * @param grWave : TGraph input output
   * @param deltaT : deltaT to translate
   */
  TGraph *translateGraph(const TGraph *grWave, Double_t deltaT);

  //! Get smoothing number - TODO GET RID OF THIS FUNCTION, LC: I don't like it
  /**
   *
   * @param deltat : sampling time
   * @param tdrift : expected drift time 
   */  
  Int_t getSmoothingNumber(double deltat, double tdrift);

  //! Green function
  /**
   *
   * @param x : time in s
   * @param par : [0] Gain x Q, [1] tau in us, [2] rise in us, [3] shift in us
   */  
  Double_t greenFunction(Double_t *x, Double_t *par);

  Double_t greenFunctionSum(Double_t *x, Double_t *par);

  Double_t funcXenon(Double_t *x, Double_t *par);

}

#endif //USEFULFUNCTIONS_H
