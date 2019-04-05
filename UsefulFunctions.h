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

namespace UsefulFunctions {
  
  TGraph *justAverage(Int_t numGraphs, TGraph **grPtrPtr);

  TGraph *getZeroedAverage(Int_t numGraphs, TGraph **graphs);

  int avgSomeGraphs(std::string filename, int nmax, TGraph **g);

  void zeroBaseline(TGraph *g);

}

#endif //USEFULFUNCTIONS_H
