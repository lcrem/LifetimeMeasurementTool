/*
 * RawDigitiser.h
 * CAEN DT5730 output file format, raw digitiser input
 * 
 */

#ifndef RAW_DIGITISER
#define RAW_DIGITISER

#include "TObject.h"

class RawDigitiser: public TObject
{
 public:
  RawDigitiser();
  ~RawDigitiser();

  Double_t timeStamp;
  Double_t ADC;
  Int_t baseline;
  Int_t trigger;
  Int_t longGate;
  Int_t shortGate;
  Int_t zero;
  ClassDef(RawDigitiser, 1);
  
};

#endif
