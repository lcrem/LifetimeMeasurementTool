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
  
  Double_t timeStamp;  ///! Timestamp (this is determined by the record length, which is ~ 300 ns for the two waveform runs we acquired. I believe it’s per nanosecond). 
  Double_t ADC;        ///! ADC Input, mV
  Int_t baseline;      ///! Baseline
  Int_t trigger;       ///! Trigger
  Int_t longGate;      ///! Long Gate
  Int_t shortGate;     ///! Short Gate
  Int_t zero;          ///! N/A (all 0s)
  
  ClassDef(RawDigitiser, 2);
 
};

#endif
