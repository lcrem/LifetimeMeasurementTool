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

  Int_t run;

  ClassDef(RawDigitiser, 1);
  
};

#endif
