//////////////////////////////////////////////////////////////////////////////
/////  LifetimeConventions.h        Lifetime Conventions                 /////
/////                                                                    /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
//////////////////////////////////////////////////////////////////////////////


#ifndef LIFETIMECONVENTIONS_H
#define LIFETIMECONVENTIONS_H
#include <string>

const std::string whereIsMyPrMdata="/Users/linda/DUNE/exampleDataFromDigitiser/";

// First Purity Monitor machined at UCL (PrM1)
Double_t PrM1distance[3] = {0.01823, 0.16424, 0.00985};
// ICARUS Purity Monitor refurbished (PrM2)
Double_t PrM2distance[3] = {0.016, 0.15, 0.085};


// pre-amplifiers tau electronics for PrM1 {K, A} and relative gain A/K TODO: APPLY CALIBRATION AND CHECK THESE NUMBERS !
Double_t PrM1preamp[3] = {91.5031*1e-6, 43.4835*1e-6, 0.89512};
// pre-amplifiers tau electronics for PrM1 {K, A} and relative gain A/K
Double_t PrM2preamp[3] = {90.*1e-6, 90.*1e-6, 1.};

#endif //LIFETIMECONVENTIONS_H
