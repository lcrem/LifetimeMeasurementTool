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
Double_t PrM2distance[3] = {0.018, 0.15, 0.01};


// pre-amplifiers tau electronics for PrM1 {K, A} and relative gain A/K , Preamp B on K, Preamp A on A
Double_t PrM1preamp[3] = {101.*1e-6, 269.*1e-6, 0.9};
// pre-amplifiers tau electronics for PrM2 {K, A} and relative gain A/K , Preamp D on K, Preamp C on A
Double_t PrM2preamp[3] = {101.*1e-6, 262.*1e-6, 1.1};

#endif //LIFETIMECONVENTIONS_H
