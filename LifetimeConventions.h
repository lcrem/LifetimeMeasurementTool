//////////////////////////////////////////////////////////////////////////////
/////  LifetimeConventions.h        Lifetime Conventions                 /////
/////                                                                    /////
/////  Author: Linda Cremonesi (l.cremonesi@ucl.ac.uk)                   /////
//////////////////////////////////////////////////////////////////////////////


#ifndef LIFETIMECONVENTIONS_H
#define LIFETIMECONVENTIONS_H
#include <string>
#include "Rtypes.h"

extern const std::string whereIsMyPrMdata;

// First Purity Monitor machined at UCL (PrM1)
extern const Double_t PrM1distance[3];
// ICARUS Purity Monitor refurbished (PrM2)
extern const Double_t PrM2distance[3];

// pre-amplifiers tau electronics for PrM1 {K, A} and relative gain A/K , Preamp B on K, Preamp A on A
extern const Double_t PrM1preamp[3];
// pre-amplifiers tau electronics for PrM2 {K, A} and relative gain A/K , Preamp D on K, Preamp C on A
extern const Double_t PrM2preamp[3];

const double gridTransparencyRenorm = 1.0; //1.16;

#endif //LIFETIMECONVENTIONS_H
