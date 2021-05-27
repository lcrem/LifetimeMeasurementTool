#include <string>
#include "LifetimeConventions.h"
#include "Rtypes.h"

const std::string whereIsMyPrMdata="/Users/linda/DUNE/exampleDataFromDigitiser/";

// First Purity Monitor machined at UCL (PrM1)
const Double_t PrM1distance[3] = {0.01823, 0.16424, 0.00985};
// ICARUS Purity Monitor refurbished (PrM2)
const Double_t PrM2distance[3] = {0.018, 0.15, 0.01};

// pre-amplifiers tau electronics for PrM1 {K, A} and relative gain A/K , Preamp B on K, Preamp A on A
const Double_t PrM1preamp[3] = {319.*1e-6, 250.*1e-6, 1.};
// pre-amplifiers tau electronics for PrM2 {K, A} and relative gain A/K , Preamp D on K, Preamp C on A
const Double_t PrM2preamp[3] = {298.*1e-6, 240.*1e-6, 1.};

