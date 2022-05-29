#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalStoichiometricInteractionPotential<double>;
template class GlobalStoichiometricInteractionPotential<float >;
} // mjolnir
