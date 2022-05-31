#include <mjolnir/forcefield/stoichiometric/StoichiometricCubicPanPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class StoichiometricCubicPanPotential<double>;
template class StoichiometricCubicPanPotential<float >;
} // mjolnir
