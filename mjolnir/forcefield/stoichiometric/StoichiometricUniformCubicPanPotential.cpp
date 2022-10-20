#include <mjolnir/forcefield/stoichiometric/StoichiometricUniformCubicPanPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class StoichiometricUniformCubicPanPotential<double>;
template class StoichiometricUniformCubicPanPotential<float >;
} // mjolnir
