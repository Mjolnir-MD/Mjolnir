#include <mjolnir/forcefield/3SPN2/ThreeSPN2BondPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2BondPotential<double>;
template class ThreeSPN2BondPotential<float>;
} // mjolnir
