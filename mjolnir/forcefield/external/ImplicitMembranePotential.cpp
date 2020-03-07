#include <mjolnir/forcefield/external/ImplicitMembranePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ImplicitMembranePotential<double>;
template class ImplicitMembranePotential<float>;
} // mjolnir
