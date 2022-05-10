#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LennardJonesPotential<double>;
template class LennardJonesPotential<float >;
} // mjolnir
