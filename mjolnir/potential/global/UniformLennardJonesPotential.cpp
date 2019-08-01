#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UniformLennardJonesPotential<double>;
template class UniformLennardJonesPotential<float>;
} // mjolnir
