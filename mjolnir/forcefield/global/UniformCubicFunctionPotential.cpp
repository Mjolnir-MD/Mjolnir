#include <mjolnir/forcefield/global/UniformCubicFunctionPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error  "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UniformCubicFunctionPotential<double>;
template class UniformCubicFunctionPotential<float>;
} // mjolnir
