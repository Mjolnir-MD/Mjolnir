#include <mjolnir/forcefield/global/UniformCubicPanPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error  "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UniformCubicPanPotential<double>;
template class UniformCubicPanPotential<float>;
} // mjolnir
