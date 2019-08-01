#include <mjolnir/potential/local/GaussianPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GaussianPotential<double>;
template class GaussianPotential<float>;
} // mjolnir
