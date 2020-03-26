#include <mjolnir/forcefield/local/PeriodicGaussianPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class PeriodicGaussianPotential<double>;
template class PeriodicGaussianPotential<float>;
} // mjolnir
