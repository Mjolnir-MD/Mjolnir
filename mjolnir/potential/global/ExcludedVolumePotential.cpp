#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
