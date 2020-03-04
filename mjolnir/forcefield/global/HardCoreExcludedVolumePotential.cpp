#include <mjolnir/potential/global/HardCoreExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HardCoreExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class HardCoreExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class HardCoreExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class HardCoreExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
