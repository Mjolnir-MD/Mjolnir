#include <mjolnir/forcefield/global/HardCoreExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class HardCoreExcludedVolumeParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class HardCoreExcludedVolumeParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class HardCoreExcludedVolumeParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class HardCoreExcludedVolumeParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
