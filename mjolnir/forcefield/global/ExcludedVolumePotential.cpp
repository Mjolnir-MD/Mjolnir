#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExcludedVolumeParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ExcludedVolumeParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ExcludedVolumeParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ExcludedVolumeParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
