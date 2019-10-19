#include <mjolnir/forcefield/3SPN2/ThreeSPN2ExcludedVolumePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
