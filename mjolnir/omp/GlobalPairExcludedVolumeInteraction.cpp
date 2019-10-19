#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
