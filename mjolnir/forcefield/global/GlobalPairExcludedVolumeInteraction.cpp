#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// exv
// CellList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
