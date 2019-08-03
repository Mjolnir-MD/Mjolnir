#include <mjolnir/interaction/global/GlobalPairExcludedVolumeInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// exv
// CellList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> >;

} // mjolnir
