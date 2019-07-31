#include <mjolnir/interaction/global/GlobalPairExcludedVolumeInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// exv
// CellList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, PeriodicGridCellList <SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , PeriodicGridCellList <SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, VerletList<SimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;

} // mjolnir
