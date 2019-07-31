#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// exv
// CellList
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, PeriodicGridCellList <OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , PeriodicGridCellList <OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, VerletList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , VerletList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, VerletList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , VerletList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename ExcludedVolumePotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ExcludedVolumePotential<float>::pair_parameter_type> >;

} // mjolnir
