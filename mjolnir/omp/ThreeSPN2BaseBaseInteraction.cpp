#include <mjolnir/omp/ThreeSPN2BaseBaseInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// CellList
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, PeriodicGridCellList <OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, PeriodicGridCellList <OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
// VerletList
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        VerletList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,                   typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        VerletList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,                   typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, VerletList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,            typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, VerletList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
// Naive
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        NaivePairCalculation<OpenMPSimulatorTraits<double, UnlimitedBoundary>,         typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        NaivePairCalculation<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,         typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, NaivePairCalculation<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,  typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, NaivePairCalculation<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;

} // mjolnir
