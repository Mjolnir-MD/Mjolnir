#include <mjolnir/omp/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// D-H
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, PeriodicGridCellList <OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , PeriodicGridCellList <OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename DebyeHuckelPotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, VerletList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,                   typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , VerletList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,                   typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, VerletList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,            typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , VerletList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename DebyeHuckelPotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, UnlimitedBoundary>,         typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,         typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,  typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename DebyeHuckelPotential<float>::pair_parameter_type> >;

} // mjolnir
