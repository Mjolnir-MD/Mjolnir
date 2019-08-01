#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// D-H
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>,        typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>,        typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, PeriodicGridCellList <SimulatorTraits<double, CuboidalPeriodicBoundary>, typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , PeriodicGridCellList <SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename DebyeHuckelPotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, VerletList<SimulatorTraits<double, UnlimitedBoundary>,                   typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , VerletList<SimulatorTraits<float,  UnlimitedBoundary>,                   typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>,            typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename DebyeHuckelPotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,         typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float> , NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,         typename DebyeHuckelPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>,  typename DebyeHuckelPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float> , NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename DebyeHuckelPotential<float>::pair_parameter_type> >;

} // mjolnir
