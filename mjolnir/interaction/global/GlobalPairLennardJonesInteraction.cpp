#include <mjolnir/interaction/global/GlobalPairLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// L-J
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>,        typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>,        typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, PeriodicGridCellList <SimulatorTraits<double, CuboidalPeriodicBoundary>, typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , PeriodicGridCellList <SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename LennardJonesPotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, VerletList<SimulatorTraits<double, UnlimitedBoundary>,                   typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , VerletList<SimulatorTraits<float,  UnlimitedBoundary>,                   typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>,            typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename LennardJonesPotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,         typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,         typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>,  typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename LennardJonesPotential<float>::pair_parameter_type> >;

} // mjolnir
