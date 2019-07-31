#include <mjolnir/omp/GlobalPairLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ============================================================================
// L-J
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, PeriodicGridCellList <OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , PeriodicGridCellList <OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename LennardJonesPotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, VerletList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,                   typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , VerletList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,                   typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, VerletList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,            typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , VerletList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename LennardJonesPotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, UnlimitedBoundary>,         typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,         typename LennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,  typename LennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename LennardJonesPotential<float>::pair_parameter_type> >;

} // mjolnir
