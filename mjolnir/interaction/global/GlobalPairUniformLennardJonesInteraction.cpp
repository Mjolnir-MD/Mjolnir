#include <mjolnir/interaction/global/GlobalPairUniformLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// ============================================================================
// Uniform L-J
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>,        typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>,        typename UniformLennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, PeriodicGridCellList <SimulatorTraits<double, CuboidalPeriodicBoundary>, typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , PeriodicGridCellList <SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename UniformLennardJonesPotential<float>::pair_parameter_type> >;
// VerletList
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, VerletList<SimulatorTraits<double, UnlimitedBoundary>,                   typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , VerletList<SimulatorTraits<float,  UnlimitedBoundary>,                   typename UniformLennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>,            typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename UniformLennardJonesPotential<float>::pair_parameter_type> >;
// Naive
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,         typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,         typename UniformLennardJonesPotential<float>::pair_parameter_type> >;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>,  typename UniformLennardJonesPotential<double>::pair_parameter_type>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename UniformLennardJonesPotential<float>::pair_parameter_type> >;

} // mjolnir
