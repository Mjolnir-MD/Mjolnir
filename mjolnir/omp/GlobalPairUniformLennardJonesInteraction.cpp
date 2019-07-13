#include <mjolnir/omp/GlobalPairUniformLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// ============================================================================
// Uniform L-J
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        typename UniformLennardJonesPotential<float>::parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, PeriodicGridCellList <OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , PeriodicGridCellList <OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, typename UniformLennardJonesPotential<float>::parameter_type> >;
// VerletList
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, VerletList<OpenMPSimulatorTraits<double, UnlimitedBoundary>,                   typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , VerletList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,                   typename UniformLennardJonesPotential<float>::parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, VerletList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,            typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , VerletList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,            typename UniformLennardJonesPotential<float>::parameter_type> >;
// Naive
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, UnlimitedBoundary>,         typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  UnlimitedBoundary>,         typename UniformLennardJonesPotential<float>::parameter_type> >;
template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, NaivePairCalculation<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>,  typename UniformLennardJonesPotential<double>::parameter_type>>;
template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> , NaivePairCalculation<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>,  typename UniformLennardJonesPotential<float>::parameter_type> >;

} // mjolnir
