#include <mjolnir/core/MolecularDynamicsSimulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// BAOAB
template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// Langevin
template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// VelVerlet
template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
