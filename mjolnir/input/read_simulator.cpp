#include <mjolnir/input/read_simulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ----------------------------------------------------------------------------
// read_molecular_dynamics_simulator

template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_simulated_annealing_simulator

template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_steepest_descent_simulator

template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_switching_forcefield_simulator

template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_simulator
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_integrator_type

template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);


} // mjolnir
