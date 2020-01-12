#include <mjolnir/input/read_integrator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_velocity_verlet_integrator(const toml::value& simulator);
template VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_velocity_verlet_integrator(const toml::value& simulator);
template VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_velocity_verlet_integrator(const toml::value& simulator);
template VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_velocity_verlet_integrator(const toml::value& simulator);

template UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_underdamped_langevin_integrator(const toml::value& simulator);
template UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_underdamped_langevin_integrator(const toml::value& simulator);
template UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_underdamped_langevin_integrator(const toml::value& simulator);
template UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_underdamped_langevin_integrator(const toml::value& simulator);

template BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_BAOAB_langevin_integrator(const toml::value& simulator);
template BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_BAOAB_langevin_integrator(const toml::value& simulator);
template BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_BAOAB_langevin_integrator(const toml::value& simulator);
template BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_BAOAB_langevin_integrator(const toml::value& simulator);

template GFWNpTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_GFW_NpT_langevin_integrator(const toml::value& simulator);
template GFWNpTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_GFW_NpT_langevin_integrator(const toml::value& simulator);
template GFWNpTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_GFW_NpT_langevin_integrator(const toml::value& simulator);
template GFWNpTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_GFW_NpT_langevin_integrator(const toml::value& simulator);
} // mjolnir
