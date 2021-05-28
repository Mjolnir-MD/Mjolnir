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

template gBAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_gBAOAB_langevin_integrator(const toml::value& simulator);
template gBAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_gBAOAB_langevin_integrator(const toml::value& simulator);
template gBAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_gBAOAB_langevin_integrator(const toml::value& simulator);
template gBAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_gBAOAB_langevin_integrator(const toml::value& simulator);

template GJFNVTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_GJFNVT_langevin_integrator(const toml::value& simulator);
template GJFNVTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_GJFNVT_langevin_integrator(const toml::value& simulator);
template GJFNVTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_GJFNVT_langevin_integrator(const toml::value& simulator);
template GJFNVTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_GJFNVT_langevin_integrator(const toml::value& simulator);

template GFWNPTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       > read_GFW_NPT_langevin_integrator(const toml::value& simulator);
template GFWNPTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       > read_GFW_NPT_langevin_integrator(const toml::value& simulator);
template GFWNPTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_GFW_NPT_langevin_integrator(const toml::value& simulator);
template GFWNPTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_GFW_NPT_langevin_integrator(const toml::value& simulator);

} // mjolnir
