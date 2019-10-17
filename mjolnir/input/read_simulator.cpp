#include <mjolnir/input/read_simulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

template RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>       > read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>       > read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_rng(const toml::value&);
} // mjolnir
