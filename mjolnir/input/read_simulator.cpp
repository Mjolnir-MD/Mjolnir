#include <mjolnir/input/read_simulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& data);

template std::unique_ptr<SimulatorBase> read_simulator_from_table<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator_from_table<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator_from_table<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
template std::unique_ptr<SimulatorBase> read_simulator_from_table<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

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

} // mjolnir
