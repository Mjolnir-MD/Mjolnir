#include <mjolnir/input/read_input_file.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<SimulatorBase> read_parallelism<double, UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
template std::unique_ptr<SimulatorBase> read_parallelism<float , UnlimitedBoundary       >(const toml::value& root, const toml::value& simulator);
template std::unique_ptr<SimulatorBase> read_parallelism<double, CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);
template std::unique_ptr<SimulatorBase> read_parallelism<float , CuboidalPeriodicBoundary>(const toml::value& root, const toml::value& simulator);

template std::unique_ptr<SimulatorBase> read_boundary<double>(const toml::value& root, const toml::value& simulator);
template std::unique_ptr<SimulatorBase> read_boundary<float >(const toml::value& root, const toml::value& simulator);
} // mjolnir
