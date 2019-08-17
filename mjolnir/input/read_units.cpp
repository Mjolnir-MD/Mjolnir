#include <mjolnir/input/read_units.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<SimulatorBase> read_units<SimulatorTraits<double, UnlimitedBoundary       >>(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_units<SimulatorTraits<float,  UnlimitedBoundary       >>(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_units<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& data);
template std::unique_ptr<SimulatorBase> read_units<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& data);
} // mjolnir
