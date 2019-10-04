#include <mjolnir/input/read_system.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template System<SimulatorTraits<double, UnlimitedBoundary>       > read_system<SimulatorTraits<double, UnlimitedBoundary       >>(const toml::value& data, std::size_t N);
template System<SimulatorTraits<float,  UnlimitedBoundary>       > read_system<SimulatorTraits<float,  UnlimitedBoundary       >>(const toml::value& data, std::size_t N);
template System<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_system<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& data, std::size_t N);
template System<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_system<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& data, std::size_t N);
} // mjolnir
