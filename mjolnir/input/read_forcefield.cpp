#include <mjolnir/input/read_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield(const toml::value& root, const std::size_t N);
template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield(const toml::value& root, const std::size_t N);
template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, const std::size_t N);
template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, const std::size_t N);
} // mjolnir
