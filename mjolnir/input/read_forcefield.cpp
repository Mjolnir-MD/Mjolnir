#include <mjolnir/input/read_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield_from_table(const toml::value& ff, const std::string& input_path);
template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield_from_table(const toml::value& ff, const std::string& input_path);

template ForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_forcefield(const toml::value& root, std::size_t N);
template ForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_forcefield(const toml::value& root, std::size_t N);
template ForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, std::size_t N);
template ForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_forcefield(const toml::value& root, std::size_t N);
} // mjolnir
