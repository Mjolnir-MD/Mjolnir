#include <mjolnir/input/read_global_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template GlobalForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_global_forcefield(const toml::array& interactions, const std::string& input_path);
template GlobalForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_global_forcefield(const toml::array& interactions, const std::string& input_path);
template GlobalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_global_forcefield(const toml::array& interactions, const std::string& input_path);
template GlobalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_global_forcefield(const toml::array& interactions, const std::string& input_path);
} // mjolnir
