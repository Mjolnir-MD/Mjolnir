#include <mjolnir/input/read_local_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template LocalForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_local_forcefield(const toml::array& interactions, const std::string& input_path);
template LocalForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_local_forcefield(const toml::array& interactions, const std::string& input_path);
template LocalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_local_forcefield(const toml::array& interactions, const std::string& input_path);
template LocalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_local_forcefield(const toml::array& interactions, const std::string& input_path);
} // mjolnir
