#include <mjolnir/input/read_external_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ExternalForceField<SimulatorTraits<double, UnlimitedBoundary>       > read_external_forcefield(const toml::array& interactions, const std::string& input_path);
template ExternalForceField<SimulatorTraits<float,  UnlimitedBoundary>       > read_external_forcefield(const toml::array& interactions, const std::string& input_path);
template ExternalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_external_forcefield(const toml::array& interactions, const std::string& input_path);
template ExternalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_external_forcefield(const toml::array& interactions, const std::string& input_path);
} // mjolnir
