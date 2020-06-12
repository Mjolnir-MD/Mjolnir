#include <mjolnir/input/read_forcefield.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_default_forcefield(const toml::value&, const std::size_t);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_default_forcefield(const toml::value&, const std::size_t);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_default_forcefield(const toml::value&, const std::size_t);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_default_forcefield(const toml::value&, const std::size_t);

template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_multiple_basin_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_multiple_basin_forcefield(const toml::value&, const toml::value&);

template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const toml::value&);
template std::unique_ptr<ForceFieldBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_forcefield(const toml::value&, const toml::value&);
} // mjolnir
