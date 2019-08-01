#include <mjolnir/input/read_spatial_partition.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);

} // mjolnir
