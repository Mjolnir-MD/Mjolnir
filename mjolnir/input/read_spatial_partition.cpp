#include <mjolnir/input/read_spatial_partition.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// D-H
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<double>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>        >> read_spatial_partition(const toml::value&);

template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<float>         >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float>         >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<float>         >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float>         >> read_spatial_partition(const toml::value&);

// EXV
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>     >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<double>     >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>     >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>     >> read_spatial_partition(const toml::value&);

template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<float>      >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float>      >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<float>      >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float>      >> read_spatial_partition(const toml::value&);

// L-J
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>       >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<double>       >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>       >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<double>       >> read_spatial_partition(const toml::value&);

template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<float>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<float>        >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float>        >> read_spatial_partition(const toml::value&);

// UL-J
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<double>>> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<double>>> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>> read_spatial_partition(const toml::value&);

template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<float> >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<float> >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> >> read_spatial_partition(const toml::value&);
template std::unique_ptr<SpatialPartitionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float> >> read_spatial_partition(const toml::value&);

} // mjolnir
