#include <mjolnir/input/read_spatial_partition.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// D-H
template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>        > read_spatial_partition<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>        , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >        > read_spatial_partition<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >        , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>        > read_spatial_partition<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>        , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >        > read_spatial_partition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >        , true>(const toml::value&);

// EXV
template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>     > read_spatial_partition<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>     , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >     > read_spatial_partition<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >     , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>     > read_spatial_partition<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>     , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >     > read_spatial_partition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >     , true>(const toml::value&);

// L-J
template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>       > read_spatial_partition<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>       , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >       > read_spatial_partition<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >       , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>       > read_spatial_partition<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>       , true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >       > read_spatial_partition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >       , true>(const toml::value&);

// UL-J
template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<double>> read_spatial_partition<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<double>, true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<float >> read_spatial_partition<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<float >, true>(const toml::value&);
template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>> read_spatial_partition<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, true>(const toml::value&);
template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >> read_spatial_partition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >, true>(const toml::value&);

} // mjolnir
