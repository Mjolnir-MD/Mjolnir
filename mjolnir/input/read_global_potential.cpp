#include <mjolnir/input/read_global_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_excluded_volume_potential(const toml::value& global);
template ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_excluded_volume_potential(const toml::value& global);
template ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_excluded_volume_potential(const toml::value& global);
template ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_excluded_volume_potential(const toml::value& global);

template HardCoreExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_hard_core_excluded_volume_potential(const toml::value& global);
template HardCoreExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_hard_core_excluded_volume_potential(const toml::value& global);
template HardCoreExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_hard_core_excluded_volume_potential(const toml::value& global);
template HardCoreExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_hard_core_excluded_volume_potential(const toml::value& global);

template LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_lennard_jones_potential(const toml::value& global);
template LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_lennard_jones_potential(const toml::value& global);
template LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_lennard_jones_potential(const toml::value& global);
template LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_lennard_jones_potential(const toml::value& global);

template UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_uniform_lennard_jones_potential(const toml::value& global);
template UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_uniform_lennard_jones_potential(const toml::value& global);
template UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_uniform_lennard_jones_potential(const toml::value& global);
template UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_uniform_lennard_jones_potential(const toml::value& global);

template DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);

template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_3spn2_excluded_volume_potential(const toml::value& global);
template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_3spn2_excluded_volume_potential(const toml::value& global);
template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_3spn2_excluded_volume_potential(const toml::value& global);
template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_3spn2_excluded_volume_potential(const toml::value& global);
}
