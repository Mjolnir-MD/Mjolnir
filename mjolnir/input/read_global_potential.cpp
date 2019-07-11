#include <mjolnir/input/read_global_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ExcludedVolumePotential<double> read_excluded_volume_potential(const toml::value& global);
template ExcludedVolumePotential<float > read_excluded_volume_potential(const toml::value& global);

template LennardJonesPotential<double> read_lennard_jones_potential(const toml::value& global);
template LennardJonesPotential<float > read_lennard_jones_potential(const toml::value& global);

template UniformLennardJonesPotential<double> read_uniform_lennard_jones_potential(const toml::value& global);
template UniformLennardJonesPotential<float > read_uniform_lennard_jones_potential(const toml::value& global);

template DebyeHuckelPotential<double> read_debye_huckel_potential(const toml::value& global);
template DebyeHuckelPotential<float > read_debye_huckel_potential(const toml::value& global);
}
