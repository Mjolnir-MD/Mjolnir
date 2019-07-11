#include <mjolnir/input/read_external_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ImplicitMembranePotential<double> read_implicit_membrane_potential(const toml::value& external);
template ImplicitMembranePotential<float > read_implicit_membrane_potential(const toml::value& external);

template LennardJonesWallPotential<double> read_lennard_jones_wall_potential(const toml::value& external);
template LennardJonesWallPotential<float > read_lennard_jones_wall_potential(const toml::value& external);

template ExcludedVolumeWallPotential<double> read_excluded_volume_wall_potential(const toml::value& external);
template ExcludedVolumeWallPotential<float > read_excluded_volume_wall_potential(const toml::value& external);
}
