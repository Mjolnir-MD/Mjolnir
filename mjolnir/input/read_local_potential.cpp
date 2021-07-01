#include <mjolnir/input/read_local_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template std::vector<std::pair<std::array<std::size_t, 2>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 3>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 3>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);

template std::vector<std::pair<std::array<std::size_t, 2>, GoContactPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GoContactPotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GoContactAttractivePotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GoContactAttractivePotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GoContactRepulsivePotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GoContactRepulsivePotential<float >>> read_local_potential(const toml::value& local);


template std::vector<std::pair<std::array<std::size_t, 4>, ClementiDihedralPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, ClementiDihedralPotential<float >>> read_local_potential(const toml::value& local);

template std::vector<std::pair<std::array<std::size_t, 2>, GaussianPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 2>, GaussianPotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 3>, GaussianPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 3>, GaussianPotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, PeriodicGaussianPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, PeriodicGaussianPotential<float >>> read_local_potential(const toml::value& local);

template std::vector<std::pair<std::array<std::size_t, 3>, FlexibleLocalAnglePotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 3>, FlexibleLocalAnglePotential<float >>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, FlexibleLocalDihedralPotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, FlexibleLocalDihedralPotential<float >>> read_local_potential(const toml::value& local);

// template std::vector<std::pair<std::array<std::size_t, 2>, ThreeSPN2BondPotential<double>>> read_local_potential(const toml::value& local);
// template std::vector<std::pair<std::array<std::size_t, 2>, ThreeSPN2BondPotential<float >>> read_local_potential(const toml::value& local);

template std::vector<std::pair<std::array<std::size_t, 4>, CosinePotential<double>>> read_local_potential(const toml::value& local);
template std::vector<std::pair<std::array<std::size_t, 4>, CosinePotential<float >>> read_local_potential(const toml::value& local);
} // mjolnir
