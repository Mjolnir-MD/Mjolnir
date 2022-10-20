#include <mjolnir/input/read_stoichiometric_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::pair<StoichiometricUniformCubicPanPotential<double>, StoichiometricEmptyCombinationRule<SimulatorTraits<double, UnlimitedBoundary>       , StoichiometricUniformCubicPanPotential<double>>> read_stoichiometric_uniform_cubic_pan_potential(const toml::value&, const std::string&, const std::string&);
template std::pair<StoichiometricUniformCubicPanPotential<float >, StoichiometricEmptyCombinationRule<SimulatorTraits<float,  UnlimitedBoundary>       , StoichiometricUniformCubicPanPotential<float >>> read_stoichiometric_uniform_cubic_pan_potential(const toml::value&, const std::string&, const std::string&);
}
