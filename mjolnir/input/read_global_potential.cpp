#include <mjolnir/input/read_global_potential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>> read_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >> read_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>> read_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >> read_excluded_volume_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , InversePowerPotential<double>> read_inverse_power_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , InversePowerPotential<float >> read_inverse_power_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, InversePowerPotential<double>> read_inverse_power_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, InversePowerPotential<float >> read_inverse_power_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , HardCoreExcludedVolumePotential<double>> read_hard_core_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , HardCoreExcludedVolumePotential<float >> read_hard_core_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<double>> read_hard_core_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<float >> read_hard_core_excluded_volume_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>> read_lennard_jones_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >> read_lennard_jones_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>> read_lennard_jones_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >> read_lennard_jones_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesAttractivePotential<double>> read_lennard_jones_attractive_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesAttractivePotential<float >> read_lennard_jones_attractive_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesAttractivePotential<double>> read_lennard_jones_attractive_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesAttractivePotential<float >> read_lennard_jones_attractive_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , WCAPotential<double>> read_wca_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , WCAPotential<float >> read_wca_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<double>> read_wca_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<float >> read_wca_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>> read_debye_huckel_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >> read_debye_huckel_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>> read_debye_huckel_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >> read_debye_huckel_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , ThreeSPN2ExcludedVolumePotential<double>> read_3spn2_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , ThreeSPN2ExcludedVolumePotential<float >> read_3spn2_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ThreeSPN2ExcludedVolumePotential<double>> read_3spn2_excluded_volume_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ThreeSPN2ExcludedVolumePotential<float >> read_3spn2_excluded_volume_potential(const toml::value&);

template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , iSoLFAttractivePotential<double>> read_isolf_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , iSoLFAttractivePotential<float >> read_isolf_potential(const toml::value&);
template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, iSoLFAttractivePotential<double>> read_isolf_potential(const toml::value&);
template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, iSoLFAttractivePotential<float >> read_isolf_potential(const toml::value&);
}
