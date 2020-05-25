#include <mjolnir/input/read_external_interaction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_interaction(const toml::value& external);

template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_position_restraint_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_position_restraint_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_position_restraint_interaction(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_position_restraint_interaction(const toml::value& external);

template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_recutangular_box_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_recutangular_box_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_recutangular_box_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_recutangular_box_interaction(const toml::value&);

template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_afm_flexible_fitting_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_afm_flexible_fitting_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_afm_flexible_fitting_interaction(const toml::value&);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_afm_flexible_fitting_interaction(const toml::value&);

template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_distance_interaction_shape(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_distance_interaction_shape(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_distance_interaction_shape(const toml::value& external);
template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_distance_interaction_shape(const toml::value& external);
} // mjolnir
