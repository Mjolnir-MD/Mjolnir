#include <mjolnir/input/read_global_interaction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_3spn2_base_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_3spn2_base_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_3spn2_base_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_3spn2_base_pair_interaction(const toml::value&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_3spn2_cross_stacking_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_3spn2_cross_stacking_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_3spn2_cross_stacking_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_3spn2_cross_stacking_interaction(const toml::value&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value&);

template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_stoichiometric_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float , UnlimitedBoundary>       >> read_stoichiometric_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_stoichiometric_interaction(const toml::value&);
template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float , CuboidalPeriodicBoundary>>> read_stoichiometric_interaction(const toml::value&);
} // mjolnir
