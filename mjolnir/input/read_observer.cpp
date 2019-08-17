#include <mjolnir/input/read_observer.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template void add_observer<SimulatorTraits<double, UnlimitedBoundary>       >(ObserverContainer<SimulatorTraits<double, UnlimitedBoundary>       >& observers, const toml::value& format, const std::string& file_prefix);
template void add_observer<SimulatorTraits<float,  UnlimitedBoundary>       >(ObserverContainer<SimulatorTraits<float,  UnlimitedBoundary>       >& observers, const toml::value& format, const std::string& file_prefix);
template void add_observer<SimulatorTraits<double, CuboidalPeriodicBoundary>>(ObserverContainer<SimulatorTraits<double, CuboidalPeriodicBoundary>>& observers, const toml::value& format, const std::string& file_prefix);
template void add_observer<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(ObserverContainer<SimulatorTraits<float,  CuboidalPeriodicBoundary>>& observers, const toml::value& format, const std::string& file_prefix);

template ObserverContainer<SimulatorTraits<double, UnlimitedBoundary>       > read_observer<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& root);
template ObserverContainer<SimulatorTraits<float,  UnlimitedBoundary>       > read_observer<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& root);
template ObserverContainer<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_observer<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& root);
template ObserverContainer<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_observer<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& root);

} // mjolnir
