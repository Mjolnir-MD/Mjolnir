#include <mjolnir/input/read_random_number_generator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>       > read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>       > read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_rng(const toml::value&);
template RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_rng(const toml::value&);
} // mjolnir
