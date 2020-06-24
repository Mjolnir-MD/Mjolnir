#include <mjolnir/forcefield/MultipleBasin/MultipleBasin2BasinUnit.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MultipleBasin2BasinUnit<SimulatorTraits<double, UnlimitedBoundary       >>;
template class MultipleBasin2BasinUnit<SimulatorTraits<float,  UnlimitedBoundary       >>;
template class MultipleBasin2BasinUnit<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MultipleBasin2BasinUnit<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}
