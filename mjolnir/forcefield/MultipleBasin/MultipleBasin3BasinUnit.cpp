#include <mjolnir/forcefield/MultipleBasin/MultipleBasin3BasinUnit.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MultipleBasin3BasinUnit<SimulatorTraits<double, UnlimitedBoundary       >>;
template class MultipleBasin3BasinUnit<SimulatorTraits<float,  UnlimitedBoundary       >>;
template class MultipleBasin3BasinUnit<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MultipleBasin3BasinUnit<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}
