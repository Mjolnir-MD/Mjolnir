#include <mjolnir/core/ExclusionList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ExclusionList<SimulatorTraits<double, UnlimitedBoundary>>;
template class ExclusionList<SimulatorTraits<float,  UnlimitedBoundary>>;
template class ExclusionList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ExclusionList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
