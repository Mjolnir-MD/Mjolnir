#include <mjolnir/core/LoaderBase.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LoaderBase<SimulatorTraits<double, UnlimitedBoundary>>;
template class LoaderBase<SimulatorTraits<float,  UnlimitedBoundary>>;
template class LoaderBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class LoaderBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
