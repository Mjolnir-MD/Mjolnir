#include <mjolnir/core/XYZLoader.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class XYZLoader<SimulatorTraits<double, UnlimitedBoundary>>;
template class XYZLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
template class XYZLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class XYZLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
