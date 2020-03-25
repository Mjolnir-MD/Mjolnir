#include <mjolnir/core/SystemMotionRemover.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class SystemMotionRemover<SimulatorTraits<double, UnlimitedBoundary>       >;
template class SystemMotionRemover<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class SystemMotionRemover<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class SystemMotionRemover<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
