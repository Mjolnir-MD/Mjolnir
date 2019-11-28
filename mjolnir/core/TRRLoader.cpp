#include <mjolnir/core/TRRLoader.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class TRRLoader<SimulatorTraits<double, UnlimitedBoundary>>;
template class TRRLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
template class TRRLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class TRRLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
