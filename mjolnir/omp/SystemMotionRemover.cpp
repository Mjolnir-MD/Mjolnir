#include <mjolnir/omp/SystemMotionRemover.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class SystemMotionRemover<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
template class SystemMotionRemover<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
template class SystemMotionRemover<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class SystemMotionRemover<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
