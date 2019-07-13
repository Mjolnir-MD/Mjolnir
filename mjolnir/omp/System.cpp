#include <mjolnir/omp/System.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class System<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class System<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class System<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class System<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
