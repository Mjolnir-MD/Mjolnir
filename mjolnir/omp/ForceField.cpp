#include <mjolnir/omp/ForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ForceField<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
template class ForceField<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
template class ForceField<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ForceField<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
