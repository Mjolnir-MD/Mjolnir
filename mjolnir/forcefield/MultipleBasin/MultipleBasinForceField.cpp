#include <mjolnir/forcefield/MultipleBasin/MultipleBasinForceField.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MultipleBasinForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
template class MultipleBasinForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
template class MultipleBasinForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class MultipleBasinForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
}
