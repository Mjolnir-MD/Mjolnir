#include <mjolnir/interaction/local/DummyInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class DummyInteraction<SimulatorTraits<double, UnlimitedBoundary>>;
template class DummyInteraction<SimulatorTraits<float,  UnlimitedBoundary>>;
template class DummyInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class DummyInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
