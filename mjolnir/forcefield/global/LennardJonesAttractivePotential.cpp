#include <mjolnir/forcefield/global/LennardJonesAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class LennardJonesAttractivePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class LennardJonesAttractivePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class LennardJonesAttractivePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class LennardJonesAttractivePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
