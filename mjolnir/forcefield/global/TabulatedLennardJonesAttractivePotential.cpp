#include <mjolnir/forcefield/global/TabulatedLennardJonesAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
