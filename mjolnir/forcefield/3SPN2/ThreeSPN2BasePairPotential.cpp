#include <mjolnir/forcefield/3SPN2/ThreeSPN2BasePairPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2BasePairParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2BasePairParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2BasePairParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2BasePairParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
