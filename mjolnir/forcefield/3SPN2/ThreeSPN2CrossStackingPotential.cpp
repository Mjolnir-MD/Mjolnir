#include <mjolnir/forcefield/3SPN2/ThreeSPN2CrossStackingPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
