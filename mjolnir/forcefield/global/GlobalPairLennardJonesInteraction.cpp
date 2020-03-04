#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// L-J
template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

} // mjolnir
